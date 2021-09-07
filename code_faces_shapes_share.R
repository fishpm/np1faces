#### SUMMARY: Analyses of faces vs. shapes contrast

#### AUTHOR: Patrick M. Fisher, 2021

####
## RELEVANT PACKAGES
####

library(ggbeeswarm) # beeswarm plots
library(parallel) # mcapply


####
## READ DATA FRAME
####

# prepared csv file with relevant values (code_dataLoad.R)
d0 <- read.csv('./pubMaterials/data_demo_aal.csv', stringsAsFactors = F)

# region name text file
aal.names <- as.vector(unlist(read.delim('./pubMaterials/atlas116.txt',
                                         sep = '\t', stringsAsFactors = F)))

# rename two AAL region names redundant with other names
aal.names[25] <- 'Frontal_Med_Orb_L' # was Frontal_Mid_Orb_L, which duplicates 9
aal.names[26] <- 'Frontal_Med_Orb_R' # was Frontal_Mid_Orb_R, which duplicates 10

####
## BASELINE ANALYSIS
####

# define baseline data frame
d.baseline <- d0[d0$exclude.baseline=='',]
d.baseline$age0 <- d.baseline$age-mean(d.baseline$age)

####
## MAIN EFFECT OF TASK - AAL
####

nregions <- 90 # exclude cerebellar regions
covar <- paste(c('age', 'sex', 'rt.faces.all', 'rt.shapes.all'),collapse='+')
covar0 <- paste(c('age0', 'sex', 'rt.faces.all0', 'rt.shapes.all0'),collapse='+')

# data frame with estimates of interest
stats.out <- data.frame(t(sapply(seq(nregions), function(i){
    
    # model with covariates
    f0 <- paste0(aal.names[i], '~', covar)
    l0 <- lm(f0, d.baseline)
    
    # extract mean-centered estimates from data subset
    l0$model$rt.faces.all0 <- l0$model$rt.faces.all-mean(l0$model$rt.faces.all)
    l0$model$rt.shapes.all0 <- l0$model$rt.shapes.all-mean(l0$model$rt.shapes.all)
    l0$model$age0 <- l0$model$age-mean(l0$model$age)
    
    # model with mean-centered covariates
    f <- paste0(aal.names[i], '~', covar0)
    l <- lm(f,data = l0$model)
    
    # responses adjusted for covariates
    res <- l$residuals + coef(l)['(Intercept)'] 
    
    # cohen's d
    coh.d <- coef(l)['(Intercept)']/sd(res) 
    
    return(c(coef(l)['(Intercept)'],
             confint(l)['(Intercept)',], 
             summary(l)$coefficients['(Intercept)','t value'],
             summary(l)$coefficients['(Intercept)','Pr(>|t|)'],
             coh.d))
    
})))

# column names
colnames(stats.out) <- c('est','lwr','upr','tstat','p','d')

# region names
stats.out$regions <- aal.names[seq(nregions)]

# bonferroni-holm adjustment
stats.out$padj <- p.adjust(stats.out$p, method='holm')

# minus log10 p (adjusted values)
stats.out$log10p <- -log10(stats.out$padj)

# color marker (useful for plots)
stats.out$color <- factor(unlist(lapply(seq(nrow(stats.out)), function(i){
    if(stats.out[i,'padj']<0.05){return('#06E406')}
    if(stats.out[i,'p']<0.05){return('#E49E09')}
    return('#000000')
})), levels = c('#06E406','#E49E09','#000000'))

####
## Data are foundation for Supplementary Table 1
####

write.csv(stats.out,'./pubMaterials/aal_maineffects.csv',row.names=F)


####
## PREDICTION OF TREATMENT RESPONSE FROM BASELINE FMRI MEASURES
####

# classification groupings
class.groups <- c('np1.binaryresponse','remitter.status','responder.status')
class.names <- c('np1','h17_rem','h17_res')

# cortical/subcortical regions
nregions <- 90

# number of resamples (to estimate variance of performance metrics)
nresample <- 100

# number of permutations (to estimate statistical significance)
nperm <- 10000

# for each classification
for(i in seq(length(class.groups))){
    
    writeLines(paste0('Working on ', class.groups[i]))
    
    
    # set of regional predictors
    voi <- aal.names[seq(nregions)]
    
    # covariates
    covar <- c('age', 'sex', 'hamd6.baseline', 'rt.faces.all', 'rt.shapes.all')
    
    # binary outcome
    y <- class.groups[i]
    
    # get appropriate subset of data
    matches <- d.baseline$exclude.treatresp=='' & d.baseline$group.status=='Case'
    dd <- d.baseline[matches,]
    row.matches <- complete.cases(dd[,c(covar,voi,y)])
    col.matches <- colnames(dd)%in%c(covar,voi,y)
    dd <- dd[row.matches,col.matches]
    
    # set as factor
    dd$sex <- as.factor(dd$sex)
    
    # fit classification model
    modelObj <- fx_modelResample(df0 = dd, # data frame
                                 cv.type = '5-fold', # type of cross-validation
                                 covar = covar, # covariate set
                                 voi = voi,  # variables of interest (i.e., brain regions)
                                 outcome = y, # class
                                 model.type = 'rf', # model type (randomForest)
                                 nresample = nresample, # number of resamples
                                 dthresh = 0.5, # threshold (not used)
                                 z.pred = F, # standardize continuous predictors
                                 balance.col = class.groups[i], # stratified cv
                                 n.cores = 10) # parallel processing
    
    # determine overall model performance
    modelPerfObj <- fx_modelResamplePerf(modelResampleObj = modelObj)
    
    # permutation testing
    permObj <- fx_perm(df0 = dd, modelObj = modelObj, nperm = nperm, n.cores = 10)
    
    # determine permutation test performance
    permPerfObj <- fx_permPerf(permObj = permObj, modelResamplePerf = modelPerfObj)
    
    # save relevant objects
    outFile <- paste0(top,'pubMaterials/pred_',
                      class.names[i],'.RData')
    save(modelObj,modelPerfObj,permPerfObj,file=outFile)
    
}

# summary outputs from model performance

# classification groupings
class.groups <- c('np1.binaryresponse','remitter.status','responder.status')
class.names <- c('np1','h17_rem','h17_res')

# for each classification
for(i in seq(length(class.groups))){
    
    writeLines(paste0('\nWorking on ', class.groups[i]))
    
    # load object with previuosly generated model performance data
    fname <- paste0(top,'pubMaterials/pred_',
                    class.names[i],'.RData')
    load(fname)
    
    # print covariate model performance
    txt1 <- paste0('ROC-AUC Covariate Model: ',
                   signif(modelPerfObj$df.summary['avg','auc.ROC.covar'],3), ' \u00B1 ',
                   signif(modelPerfObj$df.summary['std','auc.ROC.covar'],3), '; p = ',
                   signif(permPerfObj$df.summary['pval','auc.ROC.covar'],2))
    
    # print full model performance (covaritates + vois)
    txt2 <- paste0('ROC-AUC Full Model: ',
                   signif(modelPerfObj$df.summary['avg','auc.ROC.full'],3), ' \u00B1 ',
                   signif(modelPerfObj$df.summary['std','auc.ROC.full'],3), '; p = ',
                   signif(permPerfObj$df.summary['pval','auc.ROC.full'],2))
    
    # compute effect of vois on model performance
    diff.mean <- mean(modelPerfObj$df.iter$auc.ROC.full-modelPerfObj$df.iter$auc.ROC.covar)
    diff.sd <- sd(modelPerfObj$df.iter$auc.ROC.full-modelPerfObj$df.iter$auc.ROC.covar)
    perm <- permPerfObj$df.iter$auc.ROC.full-permPerfObj$df.iter$auc.ROC.covar
    
    # p-value is one sided test (count only when permutation beats observed)
    diff.p <- sum(perm>diff.mean)/length(perm)
    
    # print differential model performance
    txt3 <- paste0('ROC-AUC Difference (Full-Covariate): ',
                   signif(diff.mean,3), ' \u00B1 ',
                   signif(diff.sd,3), '; p = ',
                   signif(diff.p,2))
    
    writeLines(c(txt1,txt2,txt3))
    
    # generate roc curve (unpolished)
    fx_rocPlot(modelObj = modelObj,
               permPerfObj = permPerfObj,
               modelResamplePerf = modelPerfObj,
               title.name = 'ROC Curve',
               plot.covar = T,
               plot.full = T)
    
}


