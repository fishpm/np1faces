#### SUMMARY: Evaluation of whether treatment response is predicted by interaction between early life stress measure and amygdala reactivity (similar to Goldstein-Piekarski PNAS 2016, DOI: 10.1073/pnas.1606671113)

#### AUTHOR: Patrick M. Fisher, 2021

####
## RELEVANT PACKAGES
####

# prediction performance code
source('./fx_nruPredict.R')


####
## BASELINE ANALYSIS
####

# demographic data
d0 <- read.csv('./pubMaterials/data_demo_aal.csv', stringsAsFactors = F)

# define baseline data frame
d.baseline <- d0[d0$exclude.baseline=='',]
d.baseline$age0 <- d.baseline$age-mean(d.baseline$age)
d.baseline$remitter.status <- factor(d.baseline$remitter.status)

# emotion-specific contrast estimates
dd.aal <- read.csv('./pubMaterials/data_demo_aal_faces_order.csv',stringsAsFactors = F)

# define factors
dd.aal$emotion.factor <- factor(dd.aal$emotion.factor, levels = c('F','N','A','S'))
dd.aal$faces.block.factor <- factor(dd.aal$faces.block)
dd.aal$sex.factor <- factor(dd.aal$sex, levels = c('Female','Male'))
dd.aal$group.status <- factor(dd.aal$group.status, levels = c('Healthy Control', 'Case'))

####
### CATS
####

matches <- d.baseline$group.status=='Case' & d.baseline$exclude.treatresp==''

# aggregate data into relevant data frame
dd.pred <- do.call(rbind,
                   lapply(d.baseline[matches,'cimbi.id'], function(cimbi.id){
                       submatch <- which(d.baseline$cimbi.id == cimbi.id & matches)
                       cbind(cimbi.id,
                             d.baseline[submatch,c('age','sex','hamd6.baseline','rt.faces.all','rt.shapes.all','cats.mean','remitter.status')],
                             reshape(dd.aal[dd.aal$cimbi.id==cimbi.id,c('mr.id.full',
                                                                        'emotion.factor',
                                                                        'Amygdala_L',
                                                                        'Amygdala_R')],
                                     v.names = c('Amygdala_L', 'Amygdala_R'),
                                     timevar = 'emotion.factor',
                                     idvar = 'mr.id.full',
                                     direction = 'wide'))
                   }))

dd.pred$sex <- factor(dd.pred$sex)

# mean response across hemispheres
dd.pred$Amygdala.F <- (dd.pred$Amygdala_L.F+dd.pred$Amygdala_R.F)/2
dd.pred$Amygdala.N <- (dd.pred$Amygdala_L.N+dd.pred$Amygdala_R.N)/2
dd.pred$Amygdala.A <- (dd.pred$Amygdala_L.A+dd.pred$Amygdala_R.A)/2
dd.pred$Amygdala.S <- (dd.pred$Amygdala_L.S+dd.pred$Amygdala_R.S)/2

amy.terms <- colnames(dd.pred)[grepl('^(Amygdala.[FNAS]{1})', colnames(dd.pred))]
voi.terms <- c('cats.mean', amy.terms)

voi <- paste('cats.mean',amy.terms,sep='*')
covar <- c('age', 
           'sex', 
           'hamd6.baseline', 
           'rt.faces.all', 
           'rt.shapes.all')
y <- 'remitter.status' # binary outcome

# data frame to use contains only needed variables and complete cases
row.matches <- complete.cases(dd.pred[,c(covar,voi.terms,y)])
col.matches <- colnames(dd.pred)%in%c(covar,voi.terms,y)
dd.pred <- dd.pred[row.matches,col.matches]

nresample <- 1 # not needed with loocv
nperm <- 10000 # number of permutations

# model fit
modelObj <- fx_modelResample(df0 = dd.pred,              # data frame
                             cv.type = 'loocv',          # cross-validation type
                             covar = covar,              # covariate model
                             voi = voi,                  # variables of interest 
                             outcome = y,                # outcome
                             model.type = 'logistic',    # prediction model type
                             nresample = nresample,      # number of resamples
                             dthresh = 0.5,              # decision threshold (if needed)
                             z.pred = F,                 # standardize predictors
                             balance.col = y,            # columns to balance
                             n.cores = 10)               # parallel processing

# summarize performance metrics across resamples
modelPerfObj <- fx_modelResamplePerf(modelResampleObj = modelObj)

# run permutations
permObj <- fx_perm(df0 = dd.pred, modelObj = modelObj, nperm = nperm, n.cores = 6)

# summarize performance metrics across permutations and associated p-values
permPerfObj <- fx_permPerf(permObj = permObj, modelResamplePerf = modelPerfObj)

# save permutation results
outFile <- paste0('./pubMaterials/pred_',y,'_replication_logistic.RData')
save(modelObj,modelPerfObj,permPerfObj,file=outFile)

# read permutation results
writeLines(paste0('\nWorking on ', y))
fname <- paste0('./pubMaterials/pred_',y,'_replication_logistic.RData')
load(fname)

# organize performance metrics
txt1 <- paste0('ROC-AUC Covariate Model: ',
               signif(modelPerfObj$df.summary['avg','auc.ROC.covar'],3), ' \u00B1 ',
               signif(modelPerfObj$df.summary['std','auc.ROC.covar'],3), '; p = ',
               signif(permPerfObj$df.summary['pval','auc.ROC.covar'],2))

txt2 <- paste0('ROC-AUC Full Model: ',
               signif(modelPerfObj$df.summary['avg','auc.ROC.full'],3), ' \u00B1 ',
               signif(modelPerfObj$df.summary['std','auc.ROC.full'],3), '; p = ',
               signif(permPerfObj$df.summary['pval','auc.ROC.full'],2))

diff.mean <- mean(modelPerfObj$df.iter$auc.ROC.full-modelPerfObj$df.iter$auc.ROC.covar)
diff.sd <- sd(modelPerfObj$df.iter$auc.ROC.full-modelPerfObj$df.iter$auc.ROC.covar)
perm <- permPerfObj$df.iter$auc.ROC.full-permPerfObj$df.iter$auc.ROC.covar
diff.p <- sum(perm>diff.mean)/length(perm)

txt3 <- paste0('ROC-AUC Difference (Full-Covariate): ',
               signif(diff.mean,3), ' \u00B1 ',
               signif(diff.sd,3), '; p = ',
               signif(diff.p,2))

# write out performance metrics
writeLines(c(txt1,txt2,txt3))


## Logistic regression model without cross-validation

# four interaction terms in one model
f.covar <- paste0(y, '~', paste(covar,collapse='+'))
l.covar <- glm(f.covar, data = dd.pred, family = 'binomial')
summary(l.covar)

f.main <- paste0(y, '~', paste(c(covar,amy.terms),collapse='+'))
l.main <- glm(f.main, data = dd.pred, family = 'binomial')
summary(l.main)

f.full <- paste0(y, '~', paste(covar,collapse='+'), '+', paste(voi,collapse='+'))
l.full <- glm(f.full, data = dd.pred, family = 'binomial')
summary(l.full)

compare(l.covar,l.main)
compare(l.covar,l.full)
compare(l.main,l.full)

