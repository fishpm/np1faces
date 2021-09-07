#### SUMMARY: Publication figure code

#### AUTHOR: Patrick M. Fisher, 2021

####
## RELEVANT PACKAGES
####

library(ggplot2)

####
## FIGURE 1
####

# load file
d0 <- read.csv('./pubMaterials/data_demo_aal.csv', stringsAsFactors = F)

# define baseline data frame
d.baseline <- d0[d0$exclude.baseline=='',]
d.baseline$age0 <- d.baseline$age-mean(d.baseline$age)

# covariate sets
covar <- paste(c('age', 'sex', 'rt.faces.all', 'rt.shapes.all'),collapse='+')
covar0 <- paste(c('age0', 'sex', 'rt.faces.all0', 'rt.shapes.all0'),collapse='+')

# define model to extract info for plots
f <- paste0(aal.names[41], '~', covar)
l <- lm(f,data = d.baseline)
dd.plot <- data.frame(cimbi.id = rep(d.baseline[rownames(l$model),'cimbi.id'],2),
                      group = rep(d.baseline[rownames(l$model),'group.status'],2)
)

# adjusted residuals
dd.plot$res <- as.vector(sapply(c('Amygdala_L', 'Amygdala_R'), function(i){
    
    # model with observed covariates
    f0 <- paste0(aal.names[i], '~', covar)
    l0 <- lm(f0, d.baseline)
    
    # mean-centered values for model-specific data
    l0$model$rt.faces.all0 <- l0$model$rt.faces.all-mean(l0$model$rt.faces.all)
    l0$model$rt.shapes.all0 <- l0$model$rt.shapes.all-mean(l0$model$rt.shapes.all)
    l0$model$age0 <- l0$model$age-mean(l0$model$age)
    
    # model with mean-centered covariates
    f <- paste0(aal.names[i], '~', covar0)
    l <- lm(f,data = l0$model)
    
    # response value, adjusted for covariates
    return(l$residuals + coef(l)['(Intercept)'])}))

# region column
dd.plot$region <- rep(c('Left Amygdala','Right Amygdala'),each=nrow(dd.plot)/2)

# amygdala plots
ggplot(dd.plot, aes(x = region, y = res, group = group, color = group)) + 
    geom_beeswarm(size = 4, alpha = 0.5, cex = 1.5, dodge.width = 0.75) +
    stat_summary(fun.data = 'mean_sdl',
                 fun.args = list(mult=1),
                 position = position_dodge(0.75),
                 geom = 'errorbar',
                 size = 3,
                 show.legend = F,
                 width = 0.25) +
    stat_summary(fun = 'mean',
                 cex = 15,
                 position = position_dodge(0.75),
                 geom = 'point',
                 show.legend = F) +
    labs(x = 'Region',
         y = 'Amygdala response (arbitrary units, adjusted)',
         color = 'Group') + 
    scale_color_discrete(labels = c('Healthy','Depressed')) +
    scale_y_continuous(limits = c(-4.1,4.1), breaks = seq(-4,4)) +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          legend.title.align = 0.5)

#### this plot joined with SPM output in inkscape


####
## FIGURE 2
####

# see FÌGURE 2 sections in code_faces_order_share.R (2 parts)

####
## FIGURE 3; SUPPLEMENTARY FIGURE 2; SUPPLEMENTARY FIGURE 3
####

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
    
    # assign classes
    negClass <- modelObj[[2]]$class.levels[1]
    posClass <- modelObj[[2]]$class.levels[2]
    
    # steps
    thresh.set <- seq(0,1,by=0.05)
    
    # data framce for roc curve generation
    dd.cv <- data.frame(resample = rep(seq(length(modelObj[[1]])), each = length(thresh.set)),
                        threshold = rep(thresh.set, length(modelObj[[1]])))
    
    # compute tpr and fpr at all steps for covariate model
    covar.out <- lapply(seq(length(modelObj[[1]])), function(nfold){
        out <- t(sapply(thresh.set, function(currThresh){
            if(currThresh==0){
                tpr <- 1
                fpr <- 1
            } else if(currThresh==1){
                tpr <- 0
                fpr <- 0
            } else {
                covar.prob <- modelObj[[1]][[nfold]]$df.allfolds$pred.prob.covar
                covar.bin <- as.numeric(covar.prob > currThresh) + 1
                if (length(unique(covar.bin))==1){
                    if(all(covar.bin==1)){
                        tpr <- 0
                        fpr <- 0
                    } else if(all(covar.bin==2)){
                        tpr <- 1
                        fpr <- 1
                    }
                } else {
                    currGuess <- factor(covar.bin, labels = modelObj[[2]]$class.levels)
                    actualClass <- modelObj[[1]][[nfold]]$df.allfolds$actual.class
                    tp <- sum(currGuess == posClass & actualClass == posClass) # true positive
                    tn <- sum(currGuess == negClass & actualClass == negClass) # true negative
                    fp <- sum(currGuess == posClass & actualClass == negClass) # false positive
                    fn <- sum(currGuess == negClass & actualClass == posClass) # false negative
                    tpr <- tp/(tp+fn) # true positive rate
                    tnr <- tn/(tn+fp) # true negative rate
                    fpr <- 1 - tnr # fpr == 1 - tnr
                }
            }
            return(c(tpr,fpr))
        }))
        colnames(out) <- c('tpr','fpr')
        return(list(tpr=out[,'tpr'],fpr=out[,'fpr']))
    })
    
    # compute tpr and fpr at all steps for covariate model
    full.out <- lapply(seq(length(modelObj[[1]])), function(nfold){
        out <- t(sapply(thresh.set, function(currThresh){
            if(currThresh==0){
                tpr <- 1
                fpr <- 1
            } else if(currThresh==1){
                tpr <- 0
                fpr <- 0
            } else {
                full.prob <- modelObj[[1]][[nfold]]$df.allfolds$pred.prob.full
                full.bin <- as.numeric(full.prob > currThresh) + 1
                if (length(unique(full.bin))==1){
                    if(all(full.bin==1)){
                        tpr <- 0
                        fpr <- 0
                    } else if(all(full.bin==2)){
                        tpr <- 1
                        fpr <- 1
                    }
                } else {
                    currGuess <- factor(full.bin, labels = modelObj[[2]]$class.levels)
                    actualClass <- modelObj[[1]][[nfold]]$df.allfolds$actual.class
                    tp <- sum(currGuess == posClass & actualClass == posClass) # true positive
                    tn <- sum(currGuess == negClass & actualClass == negClass) # true negative
                    fp <- sum(currGuess == posClass & actualClass == negClass) # false positive
                    fn <- sum(currGuess == negClass & actualClass == posClass) # false negative
                    tpr <- tp/(tp+fn) # true positive rate
                    tnr <- tn/(tn+fp) # true negative rate
                    fpr <- 1 - tnr # fpr == 1 - tnr
                }
            }
            return(c(tpr,fpr))
        }))
        colnames(out) <- c('tpr','fpr')
        return(list(tpr=out[,'tpr'],fpr=out[,'fpr']))
    })
    
    # assign values to dd.cv
    dd.cv$tpr.covar <- as.vector(sapply(seq(length(covar.out)), function(nresample){
        covar.out[[nresample]]$tpr
    }))
    dd.cv$fpr.covar <- as.vector(sapply(seq(length(covar.out)), function(nresample){
        covar.out[[nresample]]$fpr
    }))
    dd.cv$tpr.full <- as.vector(sapply(seq(length(full.out)), function(nresample){
        full.out[[nresample]]$tpr
    }))
    dd.cv$fpr.full <- as.vector(sapply(seq(length(full.out)), function(nresample){
        full.out[[nresample]]$fpr
    }))
    
    # summary values for plots
    dd.cv.summary <- data.frame(threshold = thresh.set)
    summary.vals <- t(sapply(thresh.set, function(j){
        
        covar.tpr.mean <- mean(dd.cv$tpr.covar[dd.cv$threshold==j])
        full.tpr.mean <- mean(dd.cv$tpr.full[dd.cv$threshold==j])
        
        covar.tpr.lwr <- covar.tpr.mean - sd(dd.cv$tpr.covar[dd.cv$threshold==j])
        covar.tpr.upr <- covar.tpr.mean + sd(dd.cv$tpr.covar[dd.cv$threshold==j])
        full.tpr.lwr <- full.tpr.mean - sd(dd.cv$tpr.full[dd.cv$threshold==j])
        full.tpr.upr <- full.tpr.mean + sd(dd.cv$tpr.full[dd.cv$threshold==j])
        
        covar.fpr.mean <- mean(dd.cv$fpr.covar[dd.cv$threshold==j])
        full.fpr.mean <- mean(dd.cv$fpr.full[dd.cv$threshold==j])
        
        covar.fpr.lwr <- covar.fpr.mean - sd(dd.cv$fpr.covar[dd.cv$threshold==j])
        covar.fpr.upr <- covar.fpr.mean + sd(dd.cv$fpr.covar[dd.cv$threshold==j])
        full.fpr.lwr <- full.fpr.mean - sd(dd.cv$fpr.full[dd.cv$threshold==j])
        full.fpr.upr <- full.fpr.mean + sd(dd.cv$fpr.full[dd.cv$threshold==j])
        
        return(c(covar.tpr.mean,
                 covar.tpr.lwr,
                 covar.tpr.upr,
                 full.tpr.mean,
                 full.tpr.lwr,
                 full.tpr.upr,
                 covar.fpr.mean,
                 covar.fpr.lwr,
                 covar.fpr.upr,
                 full.fpr.mean,
                 full.fpr.lwr,
                 full.fpr.upr))
    }))
    colnames(summary.vals) <- c('covar.tpr.mean',
                                'covar.tpr.lwr',
                                'covar.tpr.upr',
                                'full.tpr.mean',
                                'full.tpr.lwr',
                                'full.tpr.upr',
                                'covar.fpr.mean',
                                'covar.fpr.lwr',
                                'covar.fpr.upr',
                                'full.fpr.mean',
                                'full.fpr.lwr',
                                'full.fpr.upr')
    
    summary.vals <- as.data.frame(summary.vals)
    dd.cv.summary[,c('covar.tpr.mean',
                     'covar.tpr.lwr',
                     'covar.tpr.upr',
                     'full.tpr.mean',
                     'full.tpr.lwr',
                     'full.tpr.upr',
                     'covar.fpr.mean',
                     'covar.fpr.lwr',
                     'covar.fpr.upr',
                     'full.fpr.mean',
                     'full.fpr.lwr',
                     'full.fpr.upr')] <- summary.vals
    
    # Generate ROC curves
    # FIGURE 3; SUPPLEMENTARY FIGURE 2; SUPPLEMENTARY FIGURE 3
    print(ggplot(dd.cv.summary, aes(x=covar.fpr.mean,y=covar.tpr.mean)) +
              geom_segment(aes(x=0, y=0, xend=1, yend=1), color = 'black') +
              geom_line(color = 'red', lwd = 2) +
              geom_line(aes(x=full.fpr.mean,y=full.tpr.mean),color = 'blue', lwd = 2) +
              geom_ribbon(aes(x=covar.fpr.mean,
                              ymin=covar.tpr.lwr,
                              ymax=covar.tpr.upr),
                          alpha = 0.2,
                          fill = 'red') +
              geom_ribbon(aes(x=full.fpr.mean,
                              ymin=full.tpr.lwr,
                              ymax=full.tpr.upr), 
                          alpha = 0.2,
                          fill = 'blue') +
              labs(title = class.groups[i],
                   x = '1-Specificity',
                   y = 'Sensitivity',
                   color = 'legend') +
              theme(plot.title = element_text(hjust = 0.5),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    axis.text.x = element_text(size = 16),
                    axis.text.y = element_text(size = 16)))
    
}

####
## SUPPLEMENTARY FIGURE 1
####

library(ggbeeswarm) # beeswarm plots

# region name text file
aal.names <- as.vector(unlist(read.delim('./pubMaterials/atlas116.txt',
                                         sep = '\t', stringsAsFactors = F)))

# rename two AAL region names redundant with other names
aal.names[25] <- 'Frontal_Med_Orb_L' # was Frontal_Mid_Orb_L, which duplicates 9
aal.names[26] <- 'Frontal_Med_Orb_R' # was Frontal_Mid_Orb_R, which duplicates 10

d.treatment <- read.csv('./pubMaterials/data_demo_aal_treatment.csv', stringsAsFactors = F)


nregions <- 90

idx <- data.frame(t(sapply(unique(d.treatment$cimbi.id), function(i){
    return(c(which(d.treatment$cimbi.id==i & d.treatment$follow.up==F),
             which(d.treatment$cimbi.id==i & d.treatment$follow.up==T)))
})))
colnames(idx) <- c('baseline','rescan')
rownames(idx) <- unique(d.treatment$cimbi.id)

stats.out2 <- data.frame(t(sapply(seq(nregions), function(i){
    
    vals.base <- d.treatment[idx$baseline,aal.names[i]]
    vals.rescan <- d.treatment[idx$rescan,aal.names[i]]
    vals.diff <- vals.rescan-vals.base
    ttest <- t.test(vals.diff)
    coh.d <- mean(vals.diff)/sd(vals.diff)
    
    return(c(ttest$estimate,
             ttest$conf.int,
             coh.d,
             ttest$p.value,
             mean(vals.base),
             sd(vals.base),
             mean(vals.rescan),
             sd(vals.rescan),
             mean(vals.diff),
             sd(vals.diff)))
})))

colnames(stats.out2) <- c('est','lwr','upr','d', 'p', 'base.mean', 'base.sd', 'rescan.mean', 'rescan.sd', 'diff.mean', 'diff.sd')
stats.out2$regions <- aal.names[seq(nregions)]
stats.out2$p.adj <- p.adjust(stats.out2$p, method='holm')
stats.out2$color <- factor(unlist(lapply(seq(nrow(stats.out2)), function(i){
    if(stats.out2[i,'p.adj']<0.05){return('#06E406')}
    if(stats.out2[i,'p']<0.05){return('#E49E09')}
    return('#000000')
})), levels = c('#06E406','#E49E09','#000000'))

ggplot(stats.out2, aes(x=base.mean,y=rescan.mean,color=color)) +
    geom_abline(slope = 1, intercept = 0, lwd = 2, color = 'gray', lty = 2) +
    geom_point(show.legend = F, size = 3) +
    scale_y_continuous(limits = c(-0.6, 2.1), 
                       breaks = seq(-0.6, 2.1, by = 0.3)) +
    scale_x_continuous(limits = c(-0.6, 2.1),
                       breaks = seq(-0.6, 2.1, by = 0.3)) +
    scale_color_manual(breaks = c('#06E406','#E49E09','#000000'),
                       values = c('#06E406','#E49E09','#000000')) + 
    labs(x = 'Baseline brain responses (arbitrary units)',
         y = 'Week 8 brain responses\n(arbitrary units)') +
    theme(axis.text.y = element_text(size=16),
          axis.text.x = element_text(size=16),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20))

####
## SUPPLEMENTARY FIGURE 4
####

writeLines('\nWorking on remitter.status')
fname <- paste0('./pubMaterials/pred_remitter.status_replication_logistic.RData')

load(fname)

# define classes
negClass <- modelObj[[2]]$class.levels[1]
posClass <- modelObj[[2]]$class.levels[2]

thresh.set <- seq(0,1,by=0.05)
dd.cv <- data.frame(resample = rep(seq(length(modelObj[[1]])), each = length(thresh.set)),
                    threshold = rep(thresh.set, length(modelObj[[1]])))

covar.out <- lapply(seq(length(modelObj[[1]])), function(nfold){
    out <- t(sapply(thresh.set, function(currThresh){
        if(currThresh==0){
            tpr <- 1
            fpr <- 1
        } else if(currThresh==1){
            tpr <- 0
            fpr <- 0
        } else {
            covar.prob <- modelObj[[1]][[nfold]]$df.allfolds$pred.prob.covar
            covar.bin <- as.numeric(covar.prob > currThresh) + 1
            if (length(unique(covar.bin))==1){
                if(all(covar.bin==1)){
                    tpr <- 0
                    fpr <- 0
                } else if(all(covar.bin==2)){
                    tpr <- 1
                    fpr <- 1
                }
            } else {
                currGuess <- factor(covar.bin, labels = modelObj[[2]]$class.levels)
                actualClass <- modelObj[[1]][[nfold]]$df.allfolds$actual.class
                tp <- sum(currGuess == posClass & actualClass == posClass) # true positive
                tn <- sum(currGuess == negClass & actualClass == negClass) # true negative
                fp <- sum(currGuess == posClass & actualClass == negClass) # false positive
                fn <- sum(currGuess == negClass & actualClass == posClass) # false negative
                tpr <- tp/(tp+fn) # true positive rate
                tnr <- tn/(tn+fp) # true negative rate
                fpr <- 1 - tnr # fpr == 1 - tnr
            }
        }
        return(c(tpr,fpr))
    }))
    colnames(out) <- c('tpr','fpr')
    return(list(tpr=out[,'tpr'],fpr=out[,'fpr']))
})

full.out <- lapply(seq(length(modelObj[[1]])), function(nfold){
    out <- t(sapply(thresh.set, function(currThresh){
        if(currThresh==0){
            tpr <- 1
            fpr <- 1
        } else if(currThresh==1){
            tpr <- 0
            fpr <- 0
        } else {
            full.prob <- modelObj[[1]][[nfold]]$df.allfolds$pred.prob.full
            full.bin <- as.numeric(full.prob > currThresh) + 1
            if (length(unique(full.bin))==1){
                if(all(full.bin==1)){
                    tpr <- 0
                    fpr <- 0
                } else if(all(full.bin==2)){
                    tpr <- 1
                    fpr <- 1
                }
            } else {
                currGuess <- factor(full.bin, labels = modelObj[[2]]$class.levels)
                actualClass <- modelObj[[1]][[nfold]]$df.allfolds$actual.class
                tp <- sum(currGuess == posClass & actualClass == posClass) # true positive
                tn <- sum(currGuess == negClass & actualClass == negClass) # true negative
                fp <- sum(currGuess == posClass & actualClass == negClass) # false positive
                fn <- sum(currGuess == negClass & actualClass == posClass) # false negative
                tpr <- tp/(tp+fn) # true positive rate
                tnr <- tn/(tn+fp) # true negative rate
                fpr <- 1 - tnr # fpr == 1 - tnr
            }
        }
        return(c(tpr,fpr))
    }))
    colnames(out) <- c('tpr','fpr')
    return(list(tpr=out[,'tpr'],fpr=out[,'fpr']))
})

dd.cv$tpr.covar <- as.vector(sapply(seq(length(covar.out)), function(nresample){
    covar.out[[nresample]]$tpr
}))
dd.cv$fpr.covar <- as.vector(sapply(seq(length(covar.out)), function(nresample){
    covar.out[[nresample]]$fpr
}))

dd.cv$tpr.full <- as.vector(sapply(seq(length(full.out)), function(nresample){
    full.out[[nresample]]$tpr
}))
dd.cv$fpr.full <- as.vector(sapply(seq(length(full.out)), function(nresample){
    full.out[[nresample]]$fpr
}))

dd.cv.summary <- data.frame(threshold = thresh.set)

summary.vals <- t(sapply(thresh.set, function(j){
    
    covar.tpr.mean <- mean(dd.cv$tpr.covar[dd.cv$threshold==j])
    full.tpr.mean <- mean(dd.cv$tpr.full[dd.cv$threshold==j])
    
    covar.tpr.lwr <- covar.tpr.mean - sd(dd.cv$tpr.covar[dd.cv$threshold==j])
    covar.tpr.upr <- covar.tpr.mean + sd(dd.cv$tpr.covar[dd.cv$threshold==j])
    full.tpr.lwr <- full.tpr.mean - sd(dd.cv$tpr.full[dd.cv$threshold==j])
    full.tpr.upr <- full.tpr.mean + sd(dd.cv$tpr.full[dd.cv$threshold==j])
    
    covar.fpr.mean <- mean(dd.cv$fpr.covar[dd.cv$threshold==j])
    full.fpr.mean <- mean(dd.cv$fpr.full[dd.cv$threshold==j])
    
    covar.fpr.lwr <- covar.fpr.mean - sd(dd.cv$fpr.covar[dd.cv$threshold==j])
    covar.fpr.upr <- covar.fpr.mean + sd(dd.cv$fpr.covar[dd.cv$threshold==j])
    full.fpr.lwr <- full.fpr.mean - sd(dd.cv$fpr.full[dd.cv$threshold==j])
    full.fpr.upr <- full.fpr.mean + sd(dd.cv$fpr.full[dd.cv$threshold==j])
    
    return(c(covar.tpr.mean,
             covar.tpr.lwr,
             covar.tpr.upr,
             full.tpr.mean,
             full.tpr.lwr,
             full.tpr.upr,
             covar.fpr.mean,
             covar.fpr.lwr,
             covar.fpr.upr,
             full.fpr.mean,
             full.fpr.lwr,
             full.fpr.upr))
}))
colnames(summary.vals) <- c('covar.tpr.mean',
                            'covar.tpr.lwr',
                            'covar.tpr.upr',
                            'full.tpr.mean',
                            'full.tpr.lwr',
                            'full.tpr.upr',
                            'covar.fpr.mean',
                            'covar.fpr.lwr',
                            'covar.fpr.upr',
                            'full.fpr.mean',
                            'full.fpr.lwr',
                            'full.fpr.upr')
summary.vals <- as.data.frame(summary.vals)
dd.cv.summary[,c('covar.tpr.mean',
                 'covar.tpr.lwr',
                 'covar.tpr.upr',
                 'full.tpr.mean',
                 'full.tpr.lwr',
                 'full.tpr.upr',
                 'covar.fpr.mean',
                 'covar.fpr.lwr',
                 'covar.fpr.upr',
                 'full.fpr.mean',
                 'full.fpr.lwr',
                 'full.fpr.upr')] <- summary.vals

print(ggplot(dd.cv.summary, aes(x=covar.fpr.mean,y=covar.tpr.mean)) +
          geom_segment(aes(x=0, y=0, xend=1, yend=1), color = 'black') +
          geom_line(color = 'red', lwd = 2) +
          geom_line(aes(x=full.fpr.mean,y=full.tpr.mean),color = 'blue', lwd = 2) +
          geom_ribbon(aes(x=covar.fpr.mean,
                          ymin=covar.tpr.lwr,
                          ymax=covar.tpr.upr),
                      alpha = 0.2,
                      fill = 'red') +
          geom_ribbon(aes(x=full.fpr.mean,
                          ymin=full.tpr.lwr,
                          ymax=full.tpr.upr), 
                      alpha = 0.2,
                      fill = 'blue') +
          labs(title = y,
               x = '1-Specificity',
               y = 'Sensitivity',
               color = 'legend') +
          theme(plot.title = element_text(hjust = 0.5),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16)))
