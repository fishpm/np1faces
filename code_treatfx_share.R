#### SUMMARY: Evaluation of change in fMRI responses following antidepressant treatment.

#### AUTHOR: Patrick M. Fisher, 2021

####
## RELEVANT PACKAGES
####

library(ggbeeswarm) # beeswarm plots

d.treatment <- read.csv('./pubMaterials/data_demo_aal_treatment.csv', stringsAsFactors = F)


writeLines('**** GROUP SAMPLE SIZES ****')
table(d.treatment$sex[d.treatment$follow.up==F])

# variables of interest
desc.cont <- c('age',                 # age in years
               'rt.faces.all',        # reaction time to faces (fmri task)
               'rt.shapes.all',       # reaction time to shapes (fmri task)
               'acc.faces',           # accuracy faces (fmri task)
               'acc.shapes',          # accuracy shapes (fmri task)
               'outlier.total',       # outlier volumes (fmri)
               'outlier.faces',       # outlier faces volumes (fmri)
               'outlier.shapes',      # outlier shapes volumes (fmri)
               'hamd6.baseline',      # hamd-6 at baseline
               'hamd6.w8',            # hamd-6 at week 8
               'hamd17.baseline',     # hamd-17 at baseline
               'hamd17.w8',           # hamd-17 at week 8
               'srt.rt',              # simple reaction time (ms, neuropsych)
               'srt.std',             # simple reaction time st dev (ms, neuropsych)
               'srt.acc',             # simple reaction time (accuracy, neuropsych)
               'cats.mean',           # CATS (mean)
               'cats.sum'             # CATS (total)
)

# clinical measures of interest
clin.measures <- c('hamd6.baseline', 'hamd6.w8', 'hamd17.baseline', 'hamd17.w8')

groups <- c('Baseline','Follow-Up')

####
## Data are foundation for SUPPLEMENTARY TABLE 3
####

writeLines('**** GROUP DESCRIPTIVE INFORMATION ****')
for(g in groups){
    
    if(g=='Baseline'){
        row.set <- which(d.treatment$follow.up==F)
    } else if(g=='Follow-Up'){
        row.set <- which(d.treatment$follow.up==T)
    }
    
    txt <- paste(unlist(lapply(desc.cont, function(i){
        nomit <- sum(is.na(d.treatment[row.set,i]))
        return(paste0(i, ': ', 
                      signif(mean(d.treatment[row.set,i], na.rm=T),4), ' \u00B1 ',
                      signif(sd(d.treatment[row.set,i], na.rm=T),4), ' [',
                      signif(median(d.treatment[row.set,i], na.rm=T),4), '; ',
                      signif(min(d.treatment[row.set,i], na.rm=T),4), '-',
                      signif(max(d.treatment[row.set,i], na.rm=T),4), ']',
                      ' (nomit: ', nomit, ')\n\t'))
    })),collapse='')
    
    writeLines(paste0('**** ', g, ' ****\n\t',txt))
}

txt2 <- paste(unlist(lapply(desc.cont, function(i){
    if(i%in%clin.measures){
        
        if(i%in%c('hamd17.baseline','hamd6.baseline')){
            currMeasure <- unlist(strsplit(i,'.',fixed=T))[1]
            
            diff.score <- d.treatment[d.treatment$follow.up==F,paste0(currMeasure,'.baseline')]-
                d.treatment[d.treatment$follow.up==F,paste0(currMeasure,'.w8')]
            
            ttest <- t.test(diff.score)
            
            return(paste0('Delta-', currMeasure, ': ', signif(ttest$estimate,4), '[',
                          paste(signif(ttest$conf.int,3),collapse=', '), ']; p = ',
                          signif(ttest$p.value,3), '\n\t'))
        } else {return(NULL)}
        
    } else {
        
        baseline <- sapply(unique(d.treatment$cimbi.id), function(j){
            matches <- d.treatment$cimbi.id==j & d.treatment$follow.up==F
            return(d.treatment[matches,i])
        })
        rescan <- sapply(unique(d.treatment$cimbi.id), function(j){
            matches <- d.treatment$cimbi.id==j & d.treatment$follow.up==T
            return(d.treatment[matches,i])
        })
        ttest <- t.test(rescan,baseline, paired = T)
        
        return(paste0(i, ': ', signif(ttest$estimate,4), '[',
                      paste(signif(ttest$conf.int,4),collapse=', '), ']; p = ',
                      signif(ttest$p.value,4), '\n\t'))
    }})),collapse='')

writeLines(paste0('**** GROUP DIFFERENCES: DESCRIPTIVES (PARAMETRIC) ****\n\t', txt2))

txt3 <- paste(unlist(lapply(desc.cont, function(i){
    if(i%in%c('age')){
        return(NULL)
    }
    if(i%in%clin.measures){
        
        if(i%in%c('hamd17.baseline','hamd6.baseline')){
            currMeasure <- unlist(strsplit(i,'.',fixed=T))[1]
            
            diff.score <- d.treatment[d.treatment$follow.up==F,paste0(currMeasure,'.baseline')]-
                d.treatment[d.treatment$follow.up==F,paste0(currMeasure,'.w8')]
            
            l <- wilcox.test(diff.score, mu = 0)
            
            return(paste0('Delta-', currMeasure, ': V=', l$statistic, '; p = ',
                          signif(l$p.value,3), '\n\t'))
            
        } else {return(NULL)}
        
    } else {
        
        baseline <- sapply(unique(d.treatment$cimbi.id), function(j){
            matches <- d.treatment$cimbi.id==j & d.treatment$follow.up==F
            return(d.treatment[matches,i])
        })
        rescan <- sapply(unique(d.treatment$cimbi.id), function(j){
            matches <- d.treatment$cimbi.id==j & d.treatment$follow.up==T
            return(d.treatment[matches,i])
        })
        
        diff <- rescan-baseline
        diff <- diff[which(!is.na(diff))] # keep those which are not NA
        l <- wilcox.test(diff, mu = 0)
        
        return(paste0(i, ': V=', l$statistic, '; p = ',
                      signif(l$p.value,3), '\n\t'))
    }})),collapse='')

writeLines(paste0('**** GROUP DIFFERENCES: DESCRIPTIVES (NON-PARAMETRIC) ****\n\t', txt3))


####
## TREATMENT EFFECTS - UNIVARIATE (AAL)
####

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

# Figure 2 (see also code_figures.R)
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