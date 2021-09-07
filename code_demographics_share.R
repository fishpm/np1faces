#### SUMMARY: Computation of summary measures for demographic and related data.

#### AUTHOR: Patrick M. Fisher, 2021

####
## READ DATA FRAME
####

d0 <- read.csv('./pubMaterials/data_demo_aal.csv', stringsAsFactors = F)

####
## BASELINE ANALYSIS
####

# define baseline data frame
d.baseline <- d0[d0$exclude.baseline=='',]
d.baseline$age0 <- d.baseline$age-mean(d.baseline$age)

####
## DESCRIPTIVES
####

# clinical groups
writeLines('**** GROUP SAMPLE SIZES ****')
table(d.baseline$group.status)

# sex prevalence and distribution
writeLines('**** TEST FOR DISCREPANCY IS SEX DISTRIBUTION BETWEEN GROUPS ****')
table(d.baseline$group.status, d.baseline$sex)
chisq.test(d.baseline$group.status, d.baseline$sex)

response.groups <- c('np1.treatresponse','remitter.status','responder.status')
for (i in response.groups){
    writeLines(paste0('\nTable: ', i))
    print(table(d.baseline[d.baseline$compliance.w8=='Yes',i]))
}

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

# groupings of interest
groups <- c('All Participants','Case','Healthy Control')


####
## Data are foundation for Table 1
####

writeLines('**** GROUP DESCRIPTIVE INFORMATION ****')
for(g in groups){
    if (g == 'All Participants'){
        row.set <- seq(nrow(d.baseline))
    } else {
        row.set <- which(d.baseline$group.status==g)
    }
    
    txt <- paste(unlist(lapply(desc.cont, function(i){
        nomit <- sum(is.na(d.baseline[row.set,i]))
        if(i%in%c('hamd6.w8','hamd17.w8') & g == 'Case'){
            row.set <- which(d.baseline$exclude.treatresp=='' & 
                                 d.baseline$group.status=='Case')
            nomit <- sum(d.baseline$group.status=='Case') - length(row.set)
        }
        return(paste0(i, ': ', 
                      signif(mean(d.baseline[row.set,i], na.rm=T),3), ' \u00B1 ',
                      signif(sd(d.baseline[row.set,i], na.rm=T),3), ' [',
                      signif(median(d.baseline[row.set,i], na.rm=T),3), '; ',
                      signif(min(d.baseline[row.set,i], na.rm=T),3), '-',
                      signif(max(d.baseline[row.set,i], na.rm=T),3), ']',
                      ' (nomit: ', nomit, ')\n\t'))
    })),collapse='')
    
    writeLines(paste0('**** ', g, ' ****\n\t',txt))
}

txt2 <- paste(unlist(lapply(desc.cont, function(i){
    if(i%in%clin.measures){
        
        if(i%in%c('hamd17.baseline','hamd6.baseline')){
            currMeasure <- unlist(strsplit(i,'.',fixed=T))[1]
            
            matches <- d.baseline$exclude.treatresp=='' & d.baseline$group.status=='Case'
            
            diff.score <- d.baseline[matches,paste0(currMeasure,'.baseline')]-
                d.baseline[matches,paste0(currMeasure,'.w8')]
            
            ttest <- t.test(diff.score)
            
            return(paste0('Delta-', currMeasure, ': ', signif(ttest$estimate,3), '[',
                          paste(signif(ttest$conf.int,3),collapse=', '), ']; p = ',
                          signif(ttest$p.value,3), '\n\t'))
        } else {return(NULL)}
        
    } else {
        f <- as.formula(paste0(i, '~group.status'))
        l <- t.test(f, d.baseline)
        return(paste0(i, ': ', signif(-diff(l$estimate),3), '[',
                      paste(signif(l$conf.int,3),collapse=', '), ']; p = ',
                      signif(l$p.value,3), '\n\t'))
    }})),collapse='')

writeLines(paste0('**** GROUP DIFFERENCES: DESCRIPTIVES (PARAMETRIC) ****\n\t', txt2))

txt3 <- paste(unlist(lapply(desc.cont, function(i){
    if(i%in%clin.measures){
        
        if(i%in%c('hamd17.baseline','hamd6.baseline')){
            currMeasure <- unlist(strsplit(i,'.',fixed=T))[1]
            
            diff.score <- d.baseline[,paste0(currMeasure,'.baseline')]-
                d.baseline[,paste0(currMeasure,'.w8')]
            
            l <- wilcox.test(diff.score, mu = 0)
            
            return(paste0('Delta-', currMeasure, ': V=', l$statistic, '; p = ',
                          signif(l$p.value,3), '\n\t'))
            
        } else {return(NULL)}
        
    } else {
        
        f <- as.formula(paste0(i, '~group.status'))
        l <- wilcox.test(f, d.baseline)
        
        return(paste0(i, ': W=', l$statistic, '; p = ',
                      signif(l$p.value,3), '\n\t'))
    }})),collapse='')

writeLines(paste0('**** GROUP DIFFERENCES: DESCRIPTIVES (NON-PARAMETRIC) ****\n\t', txt3))


## reaction time correlates
faces.rt <- cor.test(d.baseline$srt.rt,d.baseline$rt.faces.all, na.action = na.omit)
sum(!is.na(d.baseline$srt.rt) & !is.na(d.baseline$rt.faces.all))
shapes.rt <- cor.test(d.baseline$srt.rt,d.baseline$rt.shapes.all, na.action = na.omit)
sum(!is.na(d.baseline$srt.rt) & !is.na(d.baseline$rt.shapes.all))
table(d.baseline$group.status[!is.na(d.baseline$srt.rt) & !is.na(d.baseline$rt.shapes.all)])


