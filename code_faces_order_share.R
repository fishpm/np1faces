#### SUMMARY: Generation of data file (data_demo_aal.csv) containing organized demographic, and regional fMRI responses to faces vs. shapes contrast. This code is run only once to generate the data_demo_aal.csv file.

#### AUTHOR: Patrick M. Fisher, 2021

####
## LOAD PACKAGES
####

library(nlme)

####
## ANALYSES
####

# read in previuosly generated data frame (code_dataLoad_faces_order.R)
d0 <- read.csv('./pubMaterials/data_demo_aal_faces_order.csv', stringsAsFactors = F)

# define factors
d0$emotion.factor <- factor(d0$emotion.factor, levels = c('F','N','A','S'))
d0$faces.block.factor <- factor(d0$faces.block)
d0$sex.factor <- factor(d0$sex, levels = c('Female','Male'))
d0$group.status <- factor(d0$group.status, levels = c('Healthy Control', 'Case'))

# region name text file
aal.names <- as.vector(unlist(read.delim('./pubMaterials/atlas116.txt',sep = '\t', stringsAsFactors = F)))

# rename two AAL region names redundant with other names
aal.names[25] <- 'Frontal_Med_Orb_L' # was Frontal_Mid_Orb_L, which duplicates 9
aal.names[26] <- 'Frontal_Med_Orb_R' # was Frontal_Mid_Orb_R, which duplicates 10

# cortical and subcortical regions (cerebellum excluded)
aal.nregions <- 90

####
## MAIN EFFECT OF HABITUATION
####

covar <- c('age0','sex.factor','emotion.factor','faces.block.factor', 'rt.faces.all0','rt.shapes.all0')
out <- list()

# region-specific analyses
for(region in aal.names[seq(aal.nregions)]){
    
    writeLines(paste0('Working on ', region))
    
    # mixed-effects model
    f <- as.formula(paste0(region,'~',paste(covar,collapse='+')))
    l <- gls(f,
             correlation = corSymm(form = ~1|cimbi.id),
             weight = varIdent(form = ~1|faces.block.factor),
             data = d0)
    
    # read out model parameters
    anova.out <- anova(l)
    anova.vals <- c(anova.out['emotion.factor','F-value'],
                    anova.out['emotion.factor','p-value'],
                    anova.out['faces.block.factor','F-value'],
                    anova.out['faces.block.factor','p-value'],
                    anova.out['(Intercept)','F-value'],
                    anova.out['(Intercept)','p-value'])
    
    # build predicted outcomes
    matches <- d0[,'cimbi.id']==d0[1,'cimbi.id']
    df.predict.block <- d0[matches,covar]
    df.predict.block$age0 <- 0
    df.predict.block$rt.faces.all0 <- 0
    df.predict.block$rt.shapes.all0 <- 0
    df.predict.block$sex.factor <- factor('Female', levels = c('Female','Male')) # for a female
    df.predict.block$emotion.factor <- factor('F', levels = c('F','N','A','S')) # fear task
    df.predict.block$region <- region
    df.predict.block$predicted <- predict(l, newdata = df.predict.block) # block-specific predicted
    
    df.predict.emotion <- d0[matches,covar]
    df.predict.emotion$age0 <- 0
    df.predict.emotion$rt.faces.all0 <- 0
    df.predict.emotion$rt.shapes.all0 <- 0
    df.predict.emotion$sex.factor <- factor('Female', levels = c('Female','Male')) # for a femaile
    df.predict.emotion$faces.block.factor <- factor(1, levels = seq(4)) # first block
    df.predict.emotion$region <- region
    df.predict.emotion$predicted <- predict(l,newdata = df.predict.emotion) # emotion-specific predicted
    
    # return predicted estimates to list element
    out[[region]] <- list(df.predict.block=df.predict.block,
                          df.predict.emotion=df.predict.emotion,
                          anova.vals=anova.vals)
    
}

# collapse list elements
out.block <- do.call(rbind,lapply(seq(length(out)), function(i){out[[i]]$df.predict.block}))
out.emotion <- do.call(rbind,lapply(seq(length(out)), function(i){out[[i]]$df.predict.emotion}))

# effects across regions
df.anova <- data.frame(t(sapply(seq(length(out)), function(i){out[[i]]$anova.vals})))
colnames(df.anova) <- c('emotion.fstat', 'emotion.pval', 'block.fstat', 'block.pval', 'intercept.fstat','intercept.pval')
df.anova$regions <- aal.names[seq(aal.nregions)]
df.anova$emotion.pval.adj <- p.adjust(df.anova$emotion.pval, method = 'holm')
df.anova$block.pval.adj <- p.adjust(df.anova$block.pval, method = 'holm')
df.anova$intercept.pval.adj <- p.adjust(df.anova$intercept.pval, method = 'holm')

# main effect of faces block
p <- ggplot(out.block, aes(x=faces.block.factor,y=predicted))
p + geom_boxplot() +
    labs(title = 'Task habituation',
         x = 'Faces Block',
         y = 'Brain response (arbitrary units, predicted)') +
    scale_y_continuous(breaks = c(-0.25,0,0.25,0.50,0.75),
                       limits = c(-0.25, 0.75)) +
    theme(plot.title = element_text(hjust=0.5),
          legend.title.align = 0.5,
          axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))

# main effect of emotion block
p <- ggplot(out.emotion, aes(x=emotion.factor,y=predicted))
p + geom_boxplot() +
    labs(title = 'Emotion Effects',
         x = 'Emotion type',
         y = 'Brain response (arbitrary units, predicted)') +
    scale_y_continuous(breaks = c(-0.25,0,0.25,0.5,0.75),
                       limits = c(-0.25, 0.75)) +
    theme(plot.title = element_text(hjust=0.5),
          legend.title.align = 0.5,
          axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))

####
## Data are foundation of "Faces Block Effects" and "Emotion Type Effects" columns in Supplementary Table 2
####

write.csv(df.anova, './pubMaterials/aal_faces-block_emotion_maineffects.csv',row.names=F)


####
## HABITUATION MODERATED BY GROUP
####

covar <- c('age0','sex.factor','emotion.factor', 'rt.faces.all0','rt.shapes.all0')
int.terms <- c('group.status','faces.block.factor')
out <- list()

# region-specific analyses
for(region in aal.names[seq(aal.nregions)]){
    
    writeLines(paste0('Working on ', region))
    
    # mixed-effects interaction model
    f <- as.formula(paste0(region,'~',
                           paste(c(covar,paste(int.terms,collapse='*')),collapse='+')))
    l <- gls(f,
             correlation = corSymm(form = ~1|cimbi.id),
             weight = varIdent(form = ~1|faces.block.factor),
             data = d0)
    
    # read out model parameters for interaction effect
    anova.out <- anova(l)
    anova.vals <- c(anova.out['group.status:faces.block.factor','F-value'],
                    anova.out['group.status:faces.block.factor','p-value'])
    
    # build predicted outcomes
    df.tmp <- d0
    g1.ex <- df.tmp[df.tmp$cimbi.id==df.tmp$cimbi.id[which(df.tmp$group.status=='Healthy Control')[1]],c(covar,int.terms)]
    g2.ex <- df.tmp[df.tmp$cimbi.id==df.tmp$cimbi.id[which(df.tmp$group.status=='Case')[1]],c(covar,int.terms)]
    df.predict <- rbind(g1.ex,g2.ex)
    df.predict$age0 <- 0
    df.predict$rt.faces.all0 <- 0
    df.predict$rt.shapes.all0 <- 0
    df.predict$sex.factor <- factor('Female', levels = c('Female','Male'))
    df.predict$emotion.factor <- factor('F', levels = c('F','N','A','S'))
    df.predict$region <- region
    df.predict$predicted <- predict(l,newdata = df.predict)
    
    # return predicted estimates to list element
    out[[region]] <- list(df.predict=df.predict,
                          anova.vals=anova.vals)
}

# collapse list elements
out2 <- do.call(rbind,lapply(seq(length(out)), function(i){out[[i]]$df.predict}))

# effects across regions
df.anova <- data.frame(t(sapply(seq(length(out)), function(i){out[[i]]$anova.vals})))
colnames(df.anova) <- c('interaction.fstat', 'interaction.pval')
df.anova$regions <- aal.names[seq(aal.nregions)]
df.anova$interaction.pval.adj <- p.adjust(df.anova$interaction.pval, method = 'holm')


# FIGURE 2 PLOT (part 1 of 2)
p <- ggplot(out2, aes(x=faces.block.factor,y=predicted,color=group.status))
p +
    geom_boxplot() +
    labs(title = 'Task habituation',
         x = 'Faces Block',
         y = 'Brain response (arbitrary units, adjusted)',
         color = 'Group') +
    scale_y_continuous(breaks = c(-0.25,0,0.25,0.5,0.75),
                       limits = c(-0.25,0.75)) +
    scale_color_discrete(labels = c('Healthy','Depressed')) +
    theme(plot.title = element_text(hjust=0.5),
          legend.title.align = 0.5,
          axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))

####
## Data are foundation of "Group-x-Block Effects" columns in Supplementary Table 2
####

write.csv(df.anova, './pubMaterials/aal_faces-block_interaction.csv',row.names=F)


####
## EMOTION MODERATED BY GROUP
####

covar <- c('age0','sex.factor','faces.block.factor', 'rt.faces.all0','rt.shapes.all0')
int.terms <- c('group.status','emotion.factor')
out <- list()

# region-specific analyses
for(region in aal.names[seq(aal.nregions)]){
    
    writeLines(paste0('Working on ', region))
    
    # mixed-effects interaction model
    f <- as.formula(paste0(region,'~',
                           paste(c(covar,paste(int.terms,collapse='*')),collapse='+')))
    l <- gls(f,
             correlation = corSymm(form = ~1|cimbi.id),
             weight = varIdent(form = ~1|emotion.factor),
             data = d0)
    
    anova.out <- anova(l)
    anova.vals <- c(anova.out['group.status:emotion.factor','F-value'],
                    anova.out['group.status:emotion.factor','p-value'])
    
    # build predicted outcomes
    df.tmp <- d0
    g1.ex <- df.tmp[df.tmp$cimbi.id==df.tmp$cimbi.id[which(df.tmp$group.status=='Healthy Control')[1]],c(covar,int.terms)]
    g2.ex <- df.tmp[df.tmp$cimbi.id==df.tmp$cimbi.id[which(df.tmp$group.status=='Case')[1]],c(covar,int.terms)]
    df.predict <- rbind(g1.ex,g2.ex)
    df.predict$age0 <- 0
    df.predict$rt.faces.all0 <- 0
    df.predict$rt.shapes.all0 <- 0
    df.predict$sex.factor <- factor('Female', levels = c('Female','Male'))
    df.predict$faces.block.factor <- factor(1, levels = seq(4))
    df.predict$region <- region
    df.predict$predicted <- predict(l,newdata = df.predict)
    
    # return predicted estimates to list element
    out[[region]] <- list(df.predict=df.predict,
                          anova.vals=anova.vals)
}

# collapse list elements
out2 <- do.call(rbind,lapply(seq(length(out)), function(i){out[[i]]$df.predict}))

# effects across regions
df.anova <- data.frame(t(sapply(seq(length(out)), function(i){out[[i]]$anova.vals})))
colnames(df.anova) <- c('interaction.fstat', 'interaction.pval')
df.anova$regions <- aal.names[seq(aal.nregions)]
df.anova$interaction.pval.adj <- p.adjust(df.anova$interaction.pval, method = 'holm')


# FIGURE 2 PLOT (part 2 of 2)
p <- ggplot(out2, aes(x=factor(emotion.factor,level=c('F','A','S','N')),
                      y=predicted,color=group.status))
p +
    geom_boxplot() +
    labs(title = 'Task habituation',
         x = 'Emotion Block',
         y = 'Brain response (arbitrary units, adjusted)',
         color = 'Group') +
    scale_y_continuous(breaks = c(-0.25,0,0.25,0.5,0.75),
                       limits = c(-0.25,0.75)) +
    scale_color_discrete(labels = c('Healthy','Depressed')) +
    scale_x_discrete(breaks = c('F','N','A','S'),
                     labels = c('Fear','Neutral','Angry','Surprise')) +
    theme(plot.title = element_text(hjust=0.5),
          legend.title.align = 0.5,
          axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))

####
## Data are foundation of "Emotion-x-Block Effects" columns in Supplementary Table 2
####

write.csv(df.anova, './pubMaterials/aal_emotion_interaction.csv',row.names=F)