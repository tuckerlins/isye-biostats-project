
#BiocManager::install("curatedBreastData")
library(curatedBreastData)

data(curatedBreastDataExprSetList)
?curatedBreastDataExprSetList
exprs(curatedBreastDataExprSetList$study_17705_GPL96_MDACC_Tissue_BC_Tamoxifen)
dim(exprs(curatedBreastDataExprSetList$study_17705_GPL96_MDACC_Tissue_BC_Tamoxifen))

# name of study datasets
names(curatedBreastDataExprSetList)



# first 10 clinical variables of the first 3 patients
head(pData(curatedBreastDataExprSetList[[3]])[c(1:3),c(1:10)])

# top 5000 genes from the first 2 datasets
proc_curatedBreastDataExprSetList <- processExpressionSetList(exprSetList=curatedBreastDataExprSetList[1:2], outputFileDirectory=".",numTopVarGenes=5000)

# loading master clinical data table
data(clinicalData)
#?clinicalData
#clinicalData$clinicalVarDef
#clinicalData$clinicalTable$GEO_GSMID
data <- clinicalData$clinicalTable[c('dbUniquePatientID', 'study_ID', 'patient_ID', 'GEO_GSMID',  'age', 
                                     'path_diagnosis', 'path', 'tumor_stage_preTrt', 'hist_grade', 'tumor_size_cm_preTrt_preSurgery',
                                     'RFS', 'RFS_months_or_MIN_months_of_RFS', 'DFS', 'OS', 'dead', 'died_from_cancer_if_dead', 
                                     'ER_preTrt', 'PR_preTrt', #pre-treatment ER IHC status, pre-treatment progesterone IHC status
                                     'radiotherapyClass', 'chemotherapyClass', 'hormone_therapyClass',
                                     'chemotherapy', 'hormone_therapy', 'no_treatment', 
                                     'tamoxifen', 'anti_estrogen', 'aromatase_inhibitor', 'estrogen_receptor_blocker',
                                     'anti_HER2', 'taxaneGeneral', 'other'
)]

## looked through geo studies and found ones that involved tamoxifen
clinicalData$clinicalTable[which(clinicalData$clinicalTable$study_ID %in% c('1379', '6577','9893', '16391','17705')), ]
filtered_data <- data[which(data$study_ID %in% c('9893', '16391')),]

write.csv(filtered_data, file='/Users/Lindsey/gatech/S2024/ISYE6421/project_breastcancer_data_2.csv',)

## study 9893 has tumor size
## study 16391 has RFS months (survival curve?)
head(filtered_data[which(filtered_data$dead == 0),])

## geting the microarray data for each study
names(curatedBreastDataExprSetList)

GSE16391 <- exprs(curatedBreastDataExprSetList$study_16391_GPL570_all)
GSE16391 <- as.data.frame.integer(GSE16391)
GSE16391 <- GSE16391[,1]
rows <- curatedBreastDataExprSetList$study_16391_GPL570_all@featureData$gene_symbol
rownames(GSE16391) <- rows
head(GSE16391)


GSE9893 <- exprs(curatedBreastDataExprSetList$study_9893_GPL5049_all)
GSE9893 <- as.data.frame.integer(GSE9893)
GSE9893 <- GSE9893[,1]
rows <- curatedBreastDataExprSetList$study_9893_GPL5049_all@featureData$gene_symbol
rownames(GSE9893) <- rows
head(GSE9893)


## getting the patient IDs which have relapse free survival of 0 and 1
RFS_0 <- filtered_data$patient_ID[which(filtered_data$RFS == 0)]
RFS_1 <- filtered_data$patient_ID[which(filtered_data$RFS == 1)]
#length(RFS_0)
#length(RFS_1)

## tamoxifen
tam_0 <- filtered_data$patient_ID[which(filtered_data$tamoxifen == 0)]
tam_1 <- filtered_data$patient_ID[which(filtered_data$tamoxifen == 1)]

## chemotherapyClass
chemo_0 <- filtered_data$patient_ID[which(filtered_data$chemotherapyClass == 0)]
chemo_1 <- filtered_data$patient_ID[which(filtered_data$chemotherapyClass == 1)]

## radiotherapyClass
radio_0 <- filtered_data$patient_ID[which(filtered_data$radiotherapyClass == 0)]
radio_1 <- filtered_data$patient_ID[which(filtered_data$radiotherapyClass == 1)]

# progesterone receptor
PR_0 <- filtered_data$patient_ID[which(filtered_data$PR_preTrt == 0)]
PR_1 <- filtered_data$patient_ID[which(filtered_data$PR_preTrt == 1)]

# deaths
dead_0 <- filtered_data$patient_ID[which(filtered_data$dead == 0)]
dead_1 <- filtered_data$patient_ID[which(filtered_data$dead == 1)]


###### microarray stuff ########
head(GSE16391)
head(GSE9893)
dim(GSE9893)
dim(GSE16391)

#plot(GSE9893, GSE16391)
#plot(GSE16391)

# Calculate the confidence interval
confint(l.model, level=0.95)

### workaround for t.test error: data are essentially constant
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


#### study 9893 #####
## RFS
# subsetting the microarray data for the patients who have RFS 0,1 to compare
study9_rfs_0 <- GSE9893[,which(colnames(GSE9893) %in% RFS_0)]
study9_rfs_1 <- GSE9893[,which(colnames(GSE9893) %in% RFS_1)]
dim(study9_rfs_0)
dim(study9_rfs_1)

# t-test between RFS=0 and RFS=1 groups
pval_rfs <- c()
for (i in 1:22898) {
    x <- study9_rfs_0[i,]
    y <- study9_rfs_1[i,]
    t <- t.test(x, y, var.equal=T)
    pval_rfs <- c(pval_rfs, t$p.value)
}

# getting the genes with the smallest p values
s9_rfs_genes <- cbind.data.frame(rownames(GSE9893), pval_rfs)
s9_rfs_siggenes <- s9_rfs_genes[which(s9_rfs_genes$pval_rfs < 0.05),]
s9_rfs_siggenes <- s9_rfs_siggenes[order(s9_rfs_siggenes$pval_rfs),]
head(s9_rfs_siggenes)

# FDR corrected p-values
s9_rfs_genes <- s9_rfs_genes[order(s9_rfs_genes$pval_rfs),]
genes <- s9_rfs_genes$`rownames(GSE9893)`
pvals <- s9_rfs_genes$pval_rfs
fdr <- ( 22898 / 1:length(pvals)) * pvals
s9_rfs_data <- data.frame(genes, pvals, fdr)
dim(s9_rfs_data[which(s9_rfs_data$fdr < 0.05),])
dim(s9_rfs_data[which(s9_rfs_data$pvals < 0.05),])


## TAMOXIFEN: doesn't work, tamoxifen_0 = 0

## CEHMO: doesn't work, chemo_1 = 0

## RADIO 
# getting significant genes for patients with and without radiotherapy
study9_radio_0 <- study9_rfs_0[,which(colnames(study9_rfs_0) %in% radio_1)]
study9_radio_1 <- study9_rfs_1[,which(colnames(study9_rfs_1) %in% radio_1)]
dim(study9_radio_0)
dim(study9_radio_1)

# t-test 
pval_radio <- c()
for (i in 1:22898) {
    x <- study9_radio_0[i,]
    y <- study9_radio_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_radio <- c(pval_radio, t)
}

# getting the genes with the smallest p values
s9_radio_genes <- cbind.data.frame(rownames(GSE9893), pval_radio)
s9_radio_siggenes <- s9_radio_genes[which(s9_radio_genes$pval_radio < 0.05),]
s9_radio_siggenes <- s9_radio_siggenes[order(s9_radio_siggenes$pval_radio),]
head(s9_radio_siggenes)

# FDR
s9_radio_genes <- s9_radio_genes[order(s9_radio_genes$pval_radio),]
genes <- s9_radio_genes$`rownames(GSE9893)`
pvals <- s9_radio_genes$pval_radio
fdr <- ( 22898 / 1:length(pvals)) * pvals
s9_radio_data <- data.frame(genes, pvals, fdr)
s9_radio_data[which(s9_radio_data$fdr < 0.05),]


## PROGESTORONE RECEPTOR PR
# getting significant genes for patients with and without tamoxifen
study9_PR_0 <- study9_rfs_0[,which(colnames(study9_rfs_0) %in% PR_1)]
study9_PR_1 <- study9_rfs_1[,which(colnames(study9_rfs_1) %in% PR_1)]
dim(study9_PR_0)
dim(study9_PR_1)

# t-test 
pval_pr <- c()
for (i in 1:22898) {
    x <- study9_PR_0[i,]
    y <- study9_PR_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_pr <- c(pval_pr, t)
}

# getting the genes with the smallest p values
s9_pr_genes <- cbind.data.frame(rownames(GSE9893), pval_pr)
s9_pr_siggenes <- s9_pr_genes[which(s9_pr_genes$pval_pr < 0.05),]
s9_pr_siggenes <- s9_pr_siggenes[order(s9_pr_siggenes$pval_pr),]
head(s9_pr_siggenes)

# fdr
s9_pr_genes <- s9_pr_genes[order(s9_pr_genes$pval_pr),]
genes <- s9_pr_genes$`rownames(GSE9893)`
pvals <- s9_pr_genes$pval_pr
fdr <- ( 22898 / 1:length(pvals)) * pvals
s9_pr_data <- data.frame(genes, pvals, fdr)
sum(s9_pr_data$fdr < 0.05)
s9_pr_data[which(s9_pr_data$fdr < 0.05),]


## deaths
# getting significant genes for patients who died
study9_dead_0 <- GSE9893[,which(colnames(GSE9893) %in% dead_0)]
study9_dead_1 <- GSE9893[,which(colnames(GSE9893) %in% dead_1)]
dim(study9_dead_0)
dim(study9_dead_1)

# t-test 
pval_dead <- c()
for (i in 1:22898) {
    x <- study9_dead_0[i,]
    y <- study9_dead_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_dead <- c(pval_dead, t)
}

# getting the genes with the smallest p values
s9_dead_genes <- cbind.data.frame(rownames(GSE9893), pval_dead)
s9_dead_siggenes <- s9_dead_genes[which(s9_dead_genes$pval_dead < 0.05),]
s9_dead_siggenes <- s9_dead_siggenes[order(s9_dead_siggenes$pval_dead),]
head(s9_dead_siggenes)

# fdr
s9_dead_genes <- s9_dead_genes[order(s9_dead_genes$pval_dead),]
genes <- s9_dead_genes$`rownames(GSE9893)`
pvals <- s9_dead_genes$pval_dead
fdr <- ( 22898 / 1:length(pvals)) * pvals
s9_dead_data <- data.frame(genes, pvals, fdr)
sum(s9_dead_data$fdr < 0.05)
s9_dead_data[which(s9_dead_data$fdr < 0.05),]





##### microarray stuff for study 16391 #####
# relapse free survival: 0, 1
study1_rfs_0 <- GSE16391[,which(colnames(GSE16391) %in% RFS_0)]
study1_rfs_1 <- GSE16391[,which(colnames(GSE16391) %in% RFS_1)]
dim(study1_rfs_0)



# t-test between RFS=0 and RFS=1 groups
pval_rfs <- c()
for (i in 1:54696) {
    x <- study1_rfs_0[i,]
    y <- study1_rfs_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_rfs <- c(pval_rfs, t)
}

# getting the genes with the smallest p values
s1_rfs_genes <- cbind.data.frame(rownames(GSE16391), pval_rfs)
s1_rfs_siggenes <- s1_rfs_genes[which(s1_rfs_genes$pval_rfs < 0.05),]
s1_rfs_siggenes <- s1_rfs_siggenes[order(s1_rfs_siggenes$pval_rfs),]
head(s1_rfs_siggenes)
dim(s1_rfs_siggenes)

# fdr
s1_rfs_genes <- s1_rfs_genes[order(s1_rfs_genes$pval_rfs),]
genes <- s1_rfs_genes$`rownames(GSE16391)`
pvals <- s1_rfs_genes$pval_rfs
fdr <- ( 54696 / 1:length(pvals)) * pvals
s1_rfs_data <- data.frame(genes, pvals, fdr)
sum(s1_rfs_data$fdr < 0.05)
s1_rfs_data[which(s1_rfs_data$fdr < 0.05),]
head(s1_rfs_data[which(s1_rfs_data$pvals < 0.05),])
dim(s9_rfs_data[which(s9_rfs_data$fdr < 0.05),])


## TAMOXIFEN
# getting significant genes for patients who were treated with tamoxifen who had good vs bad outcome
study1_tam_0 <- study1_rfs_0[,which(colnames(study1_rfs_0) %in% tam_1)]
study1_tam_1 <- study1_rfs_1[,which(colnames(study1_rfs_1) %in% tam_1)]
dim(study1_tam_0)
dim(study1_tam_1)

# t-test 
pval_tam <- c()
for (i in 1:54696) {
    x <- study1_tam_0[i,]
    y <- study1_tam_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_tam <- c(pval_tam, t)
}

# getting the genes with the smallest p values
s1_tam_genes <- cbind.data.frame(rownames(GSE16391), pval_tam)
s1_tam_siggenes <- s1_tam_genes[which(s1_tam_genes$pval_tam < 0.05),]
s1_tam_siggenes <- s1_tam_siggenes[order(s1_tam_siggenes$pval_tam),]
head(s1_tam_siggenes)
dim(s1_tam_siggenes)

# fdr
s1_tam_genes <- s1_tam_genes[order(s1_tam_genes$pval_tam),]
genes <- s1_tam_genes$`rownames(GSE16391)`
pvals <- s1_tam_genes$pval_tam
fdr <- ( 54696 / 1:length(pvals)) * pvals
s1_tam_data <- data.frame(genes, pvals, fdr)
sum(s1_tam_data$fdr < 0.05)
s1_tam_data[which(s1_tam_data$fdr < 0.05),1]


## CEHMO
# getting significant genes for patients who recieved chemo who had good vs bad outcome
study1_chemo_0 <- study1_rfs_0[,which(colnames(study1_rfs_0) %in% chemo_1)]
study1_chemo_1 <- study1_rfs_1[,which(colnames(study1_rfs_1) %in% chemo_1)]
dim(study1_chemo_0)
dim(study1_chemo_1)

# t-test 
pval_chemo <- c()
for (i in 1:54696) {
    x <- study1_chemo_0[i,]
    y <- study1_chemo_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_chemo <- c(pval_chemo, t)
}

# getting the genes with the smallest p values
s1_chemo_genes <- cbind.data.frame(rownames(GSE16391), pval_chemo)
s1_chemo_siggenes <- s1_chemo_genes[which(s1_chemo_genes$pval_chemo < 0.05),]
s1_chemo_siggenes <- s1_chemo_siggenes[order(s1_chemo_siggenes$pval_chemo),]
head(s1_chemo_siggenes)

# fdr
s1_chemo_genes <- s1_chemo_genes[order(s1_chemo_genes$pval_chemo),]
genes <- s1_chemo_genes$`rownames(GSE16391)`
pvals <- s1_chemo_genes$pval_chemo
fdr <- ( 54696 / 1:length(pvals)) * pvals
s1_chemo_data <- data.frame(genes, pvals, fdr)
sum(s1_chemo_data$fdr < 0.05)
s1_chemo_data[which(s1_chemo_data$fdr < 0.05),]


## RADIO 
# getting significant genes for patients who recieved radiotherapy who had good vs bad outcome
study1_radio_0 <- study1_rfs_0[,which(colnames(study1_rfs_0) %in% radio_1)]
study1_radio_1 <- study1_rfs_1[,which(colnames(study1_rfs_1) %in% radio_1)]
dim(study1_radio_0)
dim(study1_radio_1)

# t-test 
pval_radio <- c()
for (i in 1:54696) {
    x <- study1_radio_0[i,]
    y <- study1_radio_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_radio <- c(pval_radio, t)
}

# getting the genes with the smallest p values
s1_radio_genes <- cbind.data.frame(rownames(GSE16391), pval_radio)
s1_radio_siggenes <- s1_radio_genes[which(s1_radio_genes$pval_radio < 0.05),]
s1_radio_siggenes <- s1_radio_siggenes[order(s1_radio_siggenes$pval_radio),]
head(s1_radio_siggenes)

# fdr
s1_radio_genes <- s1_radio_genes[order(s1_radio_genes$pval_radio),]
genes <- s1_radio_genes$`rownames(GSE16391)`
pvals <- s1_radio_genes$pval_radio
fdr <- ( 54696 / 1:length(pvals)) * pvals
s1_radio_data <- data.frame(genes, pvals, fdr)
sum(s1_radio_data$fdr < 0.05)
s1_radio_data[which(s1_radio_data$fdr < 0.05),]


## PROGESTORONE RECEPTOR PR
# getting significant genes for patients who were PR+ who had good vs bad outcome
study1_PR_0 <- study1_rfs_0[,which(colnames(study1_rfs_0) %in% PR_1)]
study1_PR_1 <- study1_rfs_1[,which(colnames(study1_rfs_1) %in% PR_1)]
dim(study1_PR_0)
dim(study1_PR_1)

# t-test
pval_pr <- c()
for (i in 1:54696) {
    x <- study1_PR_0[i,]
    y <- study1_PR_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_pr <- c(pval_pr, t)
}

# getting the genes with the smallest p values
s1_pr_genes <- cbind.data.frame(rownames(GSE16391), pval_pr)
s1_pr_siggenes <- s1_pr_genes[which(s1_pr_genes$pval_pr < 0.05),]
s1_pr_siggenes <- s1_pr_siggenes[order(s1_pr_siggenes$pval_pr),]
head(s1_pr_siggenes)

# fdr
s1_pr_genes <- s1_pr_genes[order(s1_pr_genes$pval_pr),]
genes <- s1_pr_genes$`rownames(GSE16391)`
pvals <- s1_pr_genes$pval_pr
fdr <- ( 54696 / 1:length(pvals)) * pvals
s1_pr_data <- data.frame(genes, pvals, fdr)
sum(s1_pr_data$fdr < 0.05)
s1_pr_data[which(s1_pr_data$fdr < 0.05),]


### significant genes ####
head(s1_tam_siggenes)

head(s1_chemo_siggenes)

head(s1_radio_siggenes)
head(s9_radio_siggenes)


rfs_sig <- Reduce(intersect, list(s1_rfs_genes$`rownames(GSE16391)`, s9_rfs_genes$`rownames(GSE9893)`))
rfs_sig <- rfs_sig[!is.na(rfs_sig)]

head(s1_rfs_data[which(s1_rfs_data$genes %in% rfs_sig),])
head(s9_rfs_data[which(s9_rfs_data$genes %in% rfs_sig),])


head(s1_rfs_siggenes[which(s1_rfs_siggenes$`rownames(GSE16391)` %in% rfs_sig),])
head(s9_rfs_siggenes[which(s9_rfs_siggenes$`rownames(GSE9893)` %in% rfs_sig),])



head(s1_pr_siggenes)
head(s9_pr_siggenes)
s1 <- s1_rfs_data[which(!is.na(s1_rfs_data$genes)),]
s1 <- s1[which(s1$fdr < 0.05),]
s9 <- s9_rfs_data[which(!is.na(s9_rfs_data$genes)),]
s9 <- s9[which(s9$fdr < 0.05),]

common_genes <- Reduce(intersect, list(s1$genes, s9$genes))
s1_sub <- s1[which(s1$genes %in% common_genes),]
s9_sub <- s9[which(s9$genes %in% common_genes),]

rfs_sig_genes <- merge(s1_sub, s9_sub, by="genes")
avg_fdr <- (rfs_sig_genes$fdr.x + rfs_sig_genes$fdr.y) / 2
rfs_sig_genes[order(avg_fdr),]

# Reduce(intersect, list(s1_rfs_data$genes, s9_rfs_data$genes))

#rfs_data <- merge(s1_rfs_data, s9_rfs_data, by="genes", )

## logistic  regression
data1 <- filtered_data[,c("study_ID", "age", "path_diagnosis", "hist_grade", "tumor_size_cm_preTrt_preSurgery", "RFS", 'RFS_months_or_MIN_months_of_RFS', "ER_preTrt", "PR_preTrt", "radiotherapyClass",
                          "chemotherapyClass", "tamoxifen", "aromatase_inhibitor", "estrogen_receptor_blocker")]
write.csv(filtered_data, file='/Users/Lindsey/isye-biostats-project/data1.csv',)

head(data1)
summary(data1)

logmodel1 <- glm(RFS ~ age + hist_grade + radiotherapyClass + chemotherapyClass + tamoxifen + aromatase_inhibitor + estrogen_receptor_blocker + 
                      ER_preTrt + PR_preTrt, family = 'binomial', data = data1)
summary(logmodel1)

# ignoring treatments
logmodel1.1 <- glm(RFS ~ age + tumor_size_cm_preTrt_preSurgery + ER_preTrt + PR_preTrt, family = 'binomial', data = data1)
summary(logmodel1.1)

# only considering sig confounding factors

logmodel2 <- glm(RFS ~ study_ID + hist_grade + PR_preTrt + radiotherapyClass + tamoxifen + radiotherapyClass * tamoxifen, family = 'binomial', data = data1)
summary(logmodel2)
anova(logmodel2, test="Chisq")

logmodel3 <- glm(RFS ~ study_ID + hist_grade + PR_preTrt + chemotherapyClass + tamoxifen + chemotherapyClass * tamoxifen, family = 'binomial', data = data1)
summary(logmodel3)
anova(logmodel3, test="Chisq")

logmodel4 <- glm(RFS ~ study_ID + hist_grade + PR_preTrt + radiotherapyClass * chemotherapyClass * tamoxifen, family = 'binomial', data=data1)
summary(logmodel4)
anova(logmodel4, test="Chisq")


logmodel5 <- glm(RFS ~ study_ID + age + hist_grade + ER_preTrt + PR_preTrt, family = 'binomial', data=data1)
summary(logmodel5)

logmodel6 <- glm(formula = RFS ~ hist_grade + PR_preTrt + chemotherapyClass * tamoxifen + radiotherapyClass * tamoxifen, family = 'binomial', data = filtered_data)
summary(logmodel6)
anova(logmodel6, test="Chisq")

logmodel7 <- glm(formula = RFS ~ hist_grade + PR_preTrt + tamoxifen, data = data1, family = 'binomial')
summary(logmodel7)

logmodel8 <- glm(formula = RFS ~ tamoxifen + radiotherapyClass + chemotherapyClass, family='binomial'(link="logit"), data=data1)
summary(logmodel8)
anova(logmodel8, test="Chisq")
model8.aov <- aov( RFS ~ tamoxifen , data=data1)
summary(model8.aov)
TukeyHSD(model8.aov, conf.level=0.95)
plot(TukeyHSD(model8.aov,conf.level=0.95))

colnames(data1)




####### summary stats ##########

## generating 2x2 tables for various tests
table(data1$RFS)
xtabs(~RFS,data=data1)
xtabs(~RFS + chemotherapyClass, data=data1)
xtabs(~RFS + tamoxifen, data=data1)
xtabs(~RFS + radiotherapyClass, data=data1)
xtabs(~RFS + study_ID, data=data1)
xtabs(~RFS + PR_preTrt, data=data1)
xtabs(~RFS + ER_preTrt, data=data1)
xtabs(~study_ID, data=data1)

## age and tumor size summary stats - box plots
data1_age_9893 <- data1$age[which(data1$study_ID == 9893)]
data1_age_16391 <- data1$age[which(data1$study_ID == 16391)]

boxplot(data1_age_9893, data1_age_16391, main="Age of patients", at=c(1,2), names=c("study 9893","study 16391"), ylab = "age")
boxplot(data1$tumor_size_cm_preTrt_preSurgery, main="Tumor size", xlab = "study 9893", ylab = "tumor size (cm)")

boxplot(data1$RFS_months_or_MIN_months_of_RFS, main="relapse free survival time", xlab = 'study 16391', ylab = 'months')

## bar graphs
barplot(table(data1$RFS), main = "Overall relapse free survival", xlab="RFS status", ylab = '# of patients')
barplot(xtabs(~tamoxifen + RFS, data=data1), beside=TRUE, 
        main = 'Tamoxifen treatment', legend = c('no tam','tamoxifen'), xlab='RFS', ylab='# of patients')

barplot(xtabs(~radiotherapyClass + study_ID, data=data1),beside=T,
        main = 'Radiotherapy across studies', legend = c('no radio','radiotherapy'), xlab="study ID", ylab='# of patients')

barplot(xtabs(~chemotherapyClass + study_ID, data=data1),beside=T,
        main = 'Chemotherapy across studies', legend = c('no chemo','chemotherapy'), xlab="study ID", ylab='# of patients')
xtabs(~radiotherapyClass + study_ID, data=data1)

barplot(xtabs(~ER_preTrt + study_ID, data=data1),beside=T,
        main = 'ER positivity prior to\ntreatment across studies', legend = c('ER-','ER+'), xlab="study ID", ylab='# of patients')

barplot(xtabs(~PR_preTrt + study_ID, data=data1),beside=T,
        main = 'PR positivity prior to\ntreatment across studies', legend = c('PR-','PR+'), xlab="study ID", ylab='# of patients')
xtabs(~ER_preTrt + study_ID, data=data1)
# i may have gone a lil wild with the graphs idk

citation("curatedBreastData")

