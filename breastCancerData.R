
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
length(RFS_0)
length(RFS_1)


## microarray stuff
head(GSE16391)
head(GSE9893)

#### study 9893 #####
# subsetting the microarray data for the patients who have RFS 0,1 to compare
study9_rfs_0 <- GSE9893[,which(colnames(GSE9893) %in% RFS_0)]
study9_rfs_1 <- GSE9893[,which(colnames(GSE9893) %in% RFS_1)]
dim(study9_rfs_0)
dim(study9_rfs_1)

study9_rfs_0[2,]

# t-test between RFS=0 and RFS=1 groups
pval_rfs <- c()
for (i in 1:22898) {
    x <- study9_rfs_0[i,]
    y <- study9_rfs_1[i,]
    t <- t.test(x, y, var.equal=T)
    pval_rfs <- c(pval_rfs, t$p.value)
}


# getting the genes with the smallest p values
s9_genes <- cbind.data.frame(rownames(GSE9893), pval_rfs)
s9_siggenes <- s9_genes[which(s9_genes$pval_rfs < 0.05),]
s9_siggenes <- s9_siggenes[order(s9_siggenes$pval_rfs),]
head(s9_siggenes)




##### same for study 16391 #####
# relapse free survival: 0, 1
study1_rfs_0 <- GSE16391[,which(colnames(GSE16391) %in% RFS_0)]
study1_rfs_1 <- GSE16391[,which(colnames(GSE16391) %in% RFS_1)]
dim(study1_rfs_0)
head(study1_rfs_1)

# 


### workaround for t.test error: data are essentially constant
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

# t-test between RFS=0 and RFS=1 groups
pval_rfs <- c()
for (i in 1:54696) {
    x <- study1_rfs_0[i,]
    y <- study1_rfs_1[i,]
    t <- my.t.test.p.value(x, y, var.equal=T)
    pval_rfs <- c(pval_rfs, t)
}

pval_rfs

# getting the genes with the smallest p values
s1_genes <- cbind.data.frame(rownames(GSE16391), pval_rfs)
s1_siggenes <- s1_genes[which(s1_genes$pval_rfs < 0.05),]
s1_siggenes <- s1_siggenes[order(s1_siggenes$pval_rfs),]
head(s1_siggenes)





## logistic  regression
data1 <- filtered_data[,c("study_ID", "age", "path_diagnosis", "hist_grade", "tumor_size_cm_preTrt_preSurgery", "RFS", "ER_preTrt", "PR_preTrt", "radiotherapyClass",
                          "chemotherapyClass", "tamoxifen", "aromatase_inhibitor", "estrogen_receptor_blocker")]
head(data1)
summary(data1)


logmodel1 <- glm(RFS ~ age + hist_grade + radiotherapyClass + chemotherapyClass + tamoxifen + aromatase_inhibitor + estrogen_receptor_blocker + 
                     ER_preTrt + PR_preTrt, family = 'binomial', data = data1)
summary(logmodel1)

logmodel2 <- glm(RFS ~ study_ID + age + hist_grade + PR_preTrt + radiotherapyClass + tamoxifen + radiotherapyClass * tamoxifen, family = 'binomial', data = data1)
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

logmodel8 <- glm(formula = RFS ~ hist_grade + PR_preTrt + tamoxifen, family='binomial'(link="logit"), data=data1)
summary(logmodel8)

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

## age and tumor size summary stats - box plots
data1_age_9893 <- data1$age[which(data1$study_ID == 9893)]
data1_age_16391 <- data1$age[which(data1$study_ID == 16391)]

boxplot(data1_age_9893, data1_age_16391, main="Age of patients", at=c(1,2), names=c("study 9893","study 16391"), ylab = "age")
boxplot(data1$tumor_size_cm_preTrt_preSurgery, main="Tumor size (study 9893)", xlab = "study 9893", ylab = "tumor size (cm)")









citation("curatedBreastData")

