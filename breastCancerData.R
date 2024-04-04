
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
?clinicalData
clinicalData$clinicalVarDef
clinicalData$clinicalTable$GEO_GSMID
data <- clinicalData$clinicalTable[c('dbUniquePatientID', 'study_ID', 'patient_ID', 'GEO_GSMID',  'age', 
                                     'path_diagnosis', 'path', 'tumor_stage_preTrt', 'hist_grade',
                                     'RFS', 'DFS', 'OS', 'dead', 'died_from_cancer_if_dead', 
                                     'ER_preTrt', 'PR_preTrt', #pre-treatment ER IHC status, pre-treatment progesterone IHC status
                                     'radiotherapyClass', 'chemotherapyClass', 'hormone_therapyClass',
                                     'chemotherapy', 'hormone_therapy', 'no_treatment', 
                                     'tamoxifen', 'anti_estrogen', 'aromatase_inhibitor', 'estrogen_receptor_blocker',
                                     'anti_HER2', 'taxaneGeneral', 'other'
)]

## looked through geo studies and found ones that involved tamoxifen
clinicalData$clinicalTable[which(clinicalData$clinicalTable$study_ID %in% c('1379', '6577','9893', '16391','17705')), ]
filtered_data <- data[which(data$study_ID %in% c('1379', '6577','9893', '16391','17705')),]
filtered_data <- data[which(data$study_ID %in% c('9893', '16391')),]

write.csv(filtered_data, file='/Users/Lindsey/gatech/S2024/ISYE6421/project_breastcancer_data.csv',)
?write.csv



## geting the microarray data for each study
names(curatedBreastDataExprSetList)

GSE16391 <- exprs(curatedBreastDataExprSetList$study_16391_GPL570_all)
GSE16391 <- as.data.frame.integer(GSE16391)
rows <- curatedBreastDataExprSetList$study_16391_GPL570_all@featureData$gene_symbol
GSE16391$gene_symbol <- rows
head(GSE16391)


GSE9893 <- exprs(curatedBreastDataExprSetList$study_9893_GPL5049_all)
GSE9893 <- as.data.frame.integer(GSE9893)
rows <- curatedBreastDataExprSetList$study_9893_GPL5049_all@featureData$gene_symbol
GSE9893$gene_symbol <- rows
head(GSE9893)


## getting the patient IDs which have relapse free survival of 0 and 1
RFS_0 <- filtered_data$patient_ID[which(filtered_data$RFS == 0)]
RFS_1 <- filtered_data$patient_ID[which(filtered_data$RFS == 1)]
length(RFS_0)
length(RFS_1)



# subsetting the microarray data for the patients who have RFS 0,1 to compare
GSE9893 <- GSE9893[,1]
study1_rfs_0 <- GSE9893[which(colnames(GSE9893) %in% RFS_0),]
study1_rfs_1 <- GSE9893[which(colnames(GSE9893) %in% RFS_1),]

# comparing the two groups
t.test(study1_rfs_0, study1_rfs_1, var.equal=T)




## same for study 2
GSE16391 <- GSE16391[,1]
study2_rfs_0 <- GSE16391[which(colnames(GSE16391) %in% RFS_0),]
study2_rfs_1 <- GSE16391[which(colnames(GSE16391) %in% RFS_1),]

t.test(study2_rfs_0, study2_rfs_1, var.equal=T)






citation("curatedBreastData")

