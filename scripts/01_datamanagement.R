### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)



library(GEOquery)
library(SummarizedExperiment)
library(stringr)
library(TCGAbiolinks)
library(sesame) #v1.19.7 
sesameDataCache()



#### 
#
# For TCGA data sets (COAD and READ)
#
####

# Downloading TCGA data using TCGAbiolinks
# Change "COAD" for "READ" to download rectal cancer samples


# Extract the clinical data from TCGA using TCGAbiolinks

query <- GDCquery(project = "TCGA-COAD",   # "TCGA-READ"
                  data.category = "Clinical",
                  data.type = "Clinical Supplement")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
clinic <- as.data.frame(clinical.BCRtab.all$clinical_patient_coad)
#write.csv(clinic, "COAD_clinic.csv", row.names = F)


library(TCGAbiolinks)

query_met <- GDCquery(
  project = "TCGA-COAD",  # "TCGA-READ"
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450")


GDCdownload(query_met)
met <- GDCprepare(query_met, summarizedExperiment=T)
#saveRDS(met,"TCGA_met_COAD.rds")



# For the TCGA samples, we implemented exclusion criteria: 
# removing cases annotated with "Item in special subset", "History of unacceptable prior treatment related to a prior/other malignancy", 
# "Case submitted is found to be a recurrence after submission", "Neoadjuvant therapy", "Synchronous malignancy", and "Pathology outside specification". 
# Samples preserved in formalin-fixed paraffin-embedded (FFPE) form were excluded to maintain consistency in sample quality. 
# Duplicate samples were identified and removed based on their plate number.

# The final list of patients and their samples are saved in "Included_samples_Patients_COAD.csv" & "Included_samples_Patients_READ.csv"

patients <- read.csv("InputData/Included_samples_Patients_COAD.csv") # "Included_samples_Patients_READ.csv"
met_mat_M  <- readRDS("TCGA_met_COAD.rds") 
met_mat_M <- met_mat_M[, colnames(met_mat_M) %in% patients$samples]
saveRDS(met_mat_M,"TCGA_COAD_SeSame.rds")



#
# Download Mutational signature data and MSI status 
#

## Download the Signatures from dcc.icgc.org

download.file("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv", destfile="TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
SBS <-  read.csv("TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv", header=T)
# select only primary tumors
SBS <- SBS[substr(SBS$Sample.Names, 14,15) == "01", ]
SBS$ID <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1", SBS$Sample.Names)
SBS <- SBS[SBS$ID %in% colnames(met_mat_M), ]


#
# MSI status was obtained from:
#

# Genetic Mechanisms of Immune Evasion in Colorectal Cancer 
# Catherine S. Grasso; Marios Giannakis; Daniel K. Wells; Tsuyoshi Hamada; Xinmeng Jasmine Mu; Michael Quist; Jonathan A. Nowak; Reiko Nishihara; Zhi Rong Qian; Kentaro Inamura;
# Teppei Morikawa; Katsuhiko Nosho; Gabriel Abril-Rodriguez; Charles Connolly; Helena Escuin-Ordinas; Milan S. Geybels; William M. Grady ;
# Li Hsu; Siwen Hu-Lieskovan; Jeroen R. Huyghe; Yeon Joo Kim; Paige Krystofinski; Mark D.M. Leiserson; Dennis J. Montoya; Brian B. Nadel; Matteo Pellegrini; Colin C. Pritchard; Cristina Puig-Saus; Elleanor H. Quist; Ben J. Raphael; Stephen J. Salipante; Daniel Sanghoon Shin; Eve Shinbrot; Brian Shirts; Sachet Shukla; Janet L. Stanford; Wei Sun; Jennifer Tsoi; Alexander Upfill-Brown; David A. Wheeler; Catherine J. Wu; Ming Yu; Syed H. Zaidi; Jesse M. Zaretsky; Stacey B. Gabriel; Eric S. Lander; Levi A. Garraway; Thomas J. Hudson; Charles S. Fuchs; Antoni Ribas; Shuji Ogino; Ulrike Peters
# Cancer Discov (2018) 8 (6): 730–749.
# https://doi.org/10.1158/2159-8290.CD-17-1327


MSS <- readr::read_tsv("215982clinic_sup_grassso_cancerDisc18.tsv")
MSS <- as.data.frame(MSS[,c("Sample", "MsiStatus")])
names(MSS)[1] <- "ID"

df_list <- list(clinic, SBS, MSS) 
clinic <- as.data.frame(Reduce(function(x, y) merge(x, y, by= "ID"), df_list))

clinic <- clinic[clinic$histologic_diagnosis == "Colon Adenocarcinoma", ]
write.csv(clinic, "Clinic_Coad.csv")

## Extract patient data information from the GEO, we only use Tumor samples in this study


my_id <- "GSE131013"   # GSE101764, GSE131013, GSE77954, GSE39958, GSE42752
gset <- getGEO(my_id, getGPL=FALSE)
gse <- gset[[1]]
sampleInfo <- pData(gse)
names(sampleInfo)
sampleInfo2 <- sampleInfo[,c("geo_accession", "age:ch1", "gender:ch1", "tissue:ch1")]
table(sampleInfo2$`tissue:ch1`)
sampleTumor <- sampleInfo2[sampleInfo2$tissue.ch1 == "Tumor", ] # we only include tumor samples
#write.csv(sampleTumor, "GSE131013_clinic.csv", row.names = T)




## Extract the DNA methylation raw data from GEO
getGEOSuppFiles("GSE131013")


#
# For GEO data sets that have idat files available (GSE101764, GSE131013, GSE77954)
#



untar("GSE131013_RAW.tar", exdir = "idat")
head(list.files("idat", pattern = "idat"))

idatFiles <- list.files("/GSE131013/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

idat <- list.files("/GSE131013/idat", pattern = "idat$", full = TRUE)
idat <- (sub('_Grn.idat','', idat))
idat <- (sub('_Red.idat','', idat))
idat <- unique(idat)
samples <- as.data.frame(idat)
samples$id <- gsub("_.*", "\\1", samples$idat)
samples$id <- str_remove(samples$id, "/GSE131013/idat/")


pdata <- read.csv("pData_GSE131013_tumor.csv")
samples <- samples[samples$id %in% pdata$geo_accession, ] # only select tumor samples

idattumor <- samples$idat

betas <- openSesame(idattumor)
saveRDS(betas, "GSE131013_SeSame.rds")


#
# For GEO data sets that have only signal file available (GSE39958, GSE42752)
#

getGEOSuppFiles("GSE39958")
path = 'GSE39958_matrix_signal.csv.gz'
Data = read.csv(gzfile(path), header=T)

Probes <- Data$ID_REF
rownames(Data) <- Data$ID_REF

library(tidyr)
library(stringr)

samples <- as.data.frame(Data[, grep("Unmethylated", colnames(Data))])
samples <- colnames(samples)
samples <-  str_remove(samples, ".Unmethylated.Signal")

unmet <- as.data.frame(Data[, grep("Unmethylated", colnames(Data))])
unmet <- sapply(unmet,as.numeric)
rownames(unmet) <- Data$ID_REF
unmet <- t(unmet)
row.names(unmet) <-  str_remove(rownames(unmet) , ".Unmethylated.Signal")

met <- as.data.frame(Data[, grep("Methylated", colnames(Data))])
met <- sapply(met,as.numeric)
rownames(met) <- Data$ID_REF
met <- t(met)
rownames(met) <-  str_remove(rownames(met) , ".Methylated.Signal")

a <- list()

for (i in samples){
  meti <- met[i, ]
  unmeti <- unmet[i, ]

  ssets <- parseGEOsignalMU(sigM =  meti, sigU= unmeti,Probe_IDs=Probes,
                            oob.mean = 500,
                            oob.sd = 300, 
                            platform = "HM450")
  
  
  b <- list(ssets)
  a <- c(a,b)
}

names(a) <- samples

sdf_preped = openSesame(a, func=getBetas)
saveRDS(sdf_preped, "GSE39958_SeSame.rds")




##### 
#
# For all data sets 
#
#####

library(SummarizedExperiment)
library(MultiAssayExperiment)
library(impute)
library(lumi)



DNAm <- readRDS("GSE39958_SeSame.rds")  # Change for each of the data sets the input file
met_mat <- as.data.frame(t(DNAm))
missCpG <- as.data.frame(colSums(is.na(met_mat)))
names(missCpG)[1] <- "missings"
missCpG$CpGs <- row.names(missCpG)

# Select more than 50 % missing 
allCpGsNA <- missCpG[missCpG$missings >= 23, ] # ADAPT based on sample size
DNAm3 <- DNAm[!rownames(DNAm) %in% allCpGsNA$CpGs, ]

# Matrix with CpGs in the rows, samples in the columns
# impute the CpGs that have less than 50% missing values.

datMethUsed2 = impute.knn(as.matrix(DNAm3))
datMethUsed <- as.data.frame(t(datMethUsed2$data))

# Select imputes as 0 for mean imputation
zero <- as.data.frame(colSums(datMethUsed==0))
data <- datMethUsed
data[data == 0] <- NA
zero$CpG <- rownames(zero)
CpGs <- zero[zero$`colSums(datMethUsed == 0)` >0, ]
data <- data[, colnames(data) %in% CpGs$CpG]


for(i in 1:ncol(data)){
  data[is.na(data[,i]), i] <- mean(data[,i], na.rm = TRUE)
}

datMeth <- datMethUsed[,! colnames(datMethUsed) %in% CpGs$CpG]
DNAm <- cbind(datMeth, data)
#saveRDS(DNAm, "GSE39958imp_knn_mean.rds")  # Change name of the cohort




#
# Calculate Horvath clock epigenetic age using the R tutorial provided in 
# https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM20_ESM.docx
# 
# Download needed files from: 
# Horvath, S. DNA methylation age of human tissues and cell types. 
# Genome Biol 14, 3156 (2013). https://doi.org/10.1186/gb-2013-14-10-r115
# 
# 1.	AdditionalFile21datMiniAnnotation27k.csv
# 2.	AdditionalFile22probeAnnotation21kdatMethUsed.csv
# 3.	AdditionalFile23predictor.csv
# 4.	AdditionalFile24NORMALIZATION.R.txt
# 5.	AdditionalFile25StepwiseAnalysis.txt



DNAm <- readRDS("GSE39958_SeSame.rds")

library(impute)
# Matrix with CpGs in the rows, samples in the columns
datMethUsed <- t(DNAm)
datMethUsed2 <- impute.knn(as.matrix(DNAm))
datMethUsed <- as.data.frame(t(datMethUsed2$data))


# Normalize the imputed Beta values using the BMIQ method from Horvarth. 
# Additional file 24: R code for normalizing the DNA methylation data.

source("AdditionalFile24NORMALIZATION.R")

## Additional file 22: Additional probe annotation file for the R tutorial. 
probeAnnotation21kdatMethUsed2 = read.table("AdditionalFile22probeAnnotation21kdatMethUsed.csv",header = T,sep=",",stringsAsFactors = F)
datMethUsed2Normalized=BMIQcalibration(datM=datMethUsed[,probeAnnotation21kdatMethUsed2$Name],goldstandard.beta= probeAnnotation21kdatMethUsed2$goldstandard2,plots=F)

trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

##Additional file 23: Coefficient values of the age predictor.
datClock=read.csv("AdditionalFile23predictor.csv")

## Code from https://horvath.genetics.ucla.edu/html/dnamage/StepwiseAnalysis.txt
selectCpGsClock=is.element(dimnames(datMethUsed2Normalized)[[2]], as.character(datClock$CpGmarker[-1]))
datMethClock0=data.frame(datMethUsed2Normalized[,selectCpGsClock])
datMethClock= data.frame(datMethClock0[ as.character(datClock$CpGmarker[-1])])
predictedAge=as.numeric(anti.trafo(datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(datClock$CoefficientTraining[-1])))
names(predictedAge) = rownames(datMethUsed2Normalized)


### From beta to M values and Removal of clock methylation
met_mat <- readRDS("GSE39958imp_knn_mean.rds")
met_mat <- as.data.frame(t(met_mat))
met_mat_beta = lumi::beta2m(met_mat)

adjusted_met_mat_beta = apply(met_mat_beta,1,function(x) lm(x~predictedAge)$residuals+mean(x)   )
adjusted_met_mat_beta = t(adjusted_met_mat_beta)
adjusted_met_mat = lumi::m2beta(adjusted_met_mat_beta)

#saveRDS(adjusted_met_mat, "GSE39958_Age_adjusted_DNAm_Bval.rds")


# We only include CpGs that are present in all data sets and remove non-specific probes and CpGs that are SNPs 
# The final included CpGs are in the InputData folder

KeepCpGs <- read.csv("InputData/overlap_CpGs.csv", header=T)
GSE39958 <- GSE39958[rownames(GSE39958) %in% KeepCpGs$CpGs, ] 

#saveRDS(adjusted_met_mat, "GSE39958_Age_adjusted_DNAm_Bval_Overlap.rds")






