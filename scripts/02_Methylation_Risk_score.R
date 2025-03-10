### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)





# We only include patients with "Colon Adenocarcinoma" and have removed [Discrepancy], [Not Available], Colon Mucinous Adenocarcinoma 
# This matrix includes all the CpGs that were included in any of the Exposome traits sets

met_mat_M  <- readRDS("Startdata/COAD_Age_adjusted_DNAm_Bval_overlap_MRSCpGs.rds")   # Do this for each input data set independent
clinic <- read.csv("Startdata/Clinic_Coad.csv", header=T)
head(clinic)

clinic$EOCRC = as.factor(ifelse(clinic$age < 50 , 1,
                                ifelse(clinic$age >=70, 0, NA)))

source("Functions/MRS_GLM_function.R")

dir.create("MRS")                              # Important to add as the MRSs will be saved in this folder
dir.create("GLM_MRS")                          # Important to add as the results of the GLMs will be saved   


GLM_MRS(datameth = met_mat_M,                  # matrix with CpGs in rows and  patients in columns
        dataclinic = clinic, 
      #  histodiag = "Colon Adenocarcinoma",    # Can be left out,  We removed [Discrepancy], [Not Available], Colon Mucinous Adenocarcinoma 
        PhenoCpG = "GW",                       # Options: "GW", "P1E5", "F01", "F005", "F001"
        datasetname = "COAD",                  # This name will be added to the output file (no other function) 
        adjust = "gender"                      # Can be left out, if not enough patients are available gender adjustment will be removed
)


#
# The MRSs obtained and used in the manuscript are located in "Results/MRS/"
# The results obtained comparing early-onset vs later-onset patients are provided in "Results/GLM_MRS/"
#



#
# Permutation CpGs
#

dir.create("GLM_MRS_perm")                     # Important to add as the results of the GLMs will be saved  

source("Functions/FunctionMRS_GLM_perm.R")

GLM_MRS_perm(datameth = met_mat_M, 
             dataclinic = clinic, 
             histodiag = "Colon Adenocarcinoma",
             PhenoCpG = "GW",               #  Options: "GW", "P1E5", "F01", "F005", "F001"
             datasetname = "COAD", 
             adjust = "gender"
)



#
# The CpG permutation results obtained comparing early-onset vs later-onset patients are provided in "Results/GLM_MRS_perm/"
#




#
# Permutation for patient categorization
#

# To use own generated MRSs use the output in directory "MRS"  
# To use the MRSs used in the manuscript use the output in directory "Results/MRS"  

df <- read.csv("Results/MRS/GWCOAD.csv", header=T) # adapt file name depending on which MRS you want to use 

# Make 1000 random combination between early- and later-onset categorization keeping the same numbers as in the discovery

B <- 1000
sample_matrix <- as.data.frame(replicate(B, sample(df$EOCRC, nrow(df), replace=F)))
rownames(sample_matrix)<- df$X
labels <- names(sample_matrix)

picsam <- cbind(df, sample_matrix)

picsam$gender[picsam$gender == "FEMALE"] <- 1
picsam$gender[picsam$gender == "MALE"] <- 0
picsam$gender <- as.factor(picsam$gender)

# Select the column names that have the permutated categorizations 
labels <- names(picsam[26:1025])

GLM <- matrix(ncol=3, nrow=length(labels))
colnames(GLM) <- c("Estimate", "Std. Error", "Pr(>|z|)")
row.names(GLM) <- labels

for (i in labels){
  newggpl <- picsam
  model <- paste(i , "~ PICLORAM + gender")
  test <- glm(model, family = binomial, data = newggpl) 
  results_df <-summary.glm(test)$coefficients
  results_df <- as.data.frame(t(results_df[2,]))
  GLM[i,1] <- results_df[1,1]
  GLM[i,2] <- results_df[1,2]
  GLM[i,3] <- results_df[1,4]
  
}

write.csv(GLM, "GLM_Patient_perm/picloram_coad_GW_perm.csv")

#
# The results obtained used in the manuscript figures are located in "Results/GLM_Patient_perm/"
#



#
#
# meta-analysis
#
# Data sets:
# Rectal: "READ", "GSE39958" 
# Colon: "GSE42752", "GSE131013"
# CRC: "READ", "GSE39958", "GSE42752", "GSE77954", "GSE101764", "GSE131013"
#
#

dir.create("meta_analysis")


source("Functions/Function_meta_analysis.R")



meta(method = "MRS", 
     PhenoCpG = "GW",                   #Options: "GW", "P1E5", "F01", "F005", "F001"
     datasets = c("READ", "GSE39958", "GSE42752", "GSE77954", "GSE101764", "GSE131013"),  
     type = "CRC"                       # This name will be added to the output file (no other function)
)


## Combine the COAD, READ, and the meta-analyses for the 5 thresholds per cancer type


files <- list.files("meta_analysis")
files2 <- files[grep("CRC", files)]  # CRC, Rectal, colon
tables <- lapply(files2, read.csv, header = TRUE)

ID <- str_remove(files2, "_MRS_CRC_") # CRC, Rectal, colon
ID <- str_remove(ID, "_.csv") 

tables <- mapply(cbind, tables, "SampleID"=ID, SIMPLIFY=F)
comb <- do.call(rbind , tables)
write.csv(comb, "meta_analysis/Comb_Meta_CRC.csv") # CRC, Rectal, colon

#
# The results obtained and used in the manuscript figures are located in "Results/meta_analysis/"
#


