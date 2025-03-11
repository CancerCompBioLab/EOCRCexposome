Scripts to reproduce results of:


# Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  


Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 




author: Silvana C.E. Maas (silvanamaas at vhio.net) See methods section for more details.



---

# Folders: 

## scripts

#### 01_datamanagement.R 
Downloading, pre-processing, imputing, and horvath clock adjustment for the DNA methylation data 

#### 02_Methylation_Risk_score.R
Making the MRS scores and association between early- and later-onset, permutation for MRS CpGs and patients, meta-analysis

#### 03_Figure_1_and_SuplFigS4.R    
Script to make the Overview Figures presented in Figure 1 and supplementary Figure S4

#### 04_Figures_MRS_eoCRC.R
The heatmap, box plot, and forest plot as presented in Figure 2 and Figure 3a,b,c
The violinplot for the permutation results as presented in Figure 3f and supplementary Figure S6
The forestplots for the additional adjustments as presented in Figure 3d and 3e.


#### 05_Figure_5_SBS1.R
Data management and download for the SBS1 analysis and the plots as presented in Figure 5.  

#### 06_SEER_Pesticide_Fig6.R          
Data management and download for the pesticide use data, the linear mixed models and the plots presented in Figure 6 

#### 07_Score_validation_FS1_FS2.R
Plots presented in Supplementary Figure S3 and S5 



## Functions
Contains the functions that are being called in the scripts in folder "scripts"


## InputData
All the information needed to generate the same start data sets used in the manuscript. 

#### CpGs included for each marker selection threshold, including their weights, p-values, and FDR. <br>
- CombGW.csv <br>
- CombP1E5.csv <br>
- CombF01.csv <br>
- CombF005.csv <br>
- CombF001.csv <br>            

#### The patient samples included after data management in TCGA-COAD and TCGA-READ <br>
- Included_samples_Patients_COAD.csv <br>
- Included_samples_Patients_READ.csv <br>

#### The CpGs that overlapped within all the available data sets, after exclusion of CpGs that are SNPs and cross-reactive probes. <br>
- overlap_CpGs.csv                  <br>


## Results
All the output files generated in ***02_Methylation_Risk_score.R*** that are used to generate all the figures included in the manuscript. 

## Startdata

#### Clinical information for the patients included in COAD <br>
COAD_clinic.csv       <br>
Clinic_Coad.csv <br>
clinic_GW_MSI.csv <br>
clinic_GW_MSS.csv <br>
clinic_GW.csv <br>

#### Age adjusted DNA methylation for the patients included in COAD <br>
COAD_Age_adjusted_DNAm_Bval_overlap_MRSCpGs.rds <br>
                         
#### Methylation risk scores- GW generated using non-age adjusted DNA methylation data in COAD <br>
COAD_GW_MRS_notAdj.csv <br>

#### Tumor mutation burden in COAD, using 38 as standard exon size <br>
TMB_per38.csv             <br>                     
<br>

## R packages and software <br>
TCGABiolinks package v2.25.0 <br>
GEOquery package v2.60.0 <br>
sesame package v1.19.7 <br>
impute package v1.66.0 <br>
lumi package v2.44.0  <br>
stats R package v4.1.3 <br>
metafor R package v4.4-0  <br>
maftools R package v2.8.5  <br>
DESeq2 R package v1.44.0  <br>
fgsea R package v1.20.0  <br>
lmerTest R package v3.1-3  <br>

Statistical analyses were performed in R v4.1.3  <br>
SEER*Stat software (version 8.4.1)  <br>
 <br>