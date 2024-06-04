Scripts to reproduce results of:


# Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  


Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 




author: Silvana C.E. Maas (silvanamaas at vhio.net) See methods section for more details.



---

# Folders: 

## scripts

#### 01_datamanagement.R 
Downloading, pre-processing, imputing, and horvath clock adjustment for the DNA methylation data 

#### 02_Methylation_Risk_score.R
Making the MRS scores and association between early- and later-onset, permutation for MRS CpGs and patients, meta-analysis

#### 03_Figure_1_and_SuplFigS3.R    
Script to make the Overview Figures presented in Figure 1 and supplementary Figure S3

#### 04_Figures_MRS_eoCRC.R
The heatmap, box plot, and forest plot as presented in Figure 2 and Figure 3a,b,c
The violinplot for the permutation results as presented in Figure 3d and supplementary Figure S4


#### 05_Figure_3_SBS1.R
Data management and download for the SBS1 analysis and the plots as presented in Figure 3e,f,g  

#### 06_SEER_Pesticide_Fig4.R          
Data management and download for the pesticide use data, the linear mixed models and the plots presented in Figure 4 

#### 07_Obesity_measured_vs_MRS_FigS1.R
Plots presented in Supplementary Figure S2 



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
