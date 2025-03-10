lib <- "/home/smaas/R/x86_64-pc-linux-gnu-library/4.1"
.libPaths(lib)
library(maftools)
library(fabricatr)

clinic <- read.csv("Startdata/COAD_clinic.csv", header=T)
clinic <- clinic [3:461,]

clinic <- clinic[,c("bcr_patient_barcode", "ethnicity",
                      "tissue_source_site", "anatomic_neoplasm_subdivision", "family_history_colorectal_cancer",
                      "weight_kg_at_diagnosis", "height_cm_at_diagnosis",
                      "ajcc_staging_edition","ajcc_tumor_pathologic_pt","ajcc_nodes_pathologic_pn","ajcc_pathologic_tumor_stage")]

names(clinic) <- c("ID", "ethnicity", "TSS", "subdivision", "family_history", "weight", "Height", "staging", "pt", "pn", "stage")

GW <- read.csv("Results/MRS/GW_COAD_full.csv", header=T)

GW$onset = as.factor(ifelse(GW$age < 50 , "eoCRC",
                              ifelse(GW$age >= 50 & GW$age < 70, "moCRC",
                                     ifelse(GW$age >=70, "loCRC", NA))))
GWcl <- merge(GW, clinic, by.x="X", by.y="ID")


# Download tumor characteristic data from:
# Genetic Mechanisms of Immune Evasion in Colorectal Cancer 
# Catherine S. Grasso; Marios Giannakis; Daniel K. Wells; Tsuyoshi Hamada; Xinmeng Jasmine Mu; Michael Quist; Jonathan A. Nowak; Reiko Nishihara; Zhi Rong Qian; Kentaro Inamura;
# Teppei Morikawa; Katsuhiko Nosho; Gabriel Abril-Rodriguez; Charles Connolly; Helena Escuin-Ordinas; Milan S. Geybels; William M. Grady ;
# Li Hsu; Siwen Hu-Lieskovan; Jeroen R. Huyghe; Yeon Joo Kim; Paige Krystofinski; Mark D.M. Leiserson; Dennis J. Montoya; Brian B. Nadel; Matteo Pellegrini; Colin C. Pritchard; Cristina Puig-Saus; Elleanor H. Quist; Ben J. Raphael; Stephen J. Salipante; Daniel Sanghoon Shin; Eve Shinbrot; Brian Shirts; Sachet Shukla; Janet L. Stanford; Wei Sun; Jennifer Tsoi; Alexander Upfill-Brown; David A. Wheeler; Catherine J. Wu; Ming Yu; Syed H. Zaidi; Jesse M. Zaretsky; Stacey B. Gabriel; Eric S. Lander; Levi A. Garraway; Thomas J. Hudson; Charles S. Fuchs; Antoni Ribas; Shuji Ogino; Ulrike Peters
# Cancer Discov (2018) 8 (6): 730–749.
# https://doi.org/10.1158/2159-8290.CD-17-1327

sub <- GWcl
df <- readr::read_tsv("215982clinic_sup_grassso_cancerDisc18.tsv")
df <- as.data.frame(df)

GWcl2 <- merge(GWcl, df,by.y="Sample", by.x="X", all.x=T)
GWcl2$EOCRC = as.factor(ifelse(GWcl2$age < 50 , 1,
                                ifelse(GWcl2$age >=70, 0, NA)))

GWcl2$onset <- factor(GWcl2$onset, levels=c("loCRC", "moCRC", "eoCRC"))
names(GWcl2)[1] <- "SampleID"
colnames(GWcl2)[which(names(GWcl2) == "Tumor Long Barcode / Filename")] <- "Tumor_Sample_Barcode"
GWcl2$Tumor_Sample_Barcode <- str_remove(GWcl2$Tumor_Sample_Barcode, "_Illumina")
GWcl2$Tumor_Sample_Barcode <- str_remove(GWcl2$Tumor_Sample_Barcode, "_hg19")
GWcl2 <- GWcl2[substr(GWcl2$Tumor_Sample_Barcode, 14,15) == "01", ]

GWcl2 <- GWcl2 %>%
mutate(piclo_Terc = split_quantile(PICLORAM, 3))

GWcl2$piclo_Terc <- as.numeric(GWcl2$piclo_Terc)

GWcl2$piclo_Terc[GWcl2$piclo_Terc == "1"] <- "Low"
GWcl2$piclo_Terc[GWcl2$piclo_Terc == "2"] <- "Mid"
GWcl2$piclo_Terc[GWcl2$piclo_Terc == "3"] <- "High"

GWcl2$piclo_Terc <- factor(GWcl2$piclo_Terc, levels = c("Low", "Mid", "High"))

#write.csv(GWcl2, "Startdata/clinic_GW.csv", row.names = F)

GWcl3 <- GWcl2[complete.cases(GWcl2$MsiStatus), ]
GWcl3 <- GWcl3[GWcl3$MsiStatus == "MSS",]

GWcl3 <- GWcl3 %>%
mutate(piclo_Terc = split_quantile(PICLORAM, 3))

#write.csv(GWcl3, "Startdata/clinic_GW_MSS.csv", row.names = F)


GWcl4 <- GWcl2[complete.cases(GWcl2$MsiStatus), ]
GWcl4 <- GWcl4[GWcl4$MsiStatus == "MSI-H",]
GWcl4 <- GWcl4 %>%
mutate(piclo_Terc = split_quantile(PICLORAM, 3))

#write.csv(GWcl4, "Startdata/clinic_GW_MSI.csv", row.names = F)

## Download mutation data using TCGAbiolinks

maf = "COAD_TCGA_snv.maf" 
clin = "Startdata/clinic_GW.csv"

library(maftools)
laml = read.maf(maf = maf, clinicalData = clin)

clinic <-read.csv("Startdata/clinic_GW.csv", header=T)

### change mutation categories

laml@data$Variant_Classification <- as.character(laml@data$Variant_Classification)
laml@data$Variant_Classification[laml@data$Variant_Classification == "Frame_Shift_Del"] <- "Truncating"        
laml@data$Variant_Classification[laml@data$Variant_Classification == "Frame_Shift_Ins"] <- "Truncating"           
laml@data$Variant_Classification[laml@data$Variant_Classification == "In_Frame_Del"] <- "In Frame"           
laml@data$Variant_Classification[laml@data$Variant_Classification == "In_Frame_Ins"] <- "In Frame"      
laml@data$Variant_Classification[laml@data$Variant_Classification == "Missense_Mutation"] <- "Missense"      
laml@data$Variant_Classification[laml@data$Variant_Classification == "Nonsense_Mutation"] <- "Missense" 
laml@data$Variant_Classification[laml@data$Variant_Classification == "Nonstop_Mutation"] <- "Truncating"            
laml@data$Variant_Classification[laml@data$Variant_Classification == "Splice_Site"] <- "Truncating" 
laml@data$Variant_Classification[laml@data$Variant_Classification == "Translation_Start_Site"] <- "Truncating" 
laml@data$Variant_Classification <- as.factor(laml@data$Variant_Classification)

clin <- subsetMaf(maf = laml, tsb = clinic$Tumor_Sample_Barcode)
colnames(clin@clinical.data)[which(names(clin@clinical.data) == "piclo_Terc")] <- "Picloram"

highclin <-  clinic[clinic$piclo_Terc == "High",]
lowclin <-  clinic[clinic$piclo_Terc == "Low",]

highpicl <- subsetMaf(maf = clin, tsb = highclin$Tumor_Sample_Barcode)
lowpicl <- subsetMaf(maf = clin, tsb = lowclin$Tumor_Sample_Barcode)

pt.vs.rt <- mafCompare(m1 = lowpicl, m2 = highpicl, m1Name = 'Low', m2Name = 'High', minMut = 0)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05)
genes2 <- pt.vs.rt$results
#write.csv(genes2, "Results/Molecular_differences/genesCompletedata.csv",row.names = F)

genessigCom <- pt.vs.rt$results[pt.vs.rt$results$pval <0.05,]
genesCom = genessigCom$Hugo_Symbol
genesComr <- rev(genesCom)

#Color coding 
fabcolors = c("#FDAE6B", "#D94801","#802A07FF")
names(fabcolors) = c("Low", "Mid", "High")
fabcolors = list(Picloram = fabcolors)

#vc_cols = mutations
vc_cols = c("#bed9ee","#5f9ed1", "#335b8c","#122858")
names(vc_cols) = c( "In Frame", "Missense", "Truncating", "Multi_Hit")

## Supplementary Figure S7a
coBarplot(m1 = lowpicl, m2 = highpicl, m1Name = 'Low', m2Name = 'High', genes = genesCom, colors=vc_cols,
          geneSize = 0.5, pctSize =0.5, axisSize = 0.5, legendTxtSize =0.5, titleSize= 0.25)

clinic <-   clinic[order(clinic$PICLORAM),]
order <-   clinic$Tumor_Sample_Barcode

## Supplementary Figure S7b

oncoplot(maf = clin, genes = genesComr, removeNonMutated = F,
         bgCol = "#EBEBEB", colors = vc_cols, 
         clinicalFeatures =  c("Picloram"),
         sortByAnnotation = T, annotationColor = fabcolors, 
         drawRowBar = F,
         SampleNamefontSize = 0.5, sepwd_genes = 0.25, sepwd_samples = 0.6, 
         fontSize = 0.5, titleFontSize = 0.5, legendFontSize = 0.5, annotationFontSize = 0.5, anno_height = 0.5, keepGeneOrder = T,
         sampleOrder = order) 

####################

######### MSI ##############

MSIclinic <-read.csv("Startdata/clinic_GW_MSI.csv", header=T)

highMSI <-  MSIclinic[MSIclinic$piclo_Terc == "High",]
lowMSI <-  MSIclinic[MSIclinic$piclo_Terc == "Low",]

MSI <- subsetMaf(maf = laml, tsb = MSIclinic$Tumor_Sample_Barcode)
colnames(MSI@clinical.data)[which(names(MSI@clinical.data) == "piclo_Terc")] <- "Picloram"

highMSIpicl <- subsetMaf(maf = MSI, tsb = highMSI$Tumor_Sample_Barcode)
lowMSIpicl <- subsetMaf(maf = MSI, tsb = lowMSI$Tumor_Sample_Barcode)

pt.vs.rt <- mafCompare(m1 = lowMSIpicl, m2 = highMSIpicl, m1Name = 'Low', m2Name = 'High', minMut = 0)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05)
genesMSI <- pt.vs.rt$results
#write.csv(genesMSI, "Results/Molecular_differences/genesMSI_ALLgenes.csv",row.names = F)

genessigMSI <- pt.vs.rt$results[pt.vs.rt$results$pval <0.05,]
genesMSI = genessigMSI$Hugo_Symbol
genesMSIr <- rev(genesMSI)

#Color coding 
fabcolors = c("#FDAE6B", "#D94801","#802A07FF")
names(fabcolors) = c("Low", "Mid", "High")
fabcolors = list(Picloram = fabcolors)

vc_cols = c("#bed9ee","#5f9ed1", "#335b8c","#122858")
names(vc_cols) = c( "In Frame", "Missense", "Truncating", "Multi_Hit")

## Figure 4. Molecular differences between high and low Picloram MRS-GW

## Figure 4c
coBarplot(m1 = lowMSIpicl, m2 = highMSIpicl, m1Name = 'MSI-Low', m2Name = 'MSI-High', genes = genesMSIr, colors=vc_cols,
          geneSize = 0.5, pctSize =0.5, axisSize = 0.5, legendTxtSize =0.5, titleSize= 0.25)

MSIclinic <-   MSIclinic[order(MSIclinic$PICLORAM),]
MSIorder <-   MSIclinic$Tumor_Sample_Barcode

## Figure 4d
oncoplot(maf = MSI, genes = genesMSI, removeNonMutated = F,
         bgCol = "#EBEBEB", colors = vc_cols, 
         clinicalFeatures =  c("Picloram"),
         sortByAnnotation = T, annotationColor = fabcolors, 
         drawRowBar = F,
         SampleNamefontSize = 0.5, sepwd_genes = 0.25, sepwd_samples = 0.6, 
         fontSize = 0.5, titleFontSize = 0.5, legendFontSize = 0.5, annotationFontSize = 0.5, anno_height = 0.5, keepGeneOrder = T,
         sampleOrder = MSIorder) 


##### MSS ####
MSSclinic <-read.csv("Startdata/clinic_GW_MSS.csv", header=T)
highMSS <-  MSSclinic[MSSclinic$piclo_Terc == "High",]
lowMSS <-  MSSclinic[MSSclinic$piclo_Terc == "Low",]

MSS <- subsetMaf(maf = laml, tsb = MSSclinic$Tumor_Sample_Barcode)
colnames(MSS@clinical.data)[which(names(MSS@clinical.data) == "piclo_Terc")] <- "Picloram"

highMSSpicl <- subsetMaf(maf = MSS, tsb = highMSS$Tumor_Sample_Barcode)
lowMSSpicl <- subsetMaf(maf = MSS, tsb = lowMSS$Tumor_Sample_Barcode)

pt.vs.rt <- mafCompare(m1 = lowMSSpicl, m2 = highMSSpicl, m1Name = 'Low', m2Name = 'High', minMut = 0)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05)
genesMSS <- pt.vs.rt$results
#write.csv(genesMSS, "Results/Molecular_differences/genesMSS_ALL.csv",row.names = F)

genessigMSS <- pt.vs.rt$results[pt.vs.rt$results$pval <0.05,]
genesMSS = genessigMSS$Hugo_Symbol
genesMSSr <- rev(genesMSS)

## Figure 4a 

coBarplot(m1 = lowMSSpicl, m2 = highMSSpicl, m1Name = 'MSS-Low', m2Name = 'MSS-High', genes = genesMSSr, colors=vc_cols,
          geneSize = 0.5, pctSize =0.5, axisSize = 0.5, legendTxtSize =0.5, titleSize= 0.25)

MSSclinic <-   MSSclinic[order(MSSclinic$PICLORAM),]
MSSorder <-   MSSclinic$Tumor_Sample_Barcode

## Figure 4b

oncoplot(maf = MSS, genes = genesMSS, removeNonMutated = F,
         bgCol = "#EBEBEB", colors = vc_cols, clinicalFeatures =  c("Picloram"),
         sortByAnnotation = T, annotationColor = fabcolors, drawRowBar = F,
         SampleNamefontSize = 0.5, sepwd_genes = 0.25, sepwd_samples = 0.6, 
         fontSize = 0.5, titleFontSize = 0.5, legendFontSize = 0.5, annotationFontSize = 0.5, anno_height = 0.5, keepGeneOrder = T,
         sampleOrder = MSSorder) 

#
# DEA high vs low picloram adjusted for tumour purity, age, MSI status, sex and tumour site
#
############

## Download TCGA gene expression data using TCGAbiolinks
# Load gene expression TCGA COAD
cts_COAD<-readRDS("COAD_TCGA_gene-expression.RDS")

#Check type of gene expression data
cts_COAD<-assays(cts_COAD)$unstranded   #Need to use raw counts for DESEQ2

# Only keep the ones with tumour, which ends with 01A = tumor 
cts_COAD<-as.data.frame(cts_COAD)
cts_COAD <- dplyr::select(cts_COAD,contains("01A"))

#Shorten sample name to 12 characters in cts
colnames(cts_COAD)<-substr(colnames(cts_COAD), 1, 12)
colData <- read.csv("Startdata/clinic_GW.csv", header=T)


# Make a new picloram tertile for MSI only
colData_MSI<-colData[colData$MsiStatus=="MSI-H",]
quantile(colData_MSI$PICLORAM, c(.33, .66, .99))

# > quantile(colData$PICLORAM, c(.33, .66, .99))
# 33%         66%         99% 
# 0.001066336 0.013617596 0.030819614 

# Add categorical variable to the data frame
colData_MSI$PICLORAM_tertile_MSI <- as.factor(ifelse(colData_MSI$PICLORAM<(0.001066336),'low',
                                                     ifelse(colData_MSI$PICLORAM<0.013617596, 'medium', 'high')))
# Remove medium Picloram part (not needed for DEA)
colData_MSI<-colData_MSI[colData_MSI$PICLORAM_tertile_MSI=="low"|colData_MSI$PICLORAM_tertile_MSI=="high",]

#Make first column to rownames
row.names(colData_MSI) <- colData_MSI$bcr_patient_barcode
colData_MSI[1] <- NULL

# Keep only samples with gene expression and colData_MSI

#Remove columns in cts which are not rows in colData_MSI
keep<-rownames(colData_MSI)
cts_COAD_MSI<-cts_COAD[,colnames(cts_COAD) %in% keep]

#Remove columns in colData_MSI which are not rows in cts_COAD_MSI
keep<-colnames(cts_COAD_MSI)
colData_MSI<-colData_MSI[rownames(colData_MSI) %in% keep,]

#Order alphabetically
colData_MSI<-colData_MSI[order(row.names(colData_MSI)), ]
cts_COAD_MSI<-cts_COAD_MSI[,order(colnames(cts_COAD_MSI)), ]

#############################################################

# Make new Picloram Tertile for MSS only
colData_MSS<-colData[colData$MsiStatus=="MSS",]

#Make a new variable Picloram low vs. high
quantile(colData_MSS$PICLORAM, c(.33, .66, .99))

# > quantile(colData_MSS$PICLORAM, c(.33, .66, .99))
# 33%          66%          99% 
# -0.002667754  0.005046180  0.026088225 

# Add categorical variable to the data frame
colData_MSS$PICLORAM_tertile_MSS <- as.factor(ifelse(colData_MSS$PICLORAM<(-0.002667754),'low',
                                                     ifelse(colData_MSS$PICLORAM<0.005046180, 'medium', 'high')))

# Remove medium Picloram part (not needed for DEA)
colData_MSS<-colData_MSS[colData_MSS$PICLORAM_tertile_MSS=="low"|colData_MSS$PICLORAM_tertile_MSS=="high",]

#Make first column to rownames
row.names(colData_MSS) <- colData_MSS$bcr_patient_barcode
colData_MSS[1] <- NULL

# Keep only samples with gene expression and colData_MSS

#Remove columns in cts which are not rows in colData_MSS
keep<-rownames(colData_MSS)
cts_COAD_MSS<-cts_COAD[,colnames(cts_COAD) %in% keep]

#Remove columns in colData_MSS which are not rows in cts_COAD_MSS
keep<-colnames(cts_COAD_MSS)
colData_MSS<-colData_MSS[rownames(colData_MSS) %in% keep,]

#Order alphabetically
colData_MSS<-colData_MSS[order(row.names(colData_MSS)), ]
cts_COAD_MSS<-cts_COAD_MSS[,order(colnames(cts_COAD_MSS)), ]

##################################

# Remove medium Picloram part (not needed for DEA)
colData<-colData[colData$PICLORAM_tertile=="low"|colData$PICLORAM_tertile=="high",]

#Make first column to rownames
row.names(colData) <- colData$bcr_patient_barcode
colData[1] <- NULL

# Keep only samples with gene expression and colData

#Remove columns in cts which are not rows in colData
keep<-rownames(colData)
cts_COAD<-cts_COAD[,colnames(cts_COAD) %in% keep]

#Remove columns in colData which are not rows in cts_COAD
keep<-colnames(cts_COAD)
colData<-colData[rownames(colData) %in% keep,]

#Order alphabetically
colData<-colData[order(row.names(colData)), ]
cts_COAD<-cts_COAD[,order(colnames(cts_COAD)), ]

############

dds0 <- DESeqDataSetFromMatrix(countData = cts_COAD,   # cts_COAD_MSS or cts_COAD_MSI
                               colData = colData,     # colData_MSS or colData_MSI  
                               design = ~ age_at_initial_pathologic_diagnosis + MsiStatus + TumorPurity + gender + Site + PICLORAM_tertile)

# in MSI tumors no adjustment for "Site" (left/right) was included due to limited sample size 


smallestGroupSize <- 66 #Needs to be adjusted to smallest group size
keep <- rowSums(counts(dds0) >= 10) >= smallestGroupSize
dds0 <- dds0[keep,]

dds0 = DESeq(dds0)
res<-results(dds0,contrast=c("PICLORAM_tertile","low","high"))
res@rownames=gsub("\\..*","",res@rownames)
res@listData$ENSEMBL<-res@rownames
res<-as.data.frame(res@listData)

#Annotate res
ens.str <- res$ENSEMBL
mart <- useMart('ENSEMBL_MART_ENSEMBL', host='www.ensembl.org')
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <-  getBM(mart = mart, attributes = c('ensembl_gene_id', 'hgnc_symbol','chromosome_name'), filter = 'ensembl_gene_id', values = ens.str, uniqueRows = TRUE)
res<-merge(res,annotLookup,by.x="ENSEMBL",by.y="ensembl_gene_id")

#Remove NA rows
res<-na.omit(res)

# save(res, file ="Results/Molecular_differences/DEA_comp_mss_msi.csv")
# save(res, file ="Results/Molecular_differences/DEA_comp_mss.csv")
# save(res, file ="Results/Molecular_differences/DEA_comp_msi.csv")

### Volcano plot  

res <- read.csv("Results/Molecular_differences/DEA_comp_mss_msi.csv", header = T)
sub <- res[order(res$pvalue),]

## Figure 4e
plot <- EnhancedVolcano(res,
                        lab = res$hgnc_symbol,
                        selectLab =sub$hgnc_symbol[1:8],
                        x = 'log2FoldChange',
                        y = 'padj',
                        ylab = bquote(~-Log[10] ~ italic(Padj)),
                        ylim= c(0,5),
                        xlim = c(-4.5,4.5),
                        axisLabSize = 5,
                        titleLabSize = 1,
                        subtitleLabSize = 1,
                        labSize = 2.0,
                        legendLabSize = 5,
                        captionLabSize=5,
                        legendIconSize = 2.5,
                        col=c("#FDAE6B", "#98A7D4","green","red"),
                        legendLabels = c("Not Significant", "Log2FC > Threshold", "Adj. P-Value < Threshold", "Significant"),
                        title = 'Differential Expression Analysis',
                        subtitle = "low vs. high Picloram with adjustm. (incl. tumour purity)",
                        pCutoff =  0.05,
                        FCcutoff = 1.2,
                        pointSize = 0.6,
                        drawConnectors=F,
                        max.overlaps=Inf,
                        legendPosition="none",
                        cutoffLineType = "longdash",
                        cutoffLineCol = "#3C3C3C",
                        cutoffLineWidth = 0.25,
                        gridlines.major = F,
                        gridlines.minor = F,
                        borderWidth = 0.4
)


plot + theme(panel.grid.major = element_line(color = "#CECECE",
                                             size = 0.5,
                                             linetype = 1),
             panel.grid.minor = element_line(color = "#DADADA",
                                             size = 0.25,
                                             linetype = 1),
             axis.ticks = element_line(size=0.25))



############## FGSEA

res <- read.csv("Results/Molecular_differences/DEA_comp_mss_msi.csv", header = T)

res<- res[order(-res$stat),]
ranks<-res$stat
names(ranks)<-res$hgnc_symbol

#Prepare pathways
# Load the pathways into a named list
download.file("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs//h.all.v2024.1.Hs.symbols.gmt",destfile ="h.all.v2024.1.Hs.symbols.gmt")
ln = readLines("h.all.v2024.1.Hs.symbols.gmt")
ln = strsplit(ln, "\t")
gs = lapply(ln, function(x) x[-(1:2)])
names(gs) = sapply(ln, function(x) x[1])
pathways.hallmark<-gs

#Run FGSEA
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

#Filter out non-significant adj. p-value + NES <1:
fgseaRes<-fgseaRes[abs(fgseaRes$NES)>1,]
fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]

#Only show top up & down-regulated pathways:
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval)), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval)), pathway]

#Make a nice plot 
topPathways <- c(topPathwaysUp,topPathwaysDown)
fgseaResTidy<-fgseaResTidy[fgseaResTidy$pathway %in% topPathways,]

#Shorten name
fgseaResTidy$pathway<- gsub("^.*?_","_",fgseaResTidy$pathway)
fgseaResTidy$pathway=str_replace_all(fgseaResTidy$pathway, "_", " ")
fgseaResTidy <- as.data.frame(fgseaResTidy)
fgseaResTidy$enr <- ifelse(fgseaResTidy$NES >0, "Up", "Down")

## Figure 4f

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill= factor(enr, levels= c("Up", "Down"))))+
  geom_col() +
  scale_fill_manual(name = "In low exposure:",values = c("#FDAE6B", "#5f9ed1")) + 
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))

