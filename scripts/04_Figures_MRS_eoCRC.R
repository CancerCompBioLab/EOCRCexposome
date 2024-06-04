### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)





library(ComplexHeatmap)
library(ggpubr)
library(ggplot2)

#
# Figure 2 per trait
# Figure 3 for A, B, and C
#
#

# Input: Matrix with CpGs in rows
met_mat_M  <- readRDS("Startdata/COAD_Age_adjusted_DNAm_Bval_overlap_MRSCpGs.rds")
met_mat_M[1:4,1:4]

clinic <- read.csv("Startdata/Clinic_Coad.csv", header=T)
head(clinic)

clinic$EOCRC = as.factor(ifelse(clinic$age < 50 , 1,
                                ifelse(clinic$age >=70, 0, NA)))

clinic <- clinic[complete.cases(clinic$EOCRC), ]
met_mat_M <- met_mat_M[, colnames(met_mat_M) %in% clinic$ID]
clinic <- clinic[clinic$ID %in% colnames(met_mat_M), ]


phenoGW = read.csv(file="InputData/CombGW.csv", header=TRUE)
head(phenoGW)
phenoGW <- phenoGW[,c(1,2,4,6)]
colnames(phenoGW) = c("CpG", "Beta","P",  "Phenotype")

# Change "education" for the trait of interest

pico <- phenoGW[phenoGW$Phenotype == "education", ]
pico <- pico[complete.cases(pico),]

# Change "education" for the trait of interest
MRS = read.csv(file="Results/MRS/GWCOAD.csv", header=TRUE)
GSVA <- MRS[,c("X", "EOCRC", "gender",  "education")]
rownames(GSVA) <- GSVA$X

DNAm <- as.data.frame(t(met_mat_M))
DNAm2 <- DNAm[, colnames(DNAm) %in% pico$CpG]


DNA <- merge(DNAm2, GSVA, by=0)
DNA <- DNA[order(DNA$education, decreasing =T), ]
pico <- pico[order(pico$Beta, decreasing =T),]
DNA2<-DNA[pico$CpG]

col_fun = colorRamp2(c(-4, -1.5, 1.5, 4), c("#FFF5EB", "#FFC08AFF", "#D94801", "#662105")) 

# Add a row annotation for the methylation risk score
row_ha = rowAnnotation(MRS = DNA$education, col = list(MRS = col_fun), show_legend = F, border = T,
                       simple_anno_size = unit(0.25, "cm"),
                       annotation_name_gp= gpar(fontsize = 5),
                       annotation_legend_param = list(MRS = list(title = "MRS"))) 


# Add the barplot above the heatmap to show the direction in the original EWAS
col_fun2 = colorRamp2(c(-0.03,-0.020,0,0.025,0.04), c("#061746","#5F9ED1FF","#F7FBFF","#D94801","#802A07FF"))
column_ha = HeatmapAnnotation(Direction =pico$Beta, col = list(Direction = col_fun2), show_legend = F, border = T,
                              simple_anno_size = unit(0.25, "cm"),
                              annotation_name_gp= gpar(fontsize = 5),
                              annotation_legend_param = list(Direction = list(title = "Direction")))

col_fun3 = colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#9ECAE1",  "#061746"))

col_fun22 = colorRamp2(c(-0.005,-0.0025,0,0.0025,0.005), c("#061746","#5F9ED1FF","#F7FBFF","#D94801","#802A07FF"))
lgd1 = Legend(at = c(-0.005,0.005), labels= c("-", "+"), title = "Direction", 
              col_fun = col_fun22, border = "black", direction = "vertical",
              labels_gp = gpar(fontsize = 5),
              legend_height = unit(1, "cm"),
              grid_width = unit(0.25, "cm"),
              title_gap = unit(1, "mm"),
              title_gp = gpar(fontsize = 5, fontface = "plain"))

lgd2 = Legend(at = c(-4,  3), labels= c("-4", "3"),title = "MRS", 
              col_fun = col_fun, border = "black", direction = "vertical",
              labels_gp = gpar(fontsize = 5),
              grid_width = unit(0.25, "cm"),
              legend_height = unit(1, "cm"),
              title_gap = unit(1, "mm"),
              title_gp = gpar(fontsize = 5, fontface = "plain")) 

lgd3 = Legend(col_fun = col_fun3, title = "Beta",  border = "black", direction = "vertical", 
              labels_gp = gpar(fontsize = 5),
              legend_height = unit(1, "cm"),
              grid_width = unit(0.25, "cm"),
              title_gap = unit(1, "mm"),
              title_gp = gpar(fontsize = 5, fontface = "plain"), at = c(0,0.5, 1), labels= c("0","0.5", "1"))

pd = packLegend(lgd1,  lgd3, lgd2, direction = "vertical", gap = unit(2, "mm"))

map <- Heatmap(data.matrix(DNA2), 
               right_annotation = row_ha, top_annotation = column_ha, col= col_fun3,
               show_row_names = F, show_row_dend = F, show_column_dend = FALSE, show_column_names = T, column_names_rot = 80,
               column_names_gp = gpar(fontsize = 5), show_parent_dend_line = FALSE, cluster_rows = F, cluster_columns = F,
               show_heatmap_legend = F, heatmap_legend_param = list(title= c("Beta")), row_dend_reorder = F, column_dend_reorder = F, 
               
)

# Figure 2 and 3 left panel: Heatmap
draw(map, heatmap_legend_list = list(pd), heatmap_legend_side = "left", annotation_legend_side = "left")



##Figure 2 and 3 middle panel: boxplot

# Load the methylation risk scores, Change "GW" for the threshold of interest
MRS = read.csv(file="Results/MRS/GWCOAD.csv", header=TRUE)

# Change "education" for the trait of interest
DNA <- MRS[,c("X", "EOCRC", "gender",  "education")]
DNA$EOCRC <- as.factor(DNA$EOCRC)
library(ggplot2)

ggplot(DNA, aes(x=EOCRC, y=education)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = gender), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Sex:",
                    values = c("#FDAE6B", "#98A7D4"), 
                    labels= c("Female", "Male"))+
  ylim(-4,3)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))+  
  ylab("Methylation Risk Score") +
  scale_x_discrete(labels=c('Later-onset', 'Early-onset'))



#
#
# META-ANALYSIS
#
#

# Change GW for the threshold and select the data sets to be included


datasets <- c("Results/GLM_MRS/GW_COAD_gender.csv", "Results/GLM_MRS/GW_GSE101764_gender.csv", "Results/GLM_MRS/GW_GSE131013NoAdjust.csv",
           "Results/GLM_MRS/GW_GSE39958_gender.csv", "Results/GLM_MRS/GW_GSE42752_gender.csv", "Results/GLM_MRS/GW_GSE77954NoAdjust.csv",                     
           "Results/GLM_MRS/GW_READ_gender.csv")   

tables <- lapply(datasets, read.csv, header = TRUE)
  comb <- do.call(rbind , tables)
  comb$N <- NULL
  names(comb) <- c("beta", "se", "z.value", "pvalue", "CI_L","CI_H", "study", "name") 
  comb <- comb[complete.cases(comb),]

  comb$EOCRC <- ifelse(comb$study == "COAD", "31/100",
                       ifelse(comb$study ==  "GSE101764",  "13/37",
                              ifelse(comb$study ==  "GSE131013", "2/58",
                                     ifelse(comb$study =="GSE39958", "12/9",
                                            ifelse(comb$study =="GSE42752", "4/7",
                                                   ifelse(comb$study =="GSE77954", "3/10", 
                                                          ifelse(comb$study =="READ", "14/30", NA_character_ )))))))
  library(dmetar)
  library(meta)
 
  # make the separation between TCGA-COAD discovery and the other datasets as replication
   comb$Tissue <-  ifelse(comb$study != "COAD", "Replication", "Discovery") 
  
    i <- "education"   # Any trait of interest 
    newggpl <- comb
    newggpl <- newggpl[which(newggpl$name == i), ]
    
    title <- paste(i, "OR",sep="_")
    
    m.gen_bin <- metagen(TE = beta,
                         seTE = se,
                         lower = CI_L,
                         upper = CI_H,
                         studlab = EOCRC,
                         data = newggpl,
                         sm = "OR",
                         method.tau = "PM",
                         fixed = FALSE,
                         random = TRUE,
                         title = title,
                         byvar= Tissue)
    
    
### Figure 2 and 3 right panel: meta-analysis Forest plot
    
    
      forest.meta(m.gen_bin, 
                sortvar = studlab,
                prediction = F, 
                print.tau2 = F,
                leftcols = c("study", "EOCRC"),
                leftlabs = c("Study", "Early/Later"),
                hetstat = FALSE,
                random = TRUE,
                rightcols = c("effect", "ci"),
                just.studlab =  c("left"),
                just.addcols =  c("left"),
                test.subgroup.random = F,
                print.subgroup.name = FALSE,
                col.by = c("black"),
                text.random = "Overall: Discovery and Replication",
                text.random.w = "Overall: Replication",
                calcwidth.random=F,
                fs.heading= 5,
                fontsize =5,
                scientific.pval =T,
                plotwidth = "2.5cm",
                colgap.left = "0.001cm",
                colgap.forest.left="0.001cm",
                colgap.forest.right="0.001cm",
                spacing=0.49,
                squaresize =1.1, 
                
                col.study = c("#5F9ED1FF", "#D94801", "#5F9ED1FF", "#807DBA", "#5F9ED1FF", "#D94801", "#807DBA"),
                col.square = c("#5F9ED1FF", "#D94801", "#5F9ED1FF", "#807DBA", "#5F9ED1FF", "#D94801", "#807DBA"), 
                col.square.lines = "black")
    


     
#Figure 3d: permutation picloram
    
PICLORAM <- read.csv("Results/GLM_MRS_perm/perm_GW_COAD_genderPICLORAM.csv", header=T)
PICLORAM <- PICLORAM[,c(7,4)]
names(PICLORAM) <- c("Label", "Pval")
PICLORAM$trait <- "CpGs"
    
library(dplyr)
    
PICLORAM <- as.data.frame(PICLORAM %>% 
                                mutate(rank = rank(Pval, ties.method = "first")))
    
picloram <- read.csv("Results/GLM_Patient_perm/picloram_coad_GW_perm.csv", header=T)
head(picloram)
picloram <- picloram[,c(1,4)]
names(picloram) <- c("Label", "Pval")
picloram$trait <- "Patients"
    
#Select the results obtained in the original categorization
GW_COAD <- read.csv("Results/GLM_MRS/GW_COAD_gender.csv", header=T)
signGW <- GW_COAD[GW_COAD$name == "PICLORAM",]
signGW$Label <- "MRS"
signGW <- signGW[,c(10,4)]
names(signGW) <- c("Label", "Pval")
signGW$trait <- "Patients"
    
comb <- rbind(picloram, signGW)
comb <- as.data.frame(comb %>% 
                            mutate(rank = rank(Pval, ties.method = "first")))
    
PICLORAM <- rbind(comb, PICLORAM)
    
library(ggplot2)


PICLORAM$Pval <- -log10(PICLORAM$Pval)
sigt <- PICLORAM[PICLORAM$Label == "MRS",]
    
ggplot(PICLORAM, aes(x= trait, y = Pval, fill= trait)) +
      geom_violin(linewidth = 0.2) +
      geom_boxplot(width=0.1, size=0.2, outlier.size = 0.1)+
      theme_bw()+
      scale_fill_manual(values=c(Patients= "#98A7D4", CpGs= "#98A7D4"), guide= "none")+
      theme(axis.title.x=element_blank(), text = element_text(size = 5),
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
      xlab("Methylation risk score Trait") + ylab("-log10 (P)")+
      geom_point(data=sigt, aes(x= trait, y = Pval, colour= 'Sepal.Width', shape= 'Sepal.Width'), size=2, stroke= 0.5)+ 
      geom_text(data=sigt, aes(label=Label), hjust=-0.5, vjust=1.0, size=1, colour= "maroon")+
      scale_color_manual(name="Sample", values = c("Sepal.Width" = "maroon"), labels = c("Model"), guide= "none")+
      scale_shape_manual(name="Sample", values = c("Sepal.Width" = 8), labels = c("Model"), guide= "none")+
      ylim(0,6)

    
    
### Figure S4 permutation 

#
#
#   Permutation CpGs Figure S4a
#
#



education <- read.csv("Results/GLM_MRS_perm/perm_GW_COAD_gendereducation.csv", header=T)
education <- education[,c(7,4)]
names(education) <- c("Label", "Pval")
education$trait <- "education"

MDS <- read.csv("Results/GLM_MRS_perm/perm_GW_COAD_genderMDS.csv", header=T)
MDS <- MDS[,c(7,4)]
names(MDS) <- c("Label", "Pval")
MDS$trait <- "MDS"


Obesity <- read.csv("Results/GLM_MRS_perm/perm_P1E5_COAD_genderObesity.csv", header=T)
Obesity <- Obesity[,c(7,4)]
names(Obesity) <- c("Label", "Pval")
Obesity$trait <- "Obesity"

smoking_sm13 <- read.csv("Results/GLM_MRS_perm/perm_GW_COAD_gendersmoking_sm13.csv", header=T)
smoking_sm13 <- smoking_sm13[,c(7,4)]
names(smoking_sm13) <- c("Label", "Pval")
smoking_sm13$trait <- "smoking_sm13"

comb <- rbind(education, MDS, Obesity, smoking_sm13)

comb <- as.data.frame(comb %>% arrange(trait, Pval) %>% 
                        group_by(trait) %>%
                        mutate(rank = rank(Pval, ties.method = "first")))


library(ggplot2)

comb$Pval <- -log10(comb$Pval)
sigt <- comb[comb$Label == "MRS",]

### Fig S4a permutation CpGs

ggplot(comb, aes(x= factor(trait, levels = c("MDS",  "Obesity", "education","smoking_sm13")), y = Pval, fill= trait)) +
  geom_violin(linewidth = 0.2) +
  geom_boxplot(width=0.1, size=0.2, outlier.size = 0.1)+
  theme_bw()+
  scale_fill_manual(values=c(MDS= "#DDE1F1",Obesity=    "#98A7D4", education= "#748CC6", smoking_sm13= "#4C72B7"), guide= "none")+
  theme(axis.title.x=element_blank(),
        text = element_text(size = 5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
  xlab("Methylation risk score traits") + ylab("-log10 (P)")+
  geom_point(data=sigt, aes(x= trait, y = Pval, colour= 'Sepal.Width', shape= 'Sepal.Width'), size=2, stroke= 0.5)+   
  geom_text(data=sigt, aes(label=Label), hjust=-0.5, vjust=1.0, size=1, colour= "maroon")+
  scale_color_manual(name="Sample", values = c("Sepal.Width" = "maroon"), labels = c("Model"), guide= "none")+
  scale_shape_manual(name="Sample", values = c("Sepal.Width" = 8), labels = c("Model"), guide= "none")+
  scale_x_discrete(labels=c("MDS \n 3 CpGs", "Obesity \n 11 CpGs", "Education \n 13 CpGs", "Smoking-Maas    \n 11 CpGs"))+
  ylim(0,6)

#
#
#   Permutation patients Figure S4b
#
#


education <- read.csv("Results/GLM_Patient_perm/education_coad_GW_perm.csv", header=T)
education <- education[,c(1,4)]
names(education) <- c("Label", "Pval")
education$trait <- "education"

MDS <- read.csv("Results/GLM_Patient_perm/MDS_coad_GW_perm.csv", header=T)
MDS <- MDS[,c(1,4)]
names(MDS) <- c("Label", "Pval")
MDS$trait <- "MDS"

Obesity <- read.csv("Results/GLM_Patient_perm/Obesity_coad_P1E5_perm.csv", header=T)
Obesity <- Obesity[,c(1,4)]
names(Obesity) <- c("Label", "Pval")
Obesity$trait <- "Obesity"

smoking_sm13 <- read.csv("Results/GLM_Patient_perm/smoking_sm13_coad_GW_perm.csv", header=T)
smoking_sm13 <- smoking_sm13[,c(1,4)]
names(smoking_sm13) <- c("Label", "Pval")
smoking_sm13$trait <- "smoking_sm13"

traits1 <- c("education", "MDS", "smoking_sm13")

traits2 <- c("Obesity")

GW_COAD <- read.csv("Results/GLM_MRS/GW_COAD_gender.csv", header=T)
signGW <- GW_COAD[GW_COAD$name %in% traits1,]
signGW$Label <- "Model"
signGW <- signGW[,c(10,4,9)]
names(signGW) <- c("Label", "Pval", "trait")


P1E5_COAD <- read.csv("Results/GLM_MRS/P1E5_COAD_gender.csv", header=T)
signP1E5 <- P1E5_COAD[P1E5_COAD$name %in% traits2,]
signP1E5$Label <- "Model"
signP1E5 <- signP1E5[,c(10,4,9)]
names(signP1E5) <- c("Label", "Pval", "trait")

comb <- rbind(education, MDS, Obesity, smoking_sm13, signGW, signP1E5)
comb <- as.data.frame(comb %>% arrange(trait, Pval) %>% 
                        group_by(trait) %>%
                        mutate(rank = rank(Pval, ties.method = "first")))

comb$Pval <- -log10(comb$Pval)
sigt <- comb[comb$Label == "Model",]


### FigS3 permutation patients

ggplot(comb, aes(x= factor(trait, levels = c("MDS",  "Obesity", "education","smoking_sm13")), y = Pval, fill= trait)) +
  geom_violin(linewidth = 0.2) +
  geom_boxplot(width=0.1, size=0.2, outlier.size = 0.1)+
  theme_bw()+
  scale_fill_manual(values=c(MDS= "#DDE1F1",Obesity=    "#98A7D4", education= "#748CC6", smoking_sm13= "#4C72B7"), guide= "none")+
  theme(axis.title.x=element_blank(), text = element_text(size = 5),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
  xlab("Methylation risk score traits") + ylab("-log10 (P)")+
  geom_point(data=sigt, aes(x= trait, y = Pval, colour= 'Sepal.Width', shape= 'Sepal.Width'), size=2, stroke= 0.5)+   
  geom_text(data=sigt, aes(label=Label), hjust=-0.5, vjust=1.0, size=1,colour= "maroon")+
  scale_color_manual(name="Sample", values = c("Sepal.Width" = "maroon"), labels = c("Model"), guide= "none")+
  scale_shape_manual(name="Sample", values = c("Sepal.Width" = 8), labels = c("Model"), guide= "none")+
  scale_x_discrete(labels=c("MDS \n 3 CpGs", "Obesity \n 11 CpGs", "Education \n 13 CpGs", "Smoking-Maas    \n 11 CpGs"))+
  ylim(0,6)












