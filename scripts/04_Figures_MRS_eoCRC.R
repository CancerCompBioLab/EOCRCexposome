### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
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
                
                col.study = c("#5F9ED1FF", "#D94801", "#5F9ED1FF","#D94801", "#5F9ED1FF", "#5F9ED1FF", "#807DBA", "#5F9ED1FF", "#D94801", "#807DBA"),
                col.square = c("#5F9ED1FF", "#D94801", "#5F9ED1FF","#D94801", "#5F9ED1FF", "#5F9ED1FF", "#807DBA", "#5F9ED1FF", "#D94801", "#807DBA"),
                col.square.lines = "black")
    


     
#Figure 3f: permutation picloram
    
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

    
    
### Figure S6 permutation 

#
#
#   Permutation CpGs Figure S6a
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

### Fig S6a permutation CpGs

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
#   Permutation patients Figure S6b
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


### FigS6 permutation patients

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



#
#
# Figure 3d: Adjustment for other methylation risk scores- GW threshold
#
#

MRS <- read.csv("Results/MRS/GWCOAD.csv", header=T)
MRS$EOCRC <- as.factor(MRS$EOCRC)
trait <- names(MRS)[2:20]
trait <- trait[-15]

pesti <- matrix(ncol= 6, nrow= length(trait))
colnames(pesti) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CI_l", "CI_H")
row.names(pesti) <- trait

test <- glm(EOCRC ~ PICLORAM + gender,  family = binomial, data = MRS) 
summary.glm(test)$coefficients

for (i in trait){
  
  model <- paste0("EOCRC ~ scale(PICLORAM) + gender + ", i)
  
  test <- glm(model,  family = binomial, data = MRS) 
  results_df <-summary.glm(test)$coefficients
  results_df <- as.data.frame(t(results_df["scale(PICLORAM)",]))
  pesti[i,1] <- results_df[1,1]
  pesti[i,2] <- results_df[1,2]
  pesti[i,3] <- results_df[1,3]
  pesti[i,4] <- results_df[1,4]
  
  CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
  CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
  
  pesti[i,5] <- as.numeric(CIl)
  pesti[i,6] <- as.numeric(CIh)
  
}

pesti <- as.data.frame(pesti)

test <- glm(EOCRC ~ scale(PICLORAM) + gender,  family = binomial, data = MRS) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df["scale(PICLORAM)",]))
pesti["Basic",1] <- results_df[1,1]
pesti["Basic",2] <- results_df[1,2]
pesti["Basic",3] <- results_df[1,3]
pesti["Basic",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
pesti["Basic",5] <- as.numeric(CIl)
pesti["Basic",6] <- as.numeric(CIh)

pesti$logP <- -log10(pesti$`Pr(>|z|)`)

pesti$name <- rownames(pesti)
mete <- pesti
mete$name[mete$name == "A_2_4_D"] <- "2,4-D"
mete$name[mete$name == "ACETOCHLOR"] <- "Acetochlor"
mete$name[mete$name == "AHEI"] <- "AHEI"
mete$name[mete$name == "Alcohol"] <- "Alcohol"
mete$name[mete$name == "ATRAZINE"] <- "Atrazine"
mete$name[mete$name == "birthweight"] <- "Birthweight"
mete$name[mete$name == "BMI"] <- "BMI"
mete$name[mete$name == "CHLORDANE"] <- "Chlordane"
mete$name[mete$name == "coffee"] <- "Coffee"
mete$name[mete$name == "DDT"] <- "DDT"
mete$name[mete$name == "dicamba"] <- "Dicamba"
mete$name[mete$name == "education"] <- "Education"
mete$name[mete$name == "GLYPHOSATE"] <- "Glyphosate"
mete$name[mete$name == "HEPTACHLOR"] <- "Heptachlor"
mete$name[mete$name == "LINDANE"] <- "Lindane"
mete$name[mete$name == "MALATHION"] <- "Malathion"
mete$name[mete$name == "MDS"] <- "MDS"
mete$name[mete$name == "MESOTRIONE"] <- "Mesotrione"
mete$name[mete$name == "METOLACHLOR"] <- "Metolachlor"
mete$name[mete$name == "NO2"] <- "NO2"
mete$name[mete$name == "Obesity"] <- "Obesity"
mete$name[mete$name == "PCB"] <- "PCB" 
mete$name[mete$name == "PICLORAM"] <- "Picloram"
mete$name[mete$name == "pm10"] <- "PM10"
mete$name[mete$name == "pm2_5"] <- "PM2.5"
mete$name[mete$name == "pm2_510"] <- "PM2.5-10"
mete$name[mete$name == "smoking"] <- "Smoking"
mete$name[mete$name == "smoking_sm13"] <- "Smoking-Maas"
mete$name[mete$name == "TOXAPHENE"] <- "Toxaphene"

mete$name[mete$name == "Basic"] <- "Baseline"
mete$cat[mete$name ==  "2,4-D"] <- "Pesticides"
mete$cat[mete$name == "Acetochlor"] <- "Pesticides"
mete$cat[mete$name == "AHEI"] <- "Lifestyle"
mete$cat[mete$name == "Alcohol"] <- "Lifestyle"
mete$cat[mete$name == "Atrazine"] <- "Pesticides"
mete$cat[mete$name == "Birthweight"] <- "Lifestyle"
mete$cat[mete$name == "BMI"] <- "Lifestyle"
mete$cat[mete$name == "Chlordane"] <- "Pesticides"
mete$cat[mete$name == "Coffee"] <- "Lifestyle"
mete$cat[mete$name == "Delivery"] <- "Lifestyle"
mete$cat[mete$name == "DDT"] <- "Pesticides"
mete$cat[mete$name == "Dicamba"] <- "Pesticides"
mete$cat[mete$name == "Education"] <- "Lifestyle"
mete$cat[mete$name == "Glyphosate"] <- "Pesticides"
mete$cat[mete$name == "Heptachlor"] <- "Pesticides"
mete$cat[mete$name == "Lindane"] <- "Pesticides"
mete$cat[mete$name == "Malathion"] <- "Pesticides"
mete$cat[mete$name == "MDS"] <- "Lifestyle"
mete$cat[mete$name == "Mesotrione"] <- "Pesticides"
mete$cat[mete$name == "Metolachlor"] <- "Pesticides"
mete$cat[mete$name == "NO2"] <- "Air pollution"
mete$cat[mete$name == "Obesity"] <- "Lifestyle"
mete$cat[mete$name == "PCB"] <- "Air pollution"
mete$cat[mete$name == "Picloram"] <- "Pesticides"
mete$cat[mete$name == "PM10"] <- "Air pollution"
mete$cat[mete$name == "PM2.5"] <- "Air pollution"
mete$cat[mete$name == "PM2.5-10"] <- "Air pollution"
mete$cat[mete$name == "Smoking"] <- "Lifestyle"
mete$cat[mete$name == "Smoking-Maas"] <- "Lifestyle"
mete$cat[mete$name == "Toxaphene"] <- "Pesticides"
mete$cat[mete$name == "Baseline"] <- "Baseline"
pesti <- mete

library(dmetar)
library(tidyverse)
library(meta)

pesti$cat <- factor(pesti$cat, levels = c("Baseline", "Air pollution", "Lifestyle", "Pesticides"))
newggpl <- pesti

names(newggpl) <- c("beta", "se", "z.value", "pvalue", "CI_L","CI_H", "logP", "sample", "cat") 

m.gen_bin <- metagen(TE = beta,
                     seTE = se,
                     lower = CI_L,
                     upper = CI_H,
                     studlab = sample,
                     data = newggpl,
                     sm = "OR",
                     method.tau = "PM",
                     fixed = FALSE,
                     random = TRUE,
                     byvar= cat
)


### Figure 3d: Adjustments for other MRS-GW


forest.meta(m.gen_bin, 
            sortvar = studlab,
            prediction = F, 
            print.tau2 = F,
            leftcols = c("sample"),
            leftlabs = c("Adjustment"),
            hetstat = FALSE,
            random = F,
            rightcols = c("effect", "ci"),
            just.studlab =  c("left"),
            just.addcols =  c("left"),
            test.subgroup.random = F,
            print.subgroup.name = FALSE,
            col.by = c("black"),
            calcwidth.random=F,
            fs.heading= 5,
            fontsize =5,
            scientific.pval =T,
            plotwidth = "2.5cm",
            colgap.left = "0.15cm",
            colgap.forest.left="0.0001cm",
            colgap.forest.right="0.001cm",
            spacing=0.49,
            squaresize =0.8, 
            weight.study = "same" 
)


#
#
# Figure 3e: Adjustment for patients or tumor-specific markers
#
#

GLM <- matrix(ncol=7, nrow= 11)
colnames(GLM) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CI_l", "CI_H", "N")
row.names(GLM) <- c("Basic_M1", "subLR_M1", "subLR_M2", "history_M1","history_M2", "race_M1","race_M2", "MSS_M1", "MSS_M2", "Puritycon_M1", "Puritycon_M2")

# Download the full clinical data from TCGAbiolinks
clinic <- read.csv("/Startdata/COAD_clinic.csv", header=T)
clinic <- clinic [3:461,]
clinic <- clinic[,c("bcr_patient_barcode", "anatomic_neoplasm_subdivision", "family_history_colorectal_cancer")]
names(clinic) <- c("ID", "subdivision", "family_history")

GW <- read.csv("Results/MRS/GWCOAD.csv", header=T)

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = GW) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["Basic_M1",1] <- results_df[1,1]
GLM["Basic_M1",2] <- results_df[1,2]
GLM["Basic_M1",3] <- results_df[1,3]
GLM["Basic_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["Basic_M1",5] <- as.numeric(CIl)
GLM["Basic_M1",6] <- as.numeric(CIh)
GLM["Basic_M1",7] <- paste(table(GW$EOCRC)[2], table(GW$EOCRC)[1], sep= "/")

GWcl <- merge(GW, clinic, by.x="X", by.y="ID")

### Tumor purity

## Tumor purity was obtained from:
# Genetic Mechanisms of Immune Evasion in Colorectal Cancer 
# Catherine S. Grasso; Marios Giannakis; Daniel K. Wells; Tsuyoshi Hamada; Xinmeng Jasmine Mu; Michael Quist; Jonathan A. Nowak; Reiko Nishihara; Zhi Rong Qian; Kentaro Inamura;
# Teppei Morikawa; Katsuhiko Nosho; Gabriel Abril-Rodriguez; Charles Connolly; Helena Escuin-Ordinas; Milan S. Geybels; William M. Grady ;
# Li Hsu; Siwen Hu-Lieskovan; Jeroen R. Huyghe; Yeon Joo Kim; Paige Krystofinski; Mark D.M. Leiserson; Dennis J. Montoya; Brian B. Nadel; Matteo Pellegrini; Colin C. Pritchard; Cristina Puig-Saus; Elleanor H. Quist; Ben J. Raphael; Stephen J. Salipante; Daniel Sanghoon Shin; Eve Shinbrot; Brian Shirts; Sachet Shukla; Janet L. Stanford; Wei Sun; Jennifer Tsoi; Alexander Upfill-Brown; David A. Wheeler; Catherine J. Wu; Ming Yu; Syed H. Zaidi; Jesse M. Zaretsky; Stacey B. Gabriel; Eric S. Lander; Levi A. Garraway; Thomas J. Hudson; Charles S. Fuchs; Antoni Ribas; Shuji Ogino; Ulrike Peters
# Cancer Discov (2018) 8 (6): 730–749.
# https://doi.org/10.1158/2159-8290.CD-17-1327

sub <- GWcl

df <- readr::read_tsv("215982clinic_sup_grassso_cancerDisc18.tsv")
df <- as.data.frame(df)
subs <- merge(sub, df,by.y="Sample", by.x="X")

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = subs) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["Puritycon_M1",1] <- results_df[1,1]
GLM["Puritycon_M1",2] <- results_df[1,2]
GLM["Puritycon_M1",3] <- results_df[1,3]
GLM["Puritycon_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["Puritycon_M1",5] <- as.numeric(CIl)
GLM["Puritycon_M1",6] <- as.numeric(CIh)
GLM["Puritycon_M1",7] <- paste(table(subs$EOCRC)[2], table(subs$EOCRC)[1], sep= "/")

# Adjusted for purity as continous
test <- glm(EOCRC ~ scale(PICLORAM) + gender + TumorPurity, family = binomial, data = subs) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["Puritycon_M2",1] <- results_df[1,1]
GLM["Puritycon_M2",2] <- results_df[1,2]
GLM["Puritycon_M2",3] <- results_df[1,3]
GLM["Puritycon_M2",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["Puritycon_M2",5] <- as.numeric(CIl)
GLM["Puritycon_M2",6] <- as.numeric(CIh)
GLM["Puritycon_M2",7] <- paste(table(subs$EOCRC)[2], table(subs$EOCRC)[1], sep= "/")

####### Tumor location left/right

sub <- GWcl
sub <- sub[! sub$subdivision == "[Not Available]"  ,]
sub <- sub[! sub$subdivision == "[Discrepancy]"  ,]
sub$subdivision[sub$subdivision == "Cecum"] <- "Left"
sub$subdivision[sub$subdivision == "Ascending Colon"] <- "Left"
sub$subdivision[sub$subdivision == "Hepatic Flexure"] <- "Left"
sub$subdivision[sub$subdivision == "Transverse Colon"] <- "Left"

sub$subdivision[sub$subdivision == "Sigmoid Colon"] <- "Right"
sub$subdivision[sub$subdivision == "Descending Colon"] <- "Right"
sub$subdivision[sub$subdivision == "Rectosigmoid Junction"] <- "Right"
sub$subdivision[sub$subdivision == "Splenic Flexure"] <- "Right"
subs <- sub[complete.cases(sub$subdivision), ] 

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = subs) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["subLR_M1",1] <- results_df[1,1]
GLM["subLR_M1",2] <- results_df[1,2]
GLM["subLR_M1",3] <- results_df[1,3]
GLM["subLR_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["subLR_M1",5] <- as.numeric(CIl)
GLM["subLR_M1",6] <- as.numeric(CIh)
GLM["subLR_M1",7] <- paste(table(subs$EOCRC)[2], table(subs$EOCRC)[1], sep= "/")


test <- glm(EOCRC ~ scale(PICLORAM) + gender + subdivision, family = binomial, data = subs) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["subLR_M2",1] <- results_df[1,1]
GLM["subLR_M2",2] <- results_df[1,2]
GLM["subLR_M2",3] <- results_df[1,3]
GLM["subLR_M2",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["subLR_M2",5] <- as.numeric(CIl)
GLM["subLR_M2",6] <- as.numeric(CIh)
GLM["subLR_M2",7] <- paste(table(subs$EOCRC)[2], table(subs$EOCRC)[1], sep= "/")

## family_history yes/no

sub <- GWcl
sub <- sub[! sub$family_history == "[Not Available]"  ,]
sub <- sub[! sub$family_history == "[Unknown]"  ,]
sub$family_history[sub$family_history >0 ] <- 1

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["history_M1",1] <- results_df[1,1]
GLM["history_M1",2] <- results_df[1,2]
GLM["history_M1",3] <- results_df[1,3]
GLM["history_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["history_M1",5] <- as.numeric(CIl)
GLM["history_M1",6] <- as.numeric(CIh)
GLM["history_M1",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")


test <- glm(EOCRC ~ scale(PICLORAM) + gender + family_history, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["history_M2",1] <- results_df[1,1]
GLM["history_M2",2] <- results_df[1,2]
GLM["history_M2",3] <- results_df[1,3]
GLM["history_M2",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["history_M2",5] <- as.numeric(CIl)
GLM["history_M2",6] <- as.numeric(CIh)
GLM["history_M2",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")


## Race white/non-white 

sub <- GWcl
sub <- sub[! sub$race == "[Not Available]"  ,]
sub$race[sub$race != "WHITE"] <- "NO"

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["race_M1",1] <- results_df[1,1]
GLM["race_M1",2] <- results_df[1,2]
GLM["race_M1",3] <- results_df[1,3]
GLM["race_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["race_M1",5] <- as.numeric(CIl)
GLM["race_M1",6] <- as.numeric(CIh)
GLM["race_M1",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")


test <- glm(EOCRC ~ scale(PICLORAM) + gender + race, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["race_M2",1] <- results_df[1,1]
GLM["race_M2",2] <- results_df[1,2]
GLM["race_M2",3] <- results_df[1,3]
GLM["race_M2",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["race_M2",5] <- as.numeric(CIl)
GLM["race_M2",6] <- as.numeric(CIh)
GLM["race_M2",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")

### MSI status

# MSI status obtained from: 
# Genetic Mechanisms of Immune Evasion in Colorectal Cancer 
# Catherine S. Grasso; Marios Giannakis; Daniel K. Wells; Tsuyoshi Hamada; Xinmeng Jasmine Mu; Michael Quist; Jonathan A. Nowak; Reiko Nishihara; Zhi Rong Qian; Kentaro Inamura;
# Teppei Morikawa; Katsuhiko Nosho; Gabriel Abril-Rodriguez; Charles Connolly; Helena Escuin-Ordinas; Milan S. Geybels; William M. Grady ;
# Li Hsu; Siwen Hu-Lieskovan; Jeroen R. Huyghe; Yeon Joo Kim; Paige Krystofinski; Mark D.M. Leiserson; Dennis J. Montoya; Brian B. Nadel; Matteo Pellegrini; Colin C. Pritchard; Cristina Puig-Saus; Elleanor H. Quist; Ben J. Raphael; Stephen J. Salipante; Daniel Sanghoon Shin; Eve Shinbrot; Brian Shirts; Sachet Shukla; Janet L. Stanford; Wei Sun; Jennifer Tsoi; Alexander Upfill-Brown; David A. Wheeler; Catherine J. Wu; Ming Yu; Syed H. Zaidi; Jesse M. Zaretsky; Stacey B. Gabriel; Eric S. Lander; Levi A. Garraway; Thomas J. Hudson; Charles S. Fuchs; Antoni Ribas; Shuji Ogino; Ulrike Peters
# Cancer Discov (2018) 8 (6): 730–749.
# https://doi.org/10.1158/2159-8290.CD-17-1327


sub <- GWcl

df <- readr::read_tsv("215982clinic_sup_grassso_cancerDisc18.tsv")
df <- as.data.frame(df)  #[,c("Sample", "TumorPurity")])
subs <- merge(sub, df,by.y="Sample", by.x="X")
subs$MsiStatus[subs$MsiStatus == "POLE"] <- NA
sub <- subs[complete.cases(subs$MsiStatus)  ,]

test <- glm(EOCRC ~ scale(PICLORAM) + gender, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["MSS_M1",1] <- results_df[1,1]
GLM["MSS_M1",2] <- results_df[1,2]
GLM["MSS_M1",3] <- results_df[1,3]
GLM["MSS_M1",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["MSS_M1",5] <- as.numeric(CIl)
GLM["MSS_M1",6] <- as.numeric(CIh)
GLM["MSS_M1",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")


test <- glm(EOCRC ~ scale(PICLORAM) + gender + MsiStatus, family = binomial, data = sub) 
results_df <-summary.glm(test)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
GLM["MSS_M2",1] <- results_df[1,1]
GLM["MSS_M2",2] <- results_df[1,2]
GLM["MSS_M2",3] <- results_df[1,3]
GLM["MSS_M2",4] <- results_df[1,4]

CIl <-  confint(test)["scale(PICLORAM)","2.5 %"]
CIh <- confint(test)["scale(PICLORAM)","97.5 %"]
GLM["MSS_M2",5] <- as.numeric(CIl)
GLM["MSS_M2",6] <- as.numeric(CIh)
GLM["MSS_M2",7] <- paste(table(sub$EOCRC)[2], table(sub$EOCRC)[1], sep= "/")


comb <- as.data.frame(GLM)
names(comb) <- c("beta", "se", "z.value", "pvalue", "CI_L","CI_H", "sample") 
comb <- comb[complete.cases(comb),]
comb$studynum <- rownames(comb)
comb$mod <- gsub(".*_","",comb$studynum)
comb$studynum <- gsub("\\_.*","",comb$studynum)

comb$studynum[comb$studynum == "subLR"] <- "Location"
comb$studynum[comb$studynum == "history"] <- "Family history"
comb$studynum[comb$studynum == "race"] <- "Race"
comb$studynum[comb$studynum == "MSS"] <- "MSS/MSI-H"
comb$studynum[comb$studynum == "Puritycon"] <- "Tumor Purity"

newggpl <- comb
newggpl$studynum[newggpl$studynum == "Basic"] <- "Baseline"
newggpl <- newggpl %>%
  mutate_at(vars(beta, se, z.value, pvalue, CI_L, CI_H), as.numeric)

m.gen_bin <- metagen(TE = beta,
                     seTE = se,
                     lower = CI_L,
                     upper = CI_H,
                     studlab = studynum,
                     data = newggpl,
                     sm = "OR",
                     method.tau = "PM",
                     fixed = FALSE,
                     random = TRUE#,
)

# Figure 3e: Adjustment association between piclom and age at onset in COAD 

forest.meta(m.gen_bin, 
            sortvar = studynum,
            prediction = F, 
            print.tau2 = F,
            leftcols = c("studynum", "mod","sample"),
            leftlabs = c("Adjustment", "Model", "Early/Late"),
            hetstat = FALSE,
            random = F,
            rightcols = c("effect", "ci"),
            just.studlab =  c("left"),
            just.addcols =  c("left"),
            test.subgroup.random = F,
            print.subgroup.name = FALSE,
            col.by = c("black"),
            calcwidth.random=F,
            fs.heading= 5,
            fontsize =5,
            scientific.pval =T,
            plotwidth = "2.5cm",
            colgap.left = "0.15cm",
            colgap.forest.left="0.0001cm",
            colgap.forest.right="0.001cm",
            spacing=0.49,
            squaresize =0.8, 
            weight.study = "same", 
            col.study = c("#802A07FF", "#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da"),
            col.square = c("#802A07FF", "#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da","#FFC08AFF", "#7eb1da"), 
            col.square.lines = "black")









