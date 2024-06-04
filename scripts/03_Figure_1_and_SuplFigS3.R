### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)






library(stringr)
library(ComplexHeatmap)
library(tidyr)
library(circlize)


# In this script we use the results obtained as described in the manuscript.

# If you followed "scripts/02_Methylation_Risk_score.R" and want to use your own obtained results, 
# please remove "Results/" from the file directory names

COAD_mrs <-read.csv("Results/meta_analysis/Comb_COAD.csv", header=T)
COAD_mrs$trait <- paste("MRS_coad", COAD_mrs$analysis, sep="_")

# set non-significant traits to 0
COAD_mrs$esi <- ifelse(COAD_mrs$Pr...z.. > 0.05, 0, COAD_mrs$Estimate)
COAD_mrs <- COAD_mrs[, c("esi", "trait", "name")]
names(COAD_mrs) <- c("Esti", "variable", "name") 


CRC_mrs <-read.csv("Results/meta_analysis/Comb_Meta_CRC.csv", header=T)
CRC_mrs$trait <- paste("CRC_mrs", CRC_mrs$SampleID, sep="_")
# set non-significant traits to 0
CRC_mrs$esi <- ifelse(CRC_mrs$pval > 0.05, 0, CRC_mrs$Beta)
CRC_mrs <- CRC_mrs[, c("esi", "trait", "name")]
names(CRC_mrs) <-  c("Esti", "variable", "name") 


new <- rbind(COAD_mrs, CRC_mrs)
mete <- spread(new, variable, Esti)

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


mete$cat <- NA

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

data_wide <- mete
rownames(data_wide) <- data_wide$name
data_wide$name <- NULL


# assign the direction of the significant associations 
data_wide[data_wide < 0] <- -1
data_wide[data_wide > 0] <- 1
data_wide[is.na(data_wide)] <- -20

text_list <- mete$cat
data_wide$cat <- NULL
names(data_wide)
names <- c("CRC F001", "CRC F005", "CRC F01", "CRC GW", "CRC P1E5", 
           "COAD F001", "COAD F005", "COAD F01", "COAD GW", "COAD P1E5") 
names(data_wide)<- names
data_wide <- data_wide[,c("COAD GW", "COAD P1E5", "COAD F01", "COAD F005", "COAD F001", 
                          "CRC GW", "CRC P1E5", "CRC F01", "CRC F005", "CRC F001")] 


met <- c("TCGA-COAD", "TCGA-COAD", "TCGA-COAD", "TCGA-COAD", "TCGA-COAD", 
         "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC")
natis <- c("GW", "P1E5","F01", "F005", "F001", "GW", "P1E5", "F01", "F005", "F001")

names(data_wide) <- natis
over <- data_wide

# Give the colors to Not available, not significant, negative and positive results
col_fun = colorRamp2(c(-20, 0, -1, 1), c("#fbfbfb", "#f2eeee" , "#98A7D4", "#FDAE6B"))

ha = rowAnnotation(foo = anno_block(gp = gpar(fill=rep(c("#FEE6CE",  "#FEE6CE", "#FEE6CE")), col="white"),
                                    labels = c("Air pollution", "Lifestyle", "Pesticides"),
                                    labels_gp = gpar(col = "black", fontsize=7)))

na = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=rep(c("#FEE6CE",  "#FEE6CE")), col="white"),
                                        labels = c("TCGA-COAD", "Meta-analysis CRC"),
                                        labels_gp = gpar(col = "black", fontsize=7)))


met2 <- sort(factor(met, levels = c("TCGA-COAD", "Meta-analysis CRC")))

# Make the custom legend
lgd1 = Legend(labels= c("Positive", "None", "Negative", "N.A"), labels_gp = gpar(fontsize = 5), 
                            title = "Association to early-onset colon and rectal cancer", nrow=1,
              title_gp = gpar(fontsize = 6), legend_gp = gpar(fill = c("#FDAE6B","#f2eeee","#98A7D4", "#fbfbfb")), 
              border = "black")


# Figure 1

map <-   Heatmap(data.matrix(over), 
         show_row_names = T, show_row_dend = FALSE, show_column_dend = FALSE, show_column_names = T,
         show_parent_dend_line = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,  row_title=NULL, 
         col = col_fun,  rect_gp = gpar(col = "white", lwd = 0.5), show_heatmap_legend = F,
         row_dend_reorder = F, column_dend_reorder = F, left_annotation = ha, top_annotation = na, row_split = text_list, row_names_gp = gpar(fontsize = 7),     
         column_split = met2, column_gap=unit(.0125, "npc"), column_names_gp = gpar(fontsize = 7), column_names_rot = 45,
         column_title_gp = gpar(fontsize = 7),  column_title = "Marker selection thresholds", column_title_side = "bottom")

draw(map, padding = unit(c(15, 2, 2, 2), "mm"))
draw(lgd1, x = unit(1, "cm"), y = unit(0.15, "cm"), just = c("left", "bottom"))


#
#  Supplementary Figure S3
#


library(stringr)
library(ComplexHeatmap)

COAD <-read.csv("Results/meta_analysis/Comb_COAD.csv", header=T)
COAD$trait <- paste("MRS_coad", COAD$analysis, sep="_")
COAD$esi <- ifelse(COAD$Pr...z.. > 0.05, 0, COAD$Estimate)
COAD <- COAD[, c("esi", "trait", "name")]
names(COAD) <- c("Esti", "variable", "name") 

READ <- read.csv("Results/meta_analysis/Comb_READ.csv", header=T)
READ$trait <- paste("MRS_READ", READ$SampleID, sep="_")
READ$esi <- ifelse(READ$pvalue > 0.05, 0, READ$beta)
READ <- READ[, c("esi", "trait", "name")]
names(READ) <- c("Esti", "variable", "name") 


rectal <-read.csv("Results/meta_analysis/Comb_Meta_Rectal.csv", header=T)
rectal$trait <- paste("MRS_rectal", rectal$SampleID, sep="_")
rectal$esi <- ifelse(rectal$pval > 0.05, 0, rectal$Beta)
rectal <- rectal[, c("esi", "trait", "name")]
names(rectal) <- c("Esti", "variable", "name") 

Colon <-read.csv("Results/meta_analysis/Comb_Meta_Colon.csv", header=T)
Colon$trait <- paste("MRS_colon", Colon$SampleID, sep="_")
Colon$esi <- ifelse(Colon$pval > 0.05, 0, Colon$Beta)
Colon <- Colon[, c("esi", "trait", "name")]
names(Colon) <- c("Esti", "variable", "name") 

CRC <-read.csv("Results/meta_analysis/Comb_Meta_CRC.csv", header=T)
CRC$trait <- paste("MRS_CRC", CRC$SampleID, sep="_")
CRC$esi <- ifelse(CRC$pval > 0.05, 0, CRC$Beta)
CRC <- CRC[, c("esi", "trait", "name")]
names(CRC) <- c("Esti", "variable", "name") 

new <- rbind(COAD,READ,rectal,Colon, CRC)

library(tidyr)

mete <- spread(new, variable, Esti)

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

mete$cat <- NA

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

data_wide <- mete
rownames(data_wide) <- data_wide$name
data_wide$name <- NULL

data_wide[data_wide < 0] <- -1
data_wide[data_wide > 0] <- 1
data_wide[is.na(data_wide)] <- -20

text_list <- mete$cat
data_wide$cat <- NULL

names(data_wide)
names <- c("COAD F001",   "COAD F005", "COAD F01", "COAD GW", "COAD P1E5",
           "Colon F001",  "Colon F005",  "Colon F01",   "Colon GW",    "Colon P1E5",  
           "CRC F001", "CRC F005", "CRC F01", "CRC GW", "CRC P1E5", 
           "READ F001", "READ F005", "READ F01", "READ GW", "READ P1E5", 
           "Rectal F001", "Rectal F005", "Rectal F01", "Rectal GW", "Rectal P1E5") 

names(data_wide)<- names

data_wide <- data_wide[,c("COAD GW", "COAD P1E5", "COAD F01", "COAD F005", "COAD F001", 
                          "READ GW", "READ P1E5", "READ F01", "READ F005", "READ F001", 
                          "Rectal GW", "Rectal P1E5", "Rectal F01", "Rectal F005", "Rectal F001", 
                          "Colon GW", "Colon P1E5", "Colon F01", "Colon F005", "Colon F001", 
                          "CRC GW", "CRC P1E5", "CRC F01", "CRC F005", "CRC F001")] 

met <- c("TCGA-COAD", "TCGA-COAD", "TCGA-COAD", "TCGA-COAD", "TCGA-COAD", 
         "TCGA-READ", "TCGA-READ", "TCGA-READ", "TCGA-READ", "TCGA-READ", 
         "Meta-analysis Rectal", "Meta-analysis Rectal", "Meta-analysis Rectal", "Meta-analysis Rectal", "Meta-analysis Rectal",
         "Meta-analysis Colon", "Meta-analysis Colon", "Meta-analysis Colon", "Meta-analysis Colon", "Meta-analysis Colon",
         "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC", "Meta-analysis CRC")

natis <- c("GW", "P1E5","F01", "F005", "F001", 
           "GW", "P1E5","F01", "F005", "F001", 
           "GW", "P1E5","F01", "F005", "F001", 
           "GW", "P1E5","F01", "F005", "F001", 
           "GW", "P1E5", "F01", "F005", "F001")

names(data_wide) <- natis
over <- data_wide

col_fun = colorRamp2(c(-20, 0, -1, 1), c("#fbfbfb", "#f2eeee" , "#98A7D4", "#FDAE6B"))

ha = rowAnnotation(foo = anno_block(gp = gpar(fill=rep(c("#FEE6CE",  "#FEE6CE", "#FEE6CE")), col="white"),
                                    labels = c("Air pollution", "Lifestyle", "Pesticides"),
                                    labels_gp = gpar(col = "black", fontsize=7)))


na = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=rep(c("#FEE6CE",  "#FEE6CE")), col="white"),
                                        labels = c("TCGA-COAD",
                                                   "TCGA-READ",
                                                   "Meta-analysis Rectal",
                                                   "Meta-analysis Colon", 
                                                   "Meta-analysis CRC"),
                                        labels_gp = gpar(col = "black", fontsize=7)))

met2 <- sort(factor(met, levels = c("TCGA-COAD",
                                    "TCGA-READ",
                                    "Meta-analysis Rectal",
                                    "Meta-analysis Colon", 
                                    "Meta-analysis CRC")))

lgd1 = Legend(labels= c("Positive", "None", "Negative", "N.A"), labels_gp = gpar(fontsize = 5), 
              title = "Association to early-onset colon and rectal cancer", nrow=1,
              title_gp = gpar(fontsize = 6), legend_gp = gpar(fill = c("#FDAE6B","#f2eeee","#98A7D4", "#fbfbfb")), 
              border = "black")

map <- 
  Heatmap(data.matrix(over), 
          show_row_names = T, 
          show_row_dend = FALSE, 
          show_column_dend = FALSE, 
          show_column_names = T,
          show_parent_dend_line = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,  row_title=NULL, 
          col = col_fun,  rect_gp = gpar(col = "white", lwd = 0.5),
          show_heatmap_legend = F,
          row_dend_reorder = F, column_dend_reorder = F, left_annotation = ha, top_annotation = na,
          row_split = text_list, 
          row_names_gp = gpar(fontsize = 7),  
          column_split = met2,
          column_gap=unit(.0125, "npc"),
          column_names_gp = gpar(fontsize = 5),
          column_names_rot = 45,
          column_title = "Marker selection thresholds", 
          column_title_side = "bottom",
          column_title_gp = gpar(fontsize = 7),)


draw(map, padding = unit(c(15, 2, 2, 2), "mm"))

draw(lgd1, x = unit(1, "cm"), y = unit(0.15, "cm"), just = c("left", "bottom"))

