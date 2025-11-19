### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)



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

clinic$height[clinic$height == "80.3"] <- "180.3"
clinic[clinic== "[Not Available]"] <- NA
clinic$weight <- as.numeric(clinic$weight)
clinic$height <- as.numeric(clinic$height)
clinic$height <- clinic$height/100
clinic$BMI <- clinic$weight/(clinic$height*clinic$height)
clinic$age <- as.numeric(clinic$age)
clinComp <- clinic[complete.cases(clinic$BMI),]

clinCompAge <- clinComp[complete.cases(clinComp$EOCRC),]
clinCompAge$obe <- ifelse(clinCompAge$BMI >30, 1,0)

MRS = read.csv(file="Results/MRS/P1E5COAD.csv", header=TRUE)
GSVA <- MRS[,c("X", "Obesity")]
rownames(GSVA) <- GSVA$X
GSVA <- GSVA[order(row.names(GSVA)), ]
cliob <- merge(clinCompAge, GSVA, by.x= "ID", by.y = "X")
cliob$obe <- as.factor(cliob$obe)



# Figure S2

ggplot(cliob, aes(x=EOCRC, y=Obesity)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = obe), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Obesity",
                    values = c("#FFC08AFF", "#9ECAE1"), 
                    labels= c("No", "Yes"))+
  ylim(-4,2.5)+
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

library(grid)

RR <- "RR: 0.67(0.26-1.76)"
grobbar <- grobTree(textGrob(RR, x=0.25,  y=0.75, hjust=0,
                          gp=gpar(col="black", fontsize=5, fontface="italic")))

# Figure S1a
ggplot(data=cliob, aes(x=EOCRC, fill= obe)) +
  geom_bar(aes(fill = obe))+
  scale_fill_manual(name = "Obesity",
                    values = c("#FFC08AFF", "#9ECAE1"), 
                    labels= c("No", "Yes"))+
  
  geom_text(stat='count', aes(label=after_stat(count)), vjust=1.75,size=1.75)+
  annotation_custom(grobbar)+
  scale_x_discrete(labels=c('Later-onset', 'Early-onset'))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))






# Smoking Maas MRS 
# Validated in TCGA data sets with additional smoking information
# TCGA data was extracted and processed following scripts/01_datamanagement.R
# the MRS was constructed following scripts/02_Methylation_Risk_score.R
# Additional clinical information was obtained from https://gdc.cancer.gov/about-data/publications/pancanatlas


ESCA <- read.csv("Results/MRS/GW_ESCA.csv")
CESC <- read.csv("Results/MRS/GW_CESC.csv")
LUSC <- read.csv("Results/MRS/GW_LUSC.csv")
PAAD <- read.csv("Results/MRS/GW_PAAD.csv")
LUAD <- read.csv("Results/MRS/GW_LUAD.csv")
HNSC <- read.csv("Results/MRS/GW_HNSC.csv")

GW <- rbind(ESCA, CESC, LUSC, PAAD, LUAD, HNSC)
#Download data from https://gdc.cancer.gov/about-data/publications/pancanatlas
clinpan <- read.csv("clinical_PANCAN_patient_with_followup.tsv", header=T, sep="\t")
dat <- clinpan[, c(2,3,4,10, 31, 56,57, 78, 79, 80, 84, 89, 206, 207, 211, 508, 519, 621, 622, 744)]
dat <- dat %>% mutate_all(na_if,"")
dat$tobacco_smoking_history[dat$tobacco_smoking_history == "[Not Available]"] <- NA
dat$tobacco_smoking_history[dat$tobacco_smoking_history == "[Unknown]"] <- NA
dat$tobacco_smoking_history[dat$tobacco_smoking_history == "[Discrepancy]"] <- NA

hc <- merge(GW, dat, by.x="X", by.y = "bcr_patient_barcode")
hc$smokhist <- hc$tobacco_smoking_history

hc$smokhist[hc$smokhist == "Current Reformed Smoker, Duration Not Specified"] <- NA
hc$smokhist[hc$smokhist == "Current reformed smoker for < or = 15 years"] <- "Former <=15y"
hc$smokhist[hc$smokhist == "Current reformed smoker for > 15 years"] <- "Former >15y"
hc$smokhist[hc$smokhist == "Lifelong Non-smoker"] <- "Never smoker"

hc <- hc[complete.cases(hc$smokhist),]

# Table S4: Smoking-Maas MRS validation number of patients per cancer type and smoking categorization in TCGA
table(hc$smokhist, hc$acronym)

hc$smokhist2 <- hc$tobacco_smoking_history

hc$smokhist2[hc$smokhist2 == "Current Reformed Smoker, Duration Not Specified"] <- NA
hc$smokhist2[hc$smokhist2 == "Current reformed smoker for < or = 15 years"] <- "Former smoker"
hc$smokhist2[hc$smokhist2 == "Current reformed smoker for > 15 years"] <- "Former smoker"
hc$smokhist2[hc$smokhist2 == "Lifelong Non-smoker"] <- "Never smoker"

hc$smokhist3 <- as.factor(ifelse(hc$smokhist2 == "Current smoker", "Ever",
                                 ifelse(hc$smokhist2 == "Never smoker", "Never", 
                                        ifelse(hc$smokhist2 == "Former smoker",  "Ever", NA))))

#Supplementary Figure S1: The relation between recorded smoking and alcohol consumption habits and the methylation risk scores  
#Supplementary Figure S1d

ggplot(data = hc, aes(x=smokhist, y=scale(smoking_sm13)))+
  geom_boxplot(notch = T, outlier.size = 0.5)+ 
  ylim(-4.5,8)+
  annotate("text",
           x = 1:length(table(hc$smokhist)),
           y = -4.5,
           label = table(hc$smokhist),
           col = "red", size = 1.5)+ 
  stat_compare_means(comparisons = list(c(1, 2),c(1, 3),c(1, 4),
                                        c(2, 3), c(2, 4), c(3, 4)), size = 1.5,
                     step.increase = 0.085,
                     tip.length = 0.02)+
  theme_bw()+
  xlab("Measured smoking habits ") + 
  ylab("Smoking-Maas MRS")+ theme(legend.position = "bottom",
                                  axis.text=element_text(size=5),
                                  legend.text=element_text(size=5),
                                  legend.title = element_text(size=5),
                                  axis.title = element_text(size = 5),
                                  legend.margin=margin(0,0,0,0),
                                  legend.box.margin=margin(-10,0,-5,0),
                                  legend.key.size = unit(0.25, 'cm'))

#Supplementary Figure S1c
ggplot(data = hc, aes(x=smokhist2, y=scale(smoking_sm13)))+
  geom_boxplot(notch = T, outlier.size = 0.5)+
  ylim(-4.5,8)+
  annotate("text",
           x = 1:length(table(hc$smokhist2)),
           y = -4.5,
           label = table(hc$smokhist2),
           col = "red", size = 1.5)+ 
  stat_compare_means(comparisons = list(c(1, 2),c(1, 3), c(2, 3)),size = 1.5,
                     step.increase = 0.085,
                     tip.length = 0.02)+
  theme_bw()+
  xlab("Measured smoking habits ") + 
  ylab("Smoking-Maas MRS") +
  theme(legend.position = "bottom",
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))




# Smoking- Maas MRS 
# Validated in blood samples from GSE50660
# Data was extracted and processed following scripts/01_datamanagement.R
# The MRS was constructed following scripts/02_Methylation_Risk_score.R

hc <- read.csv("Results/MRS/GW_GSE50660.csv")
hc <- hc[,c(1,20:23)]
names(hc) <- c("ID", "smoking_sm13", "age", "gender", "Smoking")

#From GEO: smoking (0, 1 and 2, which represent never, former and current smokers)
hc$Smoking[hc$Smoking == 0] <- "Never"
hc$Smoking[hc$Smoking == 1] <- "Former"
hc$Smoking[hc$Smoking == 2] <- "Current smokers"
table(hc$Smoking)

#Supplementary Figure S1b

ggplot(data = hc, aes(x=Smoking, y=scale(smoking_sm13)))+ 
  geom_boxplot(notch = T, outlier.size = 0.5)+
  ylim(-4.5,8)+
  annotate("text",
           x = 1:length(table(hc$Smoking)),
           y = -4.5,
           label = table(hc$Smoking),
           col = "red", size = 1.5)+ 
  stat_compare_means(comparisons = list(c(1, 2),c(1, 3), c(2, 3)),
                     size = 1.5,
                     step.increase = 0.085,
                     tip.length = 0.02)+
  theme_bw()+
  xlab("Measured smoking habits") + 
  ylab("Smoking-Maas MRS") +
  theme(legend.position = "bottom",
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))


# Alcohol MRS-GW 
# Validated in TCGA liver data
# TCGA data was extracted and processed following scripts/01_datamanagement.R
# the MRS was constructed following scripts/02_Methylation_Risk_score.R

# Alcohol information in liver tumors from TCGA was extracted from:
# Comprehensive and Integrative Genomic Characterization of Hepatocellular Carcinoma
# Ally, Adrian et al.
# Cell, Volume 169, Issue 7, 1327 - 1341.e23


LIHC <- read.csv("Results/MRS/GW_LIHC.csv")
alc <- readxl::read_xlsx("mmc1.xlsx", skip = 3) #Ally, Adrian et al. Cell, Volume 169, Issue 7, 1327 - 1341.e23
alc <-as.data.frame(alc[substr(alc$Barcode, 14,15) == "01", ])
alc$Barcode <- stringr::str_extract(alc$Barcode, "[^-]*-[^-]*-[^-]*")

liver <- merge(LIHC, alc, by.x="X", by.y="Barcode")
table(liver$`Alcoholic liver disease`)
liver$`Alcoholic liver disease`[liver$`Alcoholic liver disease` == "Alcohol"] <- "Yes"
sub <- liver[!liver$`Alcoholic liver disease`== "---",]

#Supplementary Figure S1a

ggplot(data = sub, aes(x=`Alcoholic liver disease`, y=scale(Alcohol))) +
  geom_boxplot(notch = T, outlier.size = 0.5)+
  ylim(-4.5,8)+
  annotate("text",
           x = 1:length(table(sub$`Alcoholic liver disease`)),
           y = -4.5,
           label = table(sub$`Alcoholic liver disease`),
           col = "red", size=1.5)+
  theme_bw()+
  xlab("Alcoholic liver disease") + 
  ylab("Alcohol MRS-GW")+
  stat_compare_means(size = 1.5,
                     step.increase = 0.085,
                     tip.length = 0.02)+
  theme(legend.position = "bottom",
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))


table(liver$history_hepato_carcinoma_risk_factor)

liver$history_hepato_carcinoma_risk_factor[liver$history_hepato_carcinoma_risk_factor == "Alcohol consumption"] <- "Alcohol"
liver$history_hepato_carcinoma_risk_factor[liver$history_hepato_carcinoma_risk_factor == "No History of Primary Risk Factors"] <- "None"

risk <- c("None", "Alcohol")
sub <- liver[liver$history_hepato_carcinoma_risk_factor %in% risk, ]

ggplot(data = sub, aes(x=history_hepato_carcinoma_risk_factor, y=scale(Alcohol))) +
  geom_boxplot(notch = T, outlier.size = 0.5)+
  ylim(-4.5,8)+
  annotate("text",
           x = 1:length(table(sub$history_hepato_carcinoma_risk_factor)),
           y = -4.5,
           label = table(sub$history_hepato_carcinoma_risk_factor),
           col = "red", size=1.5)+
  theme_bw()+
  
  xlab("Primary risk factor") + 
  ylab("Alcohol MRS-GW")+
  stat_compare_means(size = 1.5,
                     step.increase = 0.085,
                     tip.length = 0.02)+
  theme(legend.position = "bottom",
        #       axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))

