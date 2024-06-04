### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
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



# Figure S1 b

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
