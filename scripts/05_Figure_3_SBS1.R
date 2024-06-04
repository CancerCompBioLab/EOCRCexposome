### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)






# Input: Matrix with CpGs in rows
met_mat_M  <- readRDS("Startdata/COAD_Age_adjusted_DNAm_Bval_overlap_MRSCpGs.rds")
met_mat_M[1:4,1:4]

clinic <- read.csv("Startdata/Clinic_Coad.csv", header=T)
head(clinic)

clinic$EOCRC = ifelse(clinic$age <50, "eoCRC",
                      ifelse(clinic$age>=50 & clinic$age <70, "moCRC",
                             ifelse(clinic$age>=70, "loCRC", NA)))

clinic <- clinic[clinic$histologic_diagnosis == "Colon Adenocarcinoma", ]
clinic <- clinic[complete.cases(clinic$EOCRC), ]
met_mat_M <- met_mat_M[, colnames(met_mat_M) %in% clinic$ID]
clinic <- clinic[clinic$ID %in% colnames(met_mat_M), ]

clinic <- clinic[clinic$MsiStatus == "MSS", ]
clinic <- clinic[complete.cases(clinic$MsiStatus), ]
clinic <- clinic[complete.cases(clinic$SBS1), ]

summary(clinic$SBS1)
mean(clinic$SBS1)+ (3*sd(clinic$SBS1))

## Exclude patients with mean +/- 3*SD SBS1
clinic <- clinic[clinic$ID != "TCGA-A6-6654",]

table(clinic$EOCRC)

# Fig3e: SBS1-score distribution over the 3 age categories
library(ggplot2)

ggplot(clinic, aes(x=factor(EOCRC, level = c("eoCRC", "moCRC", "loCRC")),  y=SBS1)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = gender), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Sex:",
                    values = c("#FDAE6B", "#98A7D4"), 
                    labels= c("Female", "Male"))+
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
  scale_y_continuous(name ="SBS1 Score")+
  scale_x_discrete(name ="Age category at diagnosis", labels=c("Early-onset","Middle-onset","Later-onset"))+
  labs(x = "New x label")

# Check SBS1-score distribution across the 3 age categories 
tapply(clinic$SBS1, factor(clinic$EOCRC, level = c("eoCRC", "moCRC", "loCRC")), summary)

# The median in moCRC is 59.5 and the 1Qr in loCRC is 60 -> SBS cut off at 60

clinic$mutCRC = ifelse(clinic$SBS1 <60, "eoCRC",
                       ifelse(clinic$SBS1>=60, "loCRC", NA))



# Calculate the methylation risk scores for this subset of patients 

Pheno = read.csv(file="InputData/CombGW.csv", header=TRUE)
pheno <- Pheno[,c("CpG", "Beta", "Trait")]
pheno <- pheno[complete.cases(pheno$Trait), ]
pheno <- pheno[as.numeric(ave(pheno$Trait, pheno$Trait, FUN=length)) >= 2, ]
phenos <- unique(pheno$Trait)

Scores <- matrix(ncol= length(phenos), nrow= ncol(met_mat_M))
colnames(Scores) <- phenos
row.names(Scores) <- colnames(met_mat_M)


for (j in phenos){
    pheno1 <- pheno[pheno$Trait == j, ]
  
  # use only CpGs present in both
  DNAm <- as.data.frame(met_mat_M[rownames(met_mat_M) %in% pheno1$CpG, ])
  pheno1 <- as.data.frame(pheno1[pheno1$CpG %in% rownames(met_mat_M), ])
  df2 <- DNAm[order(rownames(DNAm),decreasing=TRUE),]
  df3 <- pheno1[order(pheno1$CpG,decreasing=TRUE),]
  scores <- as.numeric(rep(NA,ncol(df2)))
  for (i in 1:ncol(df2)){
    scores[i]<-sum(df2[,i]*df3$Beta)}  
  Scores[ ,j] <- scores
}

Scores <- as.data.frame(Scores)

library(dplyr)
Scores <- Scores %>% mutate_all(~(scale(.) %>% as.vector))
Scorescl  <- merge(Scores, clinic, by.x= 0, by.y="ID")


Scorescl$gender[Scorescl$gender == "FEMALE"] <- 1
Scorescl$gender[Scorescl$gender == "MALE"] <- 0
Scorescl$gender <- as.factor(Scorescl$gender)


## Fig3f MRS picloram vs SBS 1 score threshold of 60 


ggplot(Scorescl, aes(x=factor(mutCRC, level=c("loCRC", "eoCRC")), y=PICLORAM)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = gender), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Sex:",
                    values = c("#FDAE6B", "#98A7D4"), 
                    labels= c("Female", "Male"))+
  ylim(-2.5,3)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))+ #change legend key size)+ 
  #ggtitle("MRS Education COAD") +
  #xlab("Early or later-onset colon cancer") + 
  ylab("Methylation Risk Score") +
  scale_x_discrete(labels=c('SBS1 >60', 'SBS1 <60'))


Scorescl$age <- as.numeric(Scorescl$age)

## Fig3g Chronological age vs SBS1 

ggplot(Scorescl, aes(x=factor(mutCRC, level=c("loCRC", "eoCRC")), y=age)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = gender), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Sex:",
                    values = c("#FDAE6B", "#98A7D4"), 
                    labels= c("Female", "Male"))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))+ #change legend key size)+ 
  ylab("Age at diagnosis") +
  scale_x_discrete(labels=c('SBS1 >60', 'SBS1 <60'))



str(Scorescl)
Scorescl$gender[Scorescl$gender == "FEMALE"] <- 1
Scorescl$gender[Scorescl$gender == "MALE"] <- 0
Scorescl$gender <- as.factor(Scorescl$gender)

Scorescl$mutCRC[Scorescl$mutCRC == "eoCRC"] <- 1
Scorescl$mutCRC[Scorescl$mutCRC == "loCRC"] <- 0
Scorescl$mutCRC <- as.factor(Scorescl$mutCRC)


table(Scorescl$mutCRC)

### SBS1 categorization threshold 60 vs picloram MRS
test <- glm(mutCRC ~ PICLORAM + gender, family = binomial, data = Scorescl) 
print(summary.glm(test))

exp(0.6123)
exp(confint(test)["PICLORAM","2.5 %"])
exp(confint(test)["PICLORAM","97.5 %"])



## Picloram in Early vs later onset in this subset of patients that have SBS information, methylation data, and MSS

ggplot(Scorescl[Scorescl$EOCRC!="moCRC",], aes(x=factor(EOCRC, level=c("loCRC", "eoCRC")), y=PICLORAM)) + 
  geom_boxplot(outlier.shape = NA, notch = T)+ 
  geom_dotplot(dotsize =0.65,aes(fill = gender), binaxis='y', stackdir='center',
               position=position_dodge(0.35))+
  scale_fill_manual(name = "Sex:",
                    values = c("#FDAE6B", "#98A7D4"), 
                    labels= c("Female", "Male"))+
  ylim(-2.5,3)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.title = element_text(size=5),
        axis.title.y = element_text(size = 5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-5,0),
        legend.key.size = unit(0.25, 'cm'))+ #change legend key size)+ 
  ylab("Methylation Risk Score") +
  scale_x_discrete(labels=c('Later-onset \n N=72', 'Early-onset \n N=25'))

### Picloram MRS early vs later-onset based on age in this subset

df3sub <- Scorescl

df3sub$EOCRC[df3sub$EOCRC == "eoCRC"] <- 1
df3sub$EOCRC[df3sub$EOCRC == "loCRC"] <- 0
df3sub$EOCRC[df3sub$EOCRC == "moCRC"] <- NA
table(df3sub$EOCRC )
df3sub$EOCRC <- as.factor(df3sub$EOCRC)

test <- glm(EOCRC ~ PICLORAM + gender, family = binomial, data = df3sub) 
print(summary.glm(test))


exp(1.0969)
exp(confint(test)["PICLORAM","2.5 %"])
exp(confint(test)["PICLORAM","97.5 %"])

