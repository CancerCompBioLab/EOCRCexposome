### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


## Pesticide use estimates are downloaded from the National Water-Quality Assessment (NAWQA)

year <- seq(1992, 2012, by=1)

for (i in year){
  yearfile <- paste("https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/PesticideUseEstimates/EPest.county.estimates.", i, ".txt", sep = "")
  myfile2 <- read.table(url(yearfile), sep= "\t", header=T)
  yeardown <- paste("Est_Ann_Agricultural_Pesticide_", i, ".csv", sep = "")
  write.csv(myfile2, yeardown, row.names = F)
}


# Pesticides included in the MRS analysis but are banned in the USA and therefore not measured: 
# HEPTACHLOR  
# TOXAPHENE 
# Chlordane
# DDT


## Subset the data per pesticide and combine for all years together

pestices <- c("ACETOCHLOR", "2,4-D", "ATRAZINE", "DICAMBA", "GLYPHOSATE", "LINDANE", "MALATHION", "MESOTRIONE", "METOLACHLOR", "PICLORAM")  

overview <- matrix(ncol= 6, nrow=1)
colnames(overview) <- c("COMPOUND", "YEAR", "STATE_FIPS_CODE", "COUNTY_FIPS_CODE", "EPEST_LOW_KG", "EPEST_HIGH_KG")



for (i in year){
  yeardown <- paste("Est_Ann_Agricultural_Pesticide_", i, ".csv", sep = "")
  pesti <- read.csv(yeardown, header=T)
  overview <- rbind(overview, pesti)
}  
  
for (j in pestices){
    newoverview <-overview
    subnewpesti <- newoverview[newoverview$COMPOUND == j, ]
    if (nrow(subnewpesti) >0) {
    file <- paste(j, "1992_2012", ".csv", sep= "_")
}
  
  write.csv(subnewpesti, file, row.names = F)
}


#########

# Download the Fip codes


download.file("https://www2.census.gov/programs-surveys/popest/geographies/2018/all-geocodes-v2018.xlsx", destfile="Fipcodes.xlsx")
codes<-  openxlsx::read.xlsx("Fipcodes.xlsx", startRow = 5, colNames = TRUE)
codes_state <- codes[codes$Summary.Level == "040", ]
codes_state <- codes_state[,c(2,7)]
names(codes_state)[2] <- "State"

county <- codes[codes$Summary.Level == "050",]
code <- merge(county, codes_state, by ="State.Code.(FIPS)")
codes <- code[,c("State.Code.(FIPS)", "County.Code.(FIPS)", "Area.Name.(including.legal/statistical.area.description)", "State")]
names(codes) <- c("STATE_FIPS_CODE", "COUNTY_FIPS_CODE", "Area_Name", "State")
#write.csv(codes, "codes_state_county.csv", row.names = F)


############

# Download county level information 


download.file("https://www2.census.gov/geo/docs/reference/ua/2020_UA_COUNTY.xlsx", destfile="Countyinfo.xlsx")
counties<-  openxlsx::read.xlsx("Countyinfo.xlsx", colNames = TRUE)
countySize <- counties[,c("STATE", "COUNTY",  "STATE_NAME", "COUNTY_NAME", "ALAND_Mi²_COU")]
names(countySize) <- c("STATE", "COUNTY", "state", "county", "squaresmiles")
#write.csv(countySize, "countySize.csv", row.names = F)


#############


# Date management and combining pesticide data with eoCRC incidence rates 

library(stringr)
library(performance)


pest <- c("ACETOCHLOR", "2,4-D", "ATRAZINE", "DICAMBA", "GLYPHOSATE", "LINDANE",  "MALATHION", "MESOTRIONE", "METOLACHLOR", "PICLORAM")  

for (i in pest){

  file <- paste0(i , "_1992_2012_.csv")
  dicam <- read.csv(file, header=T)
  dicam <- dicam[-1, ]
  
  codes <- read.csv("codes_state_county.csv", header = T)
  dicams <- merge(dicam, codes, by=c("STATE_FIPS_CODE", "COUNTY_FIPS_CODE"))
  dicams$Area_Name <- sub(' County*', '', dicams$Area_Name)
  dicams$SC <- paste(dicams$State, dicams$Area_Name, sep = "-")
  
  # select columns State-County (SC), measures years, EPEST_HIGH_KG 
  dicams <- dicams[,c(9,4,6)]
  
  size <- read.csv("countySize.csv", header=T)
  size$SC <- paste(size$state,size$county, sep = "-")
  size <- size[,c(6,5)]
  names(size) <- c("SC", "Area")
  dicams <- merge(dicams, size, by = "SC")
  
  # Calculate the Pesticide use intensity by dividing the pesticide use with the county area in squaresmiles
  dicams$perArea <- dicams$EPEST_HIGH_KG/dicams$Area
  
  
  # SEER Research Plus data is extracted using the SEERstat software. 
  # For access to this data, please visit https://seer.cancer.gov/data/access.html
  # 
  # We extracted the age-Adjusted CRC incidence rates per year 
  # from SEER 8 and SEER 12, per county, 
  # patients aged 25-29, 30-34, 35-39, 40-44, 45-49 years
  # Site And Morphology: Site recode ICD-0-3-WHO 2008 - Colon and Rectum
  
  
  # *** file does not exist ***
  seer <- read.csv("SEER8.csv", header=T) # same for SEER12
  
  SeNOcompl <- merge(seer, dicams, by.x= c("SC", "Year"), by.y= c("SC", "YEAR"))
  new <- SeNOcompl[complete.cases(SeNOcompl$EPEST_HIGH_KG),]
  
  SC <- unique(new$SC)
  year <- unique(new$Year)
  age <- matrix(ncol= length(year), nrow= length(SC))
  row.names(age) <- SC
  colnames(age) <- year
  
  for (j in SC){
    data <- new
    data <- data[data$SC == j,]
    
    for (k in year){
      daty <- data
      daty <- daty[daty$Year == k, ]
      
      #combine the age categories to obtain the total eoCRC between 25 and 49, adjusted for the population
      
      Age_ad <- round((sum(daty$Count) / sum(daty$Population)) * 100000, digits = 3)
      
      m <- as.character(k)
      l <- as.character(j)
      age[l,m] <- Age_ad
      
    }
    
  }
  
  
  # only keep counties that have less than 50% of the years without a case

  p5 <- length(unique(year))*0.5
  age2 <- age[rowSums(age == 0) <= p5, ]
  age2 <- as.data.frame(age2[complete.cases(age2),])
  age2$SC <- rownames(age2)
  
  data <- reshape::melt(age2, id="SC")
  new <- merge(data, dicams, by.x= c("SC", "variable"), by.y= c("SC", "YEAR"))
  
  file3 <- paste0("SEER_pesticide/", i , "_seer.csv")
  
  write.csv(new, file3, row.names = F)
  

}



library(stringr)
library(performance)

# LINDANE not enough complete cases for analysis
# MALATHION and METOLACHLOR do not converge

pest <- c("ACETOCHLOR", "2,4-D", "ATRAZINE", "DICAMBA", "GLYPHOSATE", "MESOTRIONE", "PICLORAM") 
pesti <- matrix(ncol= 2, nrow= length(pest))
row.names(pesti) <- pest
colnames(pesti) <- c("Beta", "pval")

library(dplyr)
library(ggplot2)


for (i in pest){
  file <- paste0(i , "_seer.csv")
  new <- read.csv(file, header=T)
  new$logArea <- log(new$perArea)

  # remove the top and bottom 5% 
  pise3S <- new %>% filter(between(logArea, 
                                   quantile(logArea, 0.05), quantile(logArea, 0.95)))
  
  m3N1 <- lmerTest::lmer(value ~ logArea  + variable + (1|SC), pise3S)
  pesti[i, "Beta"] <-  summary(m3N1)$coefficients["logArea", "Estimate"]
  pesti[i, "pval"] <-  summary(m3N1)$coefficients["logArea", "Pr(>|t|)"]
  
}

pesti$pesticide <- rownames(pesti)

write.csv(pesti, "SEER_pesticide/pesticide_LMM_seer8_5percTB.csv", row.names = F) # SEER12


#
#
# Figure 4
#
#




pesti8 <- read.csv("Results/SEER_pesticide/pesti_comb_LMM_seer8_5percTB.csv",header=T)
pesti8$SEER <- "SEER 8"


pesti12 <- read.csv("Results/SEER_pesticide/pesti_comb_LMM_seer12_5percTB.csv",header=T)
pesti12$SEER <- "SEER 12"


comb2 <- rbind(pesti8, pesti12)
comb2$pval <- -log10(comb2$pval)


comb2$pesticide[comb2$pesticide == "X2.4.D"] <- "2,4-D"
comb2$pesticide[comb2$pesticide == "ACETOCHLOR"] <- "Acetochlor"
comb2$pesticide[comb2$pesticide == "ATRAZINE"] <- "Atrazine"
comb2$pesticide[comb2$pesticide == "DICAMBA"] <- "Dicamba"
comb2$pesticide[comb2$pesticide == "GLYPHOSATE"] <- "Glyphosate"
comb2$pesticide[comb2$pesticide == "MESOTRIONE"] <- "Mesotrione"
comb2$pesticide[comb2$pesticide == "PICLORAM"] <- "Picloram"



### Figure 4a

ggplot(comb2, aes(x = factor(SEER, levels=c("SEER 8", "SEER 12")), y = pval, fill = SEER)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(name = "SEER", values = c("SEER 8"   =  "#98A7D4",     "SEER 12"    =  "#FDAE6B"))+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=5),
        axis.title.y = element_text(size = 5),axis.title.x=element_blank(),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50"),  panel.spacing.x = unit(0.1, "cm"), 
        strip.text.x = element_text(face = "italic", size = 5, margin = unit(rep(10, 4), "pt")))+
  
  labs(y="-log10(P-value)")+ 
  geom_hline(yintercept=1.3, linetype="dashed", color = "maroon", linewidth=0.5)+
  ylim(0,5.1)+
  facet_wrap(nrow=1, vars(factor(pesticide, levels= c("Glyphosate", "Atrazine", "Picloram", "2,4-D", "Dicamba", "Mesotrione", "Acetochlor"))))



########## 

### Adjusting Glyphosate, Atrazine, and Picloram

# MAKE THE DF WITH RESULTS
pest <- c("ATR", "ATR_G","ATR_P", "ATR_GP", 
          "GLYP","GLYP_A", "GLYP_P", "GLYP_AP",
          "PIC", "PIC_A", "PIC_G", "PIC_AG") 

compesti <- matrix(nrow= length(pest), ncol= 2)
row.names(compesti) <- pest
colnames(compesti) <- c("Beta", "Pval")


################
#
# PICLORAM
#
###############

PICLORAM <- read.csv("PICLORAM_seer.csv", header=T)
PICLORAM$logArea <- log(PICLORAM$perArea)
PICLORAM <- PICLORAM %>% filter(between(logArea, 
                                        quantile(logArea, 0.05), quantile(logArea, 0.95)))
piclo <- PICLORAM[,c(1,2,3,7)]
head(piclo)
names(piclo)[4] <- "logPiclo"
m3N1 <- lmerTest::lmer(value ~ logPiclo  + variable + (1|SC), piclo)
summary(m3N1)


GLYPHOSATE <- read.csv("GLYPHOSATE_seer.csv", header=T)
GLYPHOSATE$logArea <- log(GLYPHOSATE$perArea)
head(GLYPHOSATE)
glyp <- GLYPHOSATE[,c(1,2,7)]
names(glyp)[3] <- "logGlyp"

ATRAZINE <- read.csv("ATRAZINE_seer.csv", header=T)
ATRAZINE$logArea <- log(ATRAZINE$perArea)
atraz <- ATRAZINE[,c(1,2,7)]
names(atraz)[3] <- "logAtraz"
head(atraz)

pesti <- merge(piclo, atraz, by=c("SC", "variable"), all = T)
pesti3 <- merge(pesti, glyp, by=c("SC", "variable"), all = T)
pesti3 <- pesti3[complete.cases(pesti3),]

PIC <- lmerTest::lmer(value ~ logPiclo  + variable + (1|SC), pesti3)
results_df <-summary(PIC)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["PIC",1] <- round(results_df[1,1], digits = 2)
compesti["PIC",2] <- format(results_df[1,5], scientific = T, digits = 3)


PICA <- lmerTest::lmer(value ~ logPiclo  + logAtraz  + variable + (1|SC), pesti3)
results_df <-summary(PICA)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["PIC_A",1] <- round(results_df[1,1], digits = 2)
compesti["PIC_A",2] <- format(results_df[1,5], scientific = T, digits = 3)


PICG <- lmerTest::lmer(value ~ logPiclo  + logGlyp  + variable + (1|SC), pesti3)
results_df <-summary(PICG)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["PIC_G",1] <- round(results_df[1,1], digits = 2)
compesti["PIC_G",2] <- format(results_df[1,5], scientific = T, digits = 3)


PICAG <- lmerTest::lmer(value ~ logPiclo  + logGlyp + logAtraz + variable + (1|SC), pesti3)
results_df <-summary(PICAG)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["PIC_AG",1] <- round(results_df[1,1], digits = 2)
compesti["PIC_AG",2] <- format(results_df[1,5], scientific = T, digits = 3)

###################
#
# GLYPHOSATE
#
###################

GLYPHOSATE <- read.csv("GLYPHOSATE_seer.csv", header=T)
GLYPHOSATE$logArea <- log(GLYPHOSATE$perArea)
GLYPHOSATE <- GLYPHOSATE %>% filter(between(logArea, 
                                            quantile(logArea, 0.05), quantile(logArea, 0.95)))
glyp <- GLYPHOSATE[,c(1,2,3,7)]
names(glyp)[4] <- "logGlyp"

PICLORAM <- read.csv("PICLORAM_seer.csv", header=T)
PICLORAM$logArea <- log(PICLORAM$perArea)
piclo <- PICLORAM[,c(1,2,7)]
names(piclo)[3] <- "logPiclo"

ATRAZINE <- read.csv("ATRAZINE_seer.csv", header=T)
ATRAZINE$logArea <- log(ATRAZINE$perArea)
atraz <- ATRAZINE[,c(1,2,7)]
names(atraz)[3] <- "logAtraz"


pesti <- merge(piclo, atraz, by=c("SC", "variable"), all = T)
pesti3 <- merge(pesti, glyp, by=c("SC", "variable"), all = T)
pesti3 <- pesti3[complete.cases(pesti3), ]


GLYP <- lmerTest::lmer(value ~ logGlyp  + variable + (1|SC), pesti3)
results_df <-summary(GLYP)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["GLYP",1] <- round(results_df[1,1], digits = 2)
compesti["GLYP",2] <- format(results_df[1,5], scientific = T, digits = 3)


GLYPA <- lmerTest::lmer(value ~ logGlyp + logAtraz  + variable + (1|SC), pesti3)
results_df <-summary(GLYPA)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["GLYP_A",1] <- round(results_df[1,1], digits = 2)
compesti["GLYP_A",2] <- format(results_df[1,5], scientific = T, digits = 3)

GLYPP <- lmerTest::lmer(value ~ logGlyp  +  logPiclo + variable + (1|SC), pesti3)
results_df <-summary(GLYPP)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["GLYP_P",1] <- round(results_df[1,1], digits = 2)
compesti["GLYP_P",2] <- format(results_df[1,5], scientific = T, digits = 3)

GLYPAP <- lmerTest::lmer(value ~  logGlyp + logPiclo  + logAtraz + variable + (1|SC), pesti3)
results_df <-summary(GLYPAP)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["GLYP_AP",1] <- round(results_df[1,1], digits = 2)
compesti["GLYP_AP",2] <- format(results_df[1,5], scientific = T, digits = 3)


###################
#
# Atrazine
#
###################

ATRAZINE <- read.csv("ATRAZINE_seer.csv", header=T)
ATRAZINE$logArea <- log(ATRAZINE$perArea)
ATRAZINE <- ATRAZINE %>% filter(between(logArea, 
                                        quantile(logArea, 0.05), quantile(logArea, 0.95)))
atraz <- ATRAZINE[,c(1,2,3,7)]
names(atraz)[4] <- "logAtraz"

GLYPHOSATE <- read.csv("GLYPHOSATE_seer.csv", header=T)
GLYPHOSATE$logArea <- log(GLYPHOSATE$perArea)
glyp <- GLYPHOSATE[,c(1,2,7)]
names(glyp)[3] <- "logGlyp"

PICLORAM <- read.csv("PICLORAM_seer.csv", header=T)
PICLORAM$logArea <- log(PICLORAM$perArea)
piclo <- PICLORAM[,c(1,2,7)]
names(piclo)[3] <- "logPiclo"

pesti <- merge(piclo, atraz, by=c("SC", "variable"), all = T)
pesti3 <- merge(pesti, glyp, by=c("SC", "variable"), all = T)
pesti3 <- pesti3[complete.cases(pesti3), ]

## glyphosate and atrazine are strongly correlated
ggplot(pesti3, aes(x=logAtraz, y=logGlyp))+
  geom_point()+
  stat_cor()+
  theme_bw()

ATR <- lmerTest::lmer(value ~ logAtraz  + variable + (1|SC), pesti3)
results_df <-summary(ATR)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["ATR",1] <- round(results_df[1,1], digits = 2)
compesti["ATR",2] <- format(results_df[1,5], scientific = T, digits = 3)

ATRP <- lmerTest::lmer(value ~ logAtraz + logPiclo  + variable + (1|SC), pesti3)
results_df <-summary(ATRP)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["ATR_P",1] <- round(results_df[1,1], digits = 2)
compesti["ATR_P",2] <- format(results_df[1,5], scientific = T, digits = 3)


ATRG <- lmerTest::lmer(value ~  logAtraz + logGlyp  + variable + (1|SC), pesti3)
results_df <-summary(ATRG)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["ATR_G",1] <- round(results_df[1,1], digits = 2)
compesti["ATR_G",2] <- format(results_df[1,5], scientific = T, digits = 3)

ATRPG <- lmerTest::lmer(value ~  logAtraz + logGlyp  + logPiclo  + variable + (1|SC), pesti3)
results_df <-summary(ATRPG)$coefficients
results_df <- as.data.frame(t(results_df[2,]))
compesti["ATR_GP",1] <- round(results_df[1,1], digits = 2)
compesti["ATR_GP",2] <- format(results_df[1,5], scientific = T, digits = 3)

#write.csv(compesti, "Results/SEER_pesticide/SEER8_adjustment.csv")

pesti <- read.csv("Results/SEER_pesticide/SEER8_adjustment.csv", header=T)
names(pesti)[1] <- "pesticide"
pesti$pesti <- gsub("\\_.*","",pesti$pesticide)
pesti$mod <- gsub(".*_","",pesti$pesticide)
pesti$Pval <- -log10(pesti$Pval)

# Figure 6: Pesticide use intensity is associated with early-onset colorectal cancer incidence.
# Figure 6b

ggplot(pesti, aes(x =factor(pesticide, levels=c("ATR",  "ATR_P", "ATR_G", "ATR_GP", 
                                                "GLYP","GLYP_P", "GLYP_A", "GLYP_AP", 
                                                "PIC", "PIC_A",  "PIC_G", "PIC_AG")) ,
                  y = Pval, fill = mod)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(name = "mod", values = c("G"   =  "#FFC08AFF",     "A"    =  "#D94801", "P"= "#802A07FF", 
                                             "AG" = "#5F9ED1FF",   "AP"= "#5F9ED1FF",      "GP"= "#5F9ED1FF",    
                                             "ATR" = "#98A7D4",   "GLYP" = "#98A7D4", "PIC"= "#98A7D4"))+
  scale_x_discrete(labels = c("ATR"    = "Baseline",     "ATR_P" = "Picloram", "ATR_G" = "Glyphosate", "ATR_GP" = "Combined",
                              "GLYP"  = "Baseline",     "GLYP_P" = "Picloram", "GLYP_A" = "Atrazine", "GLYP_AP" = "Combined",
                              "PIC"   = "Baseline",     "PIC_A" = "Atrazine", "PIC_G" = "Glyphosate", "PIC_AG" = "Combined")) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=5),
        axis.title.y = element_text(size = 5),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50"),  panel.spacing.x = unit(0.1, "cm"), 
        strip.text.x = element_text(face = "italic", size = 5, 
                                    margin = unit(rep(10, 4), "pt")))+
  labs(y="-log10(P-value)")+ 
  geom_hline(yintercept=1.3, linetype="dashed", color = "maroon", linewidth=0.5)+
  ylim(0,5.1)+
  facet_wrap(nrow=1, vars(factor(pesti, levels= c("GLYP", "ATR", "PIC"), labels = c("Glyphosate", "Atrazine",  "Picloram"))),  scales="free_x")


