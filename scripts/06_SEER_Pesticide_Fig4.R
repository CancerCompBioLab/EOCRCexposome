### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
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



### Figure 4C

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





#### Figure 4A


# Select the directory where you stored the processed SEER and pesticide data. 

file2 <-  list.files("SEER_pesticide/",  pattern =  "_seer.csv")  # Repeat same for SEE12
tables <- lapply(file2, read.csv, header = TRUE)

ID <- gsub("\\_.*","",file2)
tables <- mapply(cbind, tables, "Pesticide"=ID, SIMPLIFY=F)
comb <- do.call(rbind , tables)
comb$State <- gsub("\\-.*","",comb$SC)
comb2 <- comb[,c(1,2,7)]


#
# As the SEER data is not publicly available we have uploaded the files including the 
# needed pesticide information to replicate figures 4a
#
comb2 <- read.csv("Results/SEER_pesticide/inputFig4aSEER8.csv", header=T)


library(dplyr)
df_summary <- comb2 %>%
  group_by(Pesticide, variable) %>%
  summarize(Counties = n_distinct(SC)) %>%
  ungroup()

seer8 <- ggplot(df_summary, aes(x = variable, y = Pesticide, fill = Counties)) +
  geom_tile() + 
  scale_fill_gradient2(low ="#f3f5f9", mid = "#ccd2e5", high = "#5e6783", midpoint = 95, 
                       limits = c(55,115), 
                       breaks = c(55, 75, 95, 115),
                       labels = c(55, 75, 95, 115),
                       name = "No. of Counties") +
  facet_wrap(~ Pesticide, ncol = 1, scales = "free_y") +# Separate plot for each pesticide
  theme_minimal() +
  scale_x_continuous(name="Years measured", breaks= seq(1992, 2012, 5)) +
  
  labs(title = "SEER 8",
       x = "Year", 
       y = "Pesticides", 
       fill = "No. of Counties")+ 
  theme(legend.position = c(0.01, -0.4),
         legend.direction="horizontal",
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -2, size=5),
        panel.spacing.y = unit(0.01, "lines"),
        axis.text=element_text(size=5),
        axis.text.x = element_text(angle = 45),
        legend.text=element_text(size=5,angle = 45),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.box.margin=margin(-5,0,-5,0),
        legend.key.size = unit(0.25, 'cm'),
        plot.margin= unit(c(0,0.2,0.7,0.1), "cm")
        
  )


## Repeat for the SEER 12 pesticide files

#
# As the SEER data is not publicly available we have uploaded the files including the 
# needed pesticide information to replicate figures 4a
#
comb12 <- read.csv("Results/SEER_pesticide/inputFig4aSEER12.csv", header=T)


library(dplyr)
df12_summary <- comb12 %>%
  group_by(Pesticide, variable) %>%
  summarize(Counties = n_distinct(SC)) %>%
  ungroup()


seer12 <- 
  ggplot(df12_summary, aes(x = variable, y = Pesticide, fill = Counties)) +
  geom_tile() + 
  scale_fill_gradient2(low ="#FFF5EB", mid = "#FDAE6B", high = "#D08F58", midpoint = 105, 
                       limits = c(55,115), 
                       breaks = c(55, 75, 95, 115),
                       labels = c(55, 75, 95, 115),
                       name = "No. of Counties") +
  facet_wrap(~ Pesticide, ncol = 1, scales = "free_y") +# Separate plot for each pesticide
  theme_minimal() +
  scale_x_continuous(name="Years measured", breaks= seq(1992, 2012, 5)) +
  
  labs(title = "SEER 12",
       x = "Year", 
       y = "Pesticides", 
       fill = "No. of Counties")+ 
  theme(legend.position = c(0.01, -0.4),
          legend.direction="horizontal",
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -2, size=5),
                panel.spacing.y = unit(0.01, "lines"),
        axis.text=element_text(size=5),
        axis.text.x = element_text(angle = 45),
        legend.text=element_text(size=5,angle = 45),
        legend.title = element_text(size=5),
        axis.title = element_text(size = 5),
        legend.box.margin=margin(-5,0,-5,0),
        legend.key.size = unit(0.25, 'cm'),
        plot.margin= unit(c(0,0.2,0.7,0.1), "cm")
        
  )



## Combine the plots from SEER8 and SEER12 -> Figure 4A

plots8_12 <-cowplot::plot_grid(seer8, seer12, nrow =   1)  



#
#
#
#   Figure 4b
#
#
#


#install.packages("usmap")
library(usmap) 
library(dplyr)
library(ggplot2)

# As the SEER data used in this plot is not publicly available, we have included random number between 0 and 100 to the age adjusted rates 

i <- "PICLORAM"
file <- paste0("SEER_pesticide/", i , "_seer.csv")  # File is NOT available -> input your own generated file
new <- read.csv(file, header=T)
new$logArea <- log(new$perArea)
pise3S <- new %>% filter(between(logArea, 
                                 quantile(logArea, 0.05), quantile(logArea, 0.95)))

#pise3S$value <- as.numeric(sample(100, size = nrow(pise3S), replace = TRUE)) -> we included random number between 0 and 100 to the age adjusted rates to reproduce the figure

pise3S <- read.csv("Results/SEER_pesticide/PICLORAM_seer_Fake.csv", header=T)

pico <- as.data.frame(aggregate(pise3S$logArea, list(pise3S$SC), FUN=mean))
value <-  aggregate(pise3S$value, list(pise3S$SC), FUN=mean) 
mean <- left_join(pico, value, by = "Group.1")
names(mean) <- c("SC", "logArea", "value")
mean$State <- gsub("\\-.*","",mean$SC)
mean$county <- gsub(".*-","",mean$SC)

# Names need to match with the information in the usmap package
mean$county <- tolower(mean$county) 
mean$State <- tolower(mean$State)
SC_pico <- unique(mean$SC)
states_pico <- unique(mean$State)
mean$subregion <- mean$county
mean$subregion[mean$subregion =="dekalb"] <- "de kalb"

state <- map_data("state")
counties <- map_data("county")
names(mean)[4] <- "region"
unicoun <-   counties[,c(5,6)]   
unicoun$SC <- paste(unicoun$region, unicoun$subregion, sep = "-")
unicoun <- unicoun[!duplicated(unicoun$SC), ]
mean$SC2 <- paste(mean$region, mean$subregion, sep = "-")
Unicoun <- unicoun[! unicoun$SC %in% mean$SC2, ]
Unicoun$logArea <- 0
Unicoun$value <- 0
count <- Unicoun[,c(3,4,5,1,2)]
mean2 <- mean[,c(7,2,3,4,6)]
names(mean2)[1] <- "SC" 

Counties_p <- rbind(count, mean2)

States_sub <- subset(state, region %in% states_pico[2])
counties_sub <- subset(counties, region %in% states_pico[2])

pico.map2 <- merge(counties_sub, Counties_p, by=c("region","subregion"))
pico.map[pico.map == 0] <- NA

map1 <-  ggplot(data=States_sub, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.5) + 
  geom_polygon(color="black", fill="white") + 
  geom_polygon(data=pico.map, 
               aes(fill = value), 
               color="white") + 
  geom_polygon(color="black", fill=NA) +
   # note that these number are not the same as used in the Figure in the manuscript
  
    scale_fill_gradient2(name = "IR eoCRC",  low = "#F7FBFF", mid = "#9ECAE1",  high = "#061746",na.value = "grey95", limits = c(35,65), midpoint = 50, 
                       breaks = c(35, 50, 65),
                       labels = c(35, 50, 65))+
  theme_void()+
  theme(legend.position = "left", legend.text=element_text(size=5),
        legend.title = element_text(size=5),legend.key.size = unit(0.15, 'cm'),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.box.spacing = unit(-2.5, "pt"))
map1 

map2 <-   ggplot(data=States_sub, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.5) + 
  geom_polygon(color="black", fill="white") + 
  geom_polygon(data=pico.map, aes(fill = logArea), color="white") + 
  geom_polygon(color="black", fill=NA) +
  scale_fill_gradient2(name= "Picloram",low = "#FFF5EB", mid = "#D94801", high = "#662105",na.value = "grey95", limits = c(-5,0), midpoint = -2.7,
                       breaks=c(-5,0),labels=c("low","high"))+
  theme_void()+ 
  theme(legend.position = "left", legend.text=element_text(size=5),
        legend.title = element_text(size=5),legend.key.size = unit(0.15, 'cm'),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.box.spacing = unit(-2.5, "pt"))#

map2

plot <-cowplot::plot_grid(map1, map2, ncol =  1)  

# Figure 4b
plot


