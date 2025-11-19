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


overview <- matrix(ncol= 6, nrow=1)
colnames(overview) <- c("COMPOUND", "YEAR", "STATE_FIPS_CODE", "COUNTY_FIPS_CODE", "EPEST_LOW_KG", "EPEST_HIGH_KG")


year <- seq(1992, 2012, by=1)
for (i in year){
  #  i <- "2006"
  yeardown <- paste("Est_Ann_Agricultural_Pesticide_", i, ".csv", sep = "")
  pesti <- read.csv(yeardown, header=T)
  
  overview <- rbind(overview, pesti)
  
}  


colSums(is.na(overview))
df <- overview[complete.cases(overview$COMPOUND), ]

head(df)
library(dplyr)
library(tidyr)

df_wide <- df %>%
  pivot_wider(
    id_cols = c(YEAR, STATE_FIPS_CODE, COUNTY_FIPS_CODE),
    names_from = COMPOUND,
    values_from = EPEST_HIGH_KG
  )

head(df_wide)
colSums(is.na(df_wide))

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

codes <- read.csv("codes_state_county.csv", header = T)
head(codes)
names(codes) <- c("STATE_FIPS_CODE", "COUNTY_FIPS_CODE", "Area_Name", "State")


df_wide <- merge(codes, df_wide, by=c("STATE_FIPS_CODE", "COUNTY_FIPS_CODE"))
head(df_wide)

df_wide$Area_Name <- sub(' County*', '', df_wide$Area_Name)

df_wide$SC <- paste(df_wide$State, df_wide$Area_Name, sep = "-")
length(unique(df_wide$SC))

size <- read.csv("countySize.csv", header=T)
head(size)
size$SC <- paste(size$state,size$county, sep = "-")

size <- size[,c(6,5)]
names(size) <- c("SC", "Area")

Pesti <- merge(size, df_wide, by = "SC")

#File not saved in github
saveRDS(Pesti, "pesti_TOT_1992_2019.rds") 

#
#
# 
# SEER Research Plus data is extracted using the SEERstat software. 
# For access to this data, please visit https://seer.cancer.gov/data/access.html
# 
# We extracted the age-Adjusted CRC incidence rates per year 
# from SEER 8 and SEER 12, per county, 
# patients aged 25-29, 30-34, 35-39, 40-44, 45-49 years
# Site And Morphology: Site recode ICD-0-3-WHO 2008 - Colon and Rectum
#
#
#


# *** The SEER data is not available on github ***

seer <- read.csv("SEER_CRC_75_19_county_Oct2023_EO.csv", header=T)

#Keep only overlapping years between the pesticide data and CRC incidence
new <- seer[seer$Year > 1991 & seer$Year < 2013, ]

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
    l <-as.character(j)
    age[l,m] <- Age_ad
    
  }
  
}



# only keep counties that have less than 50%  (10 year) without a case

p5 <- length(unique(year))*0.5
age2 <- age[rowSums(age == 0) <= p5, ]
age2 <- as.data.frame(age2[complete.cases(age2),])
age2$SC <- rownames(age2)
data2 <- reshape::melt(age2, id="SC")

Pesti <- readRDS("/CCBdata/projects/eoCRC_sil/popdata/pesti_TOT_1992_2019.rds") 
tots <- merge(data2, Pesti, by.x= c("SC", "variable"), by.y= c("SC", "YEAR"), all.x=T)

saveRDS(tots, "SEER_pesti1992_2019_50perc_cases.rds") 

####



library(lmerTest)
library(dplyr)
library(tidyr)
library(progress)

# Data not available in github 
DF <- readRDS("SEER_pesti1992_2019_50perc_cases.rds")

ses_vars <- 3:9  # SES variables are columns 3-9

DF <- DF[,colSums(is.na(DF))<nrow(DF)]

pesticide_vars <- names(DF)[14:442]  # This should be your pesticide columns

### Do SES variables impact the association between pesticides and eoCRC?

source("Functions/pesti_ses_adjustments.R")

results <- pesticide_ses_analysis(data = DF, pesticides = pesticide_vars, ses_vars = ses_vars, 
                                  outcome = "value", time_var = "variable", group_var = "SC",
                                  filter_group = TRUE,    group_min = 10,
                                  filter_top = TRUE,      top_percent = 0.05,
                                  filter_bottom = FALSE,  bottom_percent = 0.05)

adjustment_results <- results$adjustment_results  # the p-value of the pesticide with eoCRC including SES adjustment
SES_results <- results$SES_results                # the p-value of SES variable in the model testing the pesticide with eoCRC
summary <- results$summary_table                  # which pesticides are significant in the NULL model and in how many stay significant after adjustments  


# write.csv(adjustment_results, "Results/SEER_pesticide/AllPestiAdjSES.csv", row.names = F)
# write.csv(SES_results, "Results/SEER_pesticide/SESinAdjAllPesti.csv", row.names = F)
# write.csv(summary, "Results/SEER_pesticide/Summary_allPesti_SES.csv", row.names = F)



################

# Adjust the association between pesticides and eoCRC for the other pesticides. 
# Assess whether significant pesticides lose their association when adjusted for each other, 
# and test the correlation between the two tested pesticides

# Pesticides significant in the null model
sigP <- summary[summary$NULL_Model_Significant == T, ] 
# Pesticides significant in the null model and after independent SES adjustment 
sig <- summary[summary$NULL_Model_Significant == T & summary$Pesticide_Significant_After_SES_Adjustment == 7 ,]

significant_pesticides <- sig$Pesticide


# Table S12: Pesticide use association with early-onset CRC incidence rates adjusting for socioeconomic factors
# TableS12 <- adjustment_results[adjustment_results$Pesticide %in% significant_pesticides, ]


source("Functions/pesti_pesti_adjustments.R")
# Run analysis
results_adjustmentPesti <- mutual_adjustment_analysis(DF, significant_pesticides, outcome = "value",  time_var = "variable",
                                                      group_var = "SC", filter_group = TRUE, group_min = 10,
                                                      filter_top = TRUE,  top_percent = 0.05,
                                                      filter_bottom = FALSE,    bottom_percent = 0.05)

#write.csv(results_adjustment, "Results/SEER_pesticide/Pesti_Pestiadj_SESsign.csv", row.names = F)

summary_tablepesti2 <- results_adjustmentPesti %>%
  group_by(Pesticide_2) %>%
  summarise(
    N_tests = n(),
    Made_non_significant = sum(Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj >= 0.05, na.rm = TRUE),
    Still_significant = sum(Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj < 0.05, na.rm = TRUE),
    Never_significant = sum(Pvalue_Pest1_unadj >= 0.05, na.rm = TRUE)
  ) %>%
  mutate(
    Proportion_non_significant = Never_significant / N_tests
  )

high<- c("XFENBUCONAZOLE", "XPROPARGITE", "XFLUMIOXAZIN") # Pesticides with low overlapping data points

results_adjustmentPesti2 <- results_adjustmentPesti[! results_adjustmentPesti$Pesticide_1 %in% high & !results_adjustmentPesti$Pesticide_2 %in% high, ] 


summary_tablepesti_2 <- results_adjustmentPesti2 %>%
  group_by(Pesticide_1) %>%
  summarise(
    N_tests = n(),
    Unadjusted_signif = sum(Pvalue_Pest1_unadj < 0.05, na.rm = TRUE),
    Adjusted_signif = sum(Pvalue_Pest1_adj < 0.05, na.rm = TRUE),
    Still_signif_after_adjustment = sum(Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj < 0.05, na.rm = TRUE),
    Lost_significance = sum(Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj >= 0.05, na.rm = TRUE),
    Gained_significance = sum(Pvalue_Pest1_unadj >= 0.05 & Pvalue_Pest1_adj < 0.05, na.rm = TRUE)
  )


significant_pesticides2 <- summary_tablepesti_2[summary_tablepesti_2$Adjusted_signif >10, ]$Pesticide_1


################

# Adjust the association between pesticides and eoCRC for the other pesticides and the social economic factors. 

source("Functions/pesticide_pesti_SES_adjustment.R")

results <- pesticide_interaction_adjustment(data = DF, pesticides = significant_pesticides2, ses_vars = 3:9, outcome = "value",
                                            time_var = "variable", group_var = "SC",
                                            SC_filter = TRUE, SC_n_min = 10,
                                            top_filter = TRUE, top_percent = 0.05,
                                            bottom_filter = FALSE, bottom_percent = 0.05) 

PestiSESadj_results <- results$adjustment_results  # the p-value of the pesticide with eoCRC including SES adjustment
SES_PestiSESadj_results <- results$SES_results                # the p-value of SES variable in the model testing the pesticide with eoCRC
summary_PestiSESadj <- results$summary_table                  # which pesticides are significant in the NULL model and in how many stay significant after adjustments  


# write.csv(PestiSESadj_results, "Results/SEER_pesticide/AllPestiAdjSESpesti.csv", row.names = F)
# write.csv(SES_PestiSESadj_results, "Results/SEER_pesticide/SESinAdjAllPestiSES.csv", row.names = F)
# write.csv(summary_PestiSESadj, "Results/SEER_pesticide/Summary_allPesti_SESpesti.csv", row.names = F)



# Table S13: Pesticide use association with early-onset CRC incidence rates adjusting for other pesticide usage and Socioeconomic factors
TableS13 <- PestiSESadj_results 


# Define SES variables
ses_vars <- c("Education.index", "Income", "home_value", "rent", "Poverty_150", "Unemployed", "Working_class")

# Filter relevant rows
df <- PestiSESadj_results
df_pesti2 <- df[df$SES_Adjusted == "NULL+Pesti2", ]
df_ses <- df[df$SES_Adjusted %in% ses_vars, ]

# Flag significant (P < 0.05) pesti2 adjustments
df_pesti2$signif_pesti2 <- df_pesti2$P_value < 0.05

# Build summary per Pesticide1
summary_list <- lapply(unique(df$Pesticide1), function(pest1) {
  df_pest1 <- df_pesti2[df_pesti2$Pesticide1 == pest1, ]
  total_tests <- nrow(df_pest1)
  signif_tests <- sum(df_pest1$signif_pesti2)
  
  # For each significant Pesticide2, count how often each SES is significant
  ses_counts <- sapply(ses_vars, function(ses) {
    sum(df_ses$Pesticide1 == pest1 &
          df_ses$Pesticide2 %in% df_pest1$Pesticide2[df_pest1$signif_pesti2] &
          df_ses$SES_Adjusted == ses &
          df_ses$P_value < 0.05)
  })
  
  # Combine results
  c(Pesticide1 = pest1,
    n_pesti2_tests = total_tests,
    n_signif_pesti2 = signif_tests,
    setNames(ses_counts, paste0("SES_", ses_vars, "_count")))
})

# Combine list into data frame
summary_df <- as.data.frame(do.call(rbind, summary_list))

# Convert numeric columns
summary_df[, -1] <- lapply(summary_df[, -1], as.numeric)

summary_df

#
#
#  Figure 6: Pesticide use intensity is associated with early-onset colorectal cancer incidence
#
#

# Load packages
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(stringr)
library(tibble) 


## Organize SES results

adjustment_results <- read.csv("Results/SEER_pesticide/AllPestiAdjSES.csv", header=T)

sig_pestis <- adjustment_results %>%
  filter(SES_Adjusted == "NULL", P_value < 0.05) %>%
  pull(Pesticide) %>%
  unique()
length(sig_pestis)


adjustment_results$logP <- -log10(adjustment_results$P_value)


# Filter only SES-adjusted models (exclude NULL)
ses_mat <- adjustment_results %>%
  filter(Pesticide %in% sig_pestis, SES_Adjusted != "NULL") %>%
  select(Pesticide, SES_Adjusted, logP) %>%
  pivot_wider(names_from = SES_Adjusted, values_from = logP)

ses_mat$Combined_SES <- NULL

ses_mat2 <- as.data.frame(ses_mat)
rownames(ses_mat2) <- ses_mat2$Pesticide
ses_mat2 <- ses_mat2[ , -1]


# threshold for significance: -log10(0.05) ≈ 1.3
threshold <- 1.3
# Count number of significant SES adjustments per pesticide
sig_counts <- rowSums(ses_mat2 > threshold, na.rm = TRUE)

# Overall strength (sum of -log10(p))
sig_strength <- rowSums(ses_mat2, na.rm = TRUE)

# Create ranking: first by count, then by strength
row_order <- order(sig_counts, sig_strength, decreasing = TRUE)

# Identify rows where all variables < threshold
bold_rows <- rownames(ses_mat2)[rowSums(ses_mat2 > threshold, na.rm = TRUE) == ncol(ses_mat2)]

# Set row name styles dynamically
row_fontfaces <- ifelse(rownames(ses_mat2) %in% bold_rows, "bold", "plain")

# Define colors: p >= 0.05 (<=1.3) = gray, p < 0.05 (>1.3) = blue gradient by strength
col_fun <- colorRamp2(c(0, 1.2999999999999999, threshold, max(ses_mat2, na.rm = TRUE)),
                      c("#f2eeee", "#f2eeee", "#c2cae5", "#536cb7"))

rownames(ses_mat2) <- str_remove(rownames(ses_mat2), "X")
rownames(ses_mat2) <- str_to_title(rownames(ses_mat2))
names(ses_mat2) <- c("Education index",   "Median household income", "Median house value", "Median rent", "Below 150% poverty (%)", "Unemployed (%)", "Working class (%)")

lgd = Legend(col_fun = col_fun, title = "-log10(p)",
             labels_gp = gpar(fontsize = 6),
             title_gp  = gpar(fontsize = 8, fontface = "bold"),
             grid_width = unit(0.25, "cm"))


# Plot heatmap with custom ordering
ht<-Heatmap(as.matrix(ses_mat2),
            name = "-log10(p)",
            col = col_fun,
            row_order = row_order,
            cluster_columns = TRUE,
            cluster_rows = FALSE,   # disable row clustering, since we use our own order
            rect_gp = gpar(col = "white"),
            row_names_gp = gpar(fontsize = 7, fontface = row_fontfaces),
            column_names_gp = gpar(fontsize = 7),
            column_names_rot = 75, show_column_dend = FALSE,
            column_title = "County-level socioeconomic status indicators",
            column_title_side = "bottom",
            column_title_gp = gpar(fontsize = 7),
            show_heatmap_legend = F)

#pdf("Fig6A.pdf", width = 3.3 , height = 8.2)

draw(ht, show_heatmap_legend = FALSE, newpage = TRUE)
draw(lgd, x = unit(0.75, "npc"), y = unit(0.1, "npc"), just = c("left", "center"))

#dev.off()



## Organize pesticide adjustment results

df <- read.csv("Results/SEER_pesticide/Pesti_Pestiadj_SESsign.csv", header=T)
df$Pvalues_Pest1 <- -log10(df$Pvalue_Pest1_adj)

df <- df %>%
  mutate(
    Pvalues_Pest1_mod = case_when(
      Pvalue_Pest1_unadj >= 0.05 ~ -20,  # overwrite with -20
      Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj >= 0.05 ~ -10,  # overwrite with -10
      Pvalue_Pest1_unadj < 0.05 & Pvalue_Pest1_adj < 0.05 ~ Pvalues_Pest1  # keep original -log10 value
    )
  )

mat_logp <-  as.data.frame(df %>%
                             select(Pesticide_1, Pesticide_2, Pvalues_Pest1_mod) %>%
                             pivot_wider(names_from = Pesticide_2, values_from = Pvalues_Pest1_mod))

rownames(mat_logp) <- mat_logp$Pesticide_1
mat_logp$Pesticide_1 <- NULL
threshold <- 1.3

# Loop over all cells and set structural NA (row == col) to "white"
for(i in seq_len(nrow(mat_logp))){
  for(j in seq_len(ncol(mat_logp))){
    if(rownames(mat_logp)[i] == colnames(mat_logp)[j] && is.na(mat_logp[i,j])){
      mat_logp[i,j] <- -7
    }
  }
}

mat_logp[is.na(mat_logp)] <- -5
rownames(mat_logp) <- str_remove(rownames(mat_logp), "X")
colnames(mat_logp) <- str_remove(colnames(mat_logp), "X")
rownames(mat_logp) <- str_to_title(rownames(mat_logp))
colnames(mat_logp) <- str_to_title(colnames(mat_logp))
ordered_names <- sort(rownames(mat_logp))

# Apply to both rows and columns
mat_logp <- mat_logp[ordered_names, ordered_names]


# -20 = loss null model
# -10 = Pesticide adjustment
# -7 = same pesti
# -5 = Not tested


col_fun <- colorRamp2(c(-20, -10, -7,-5, threshold, max(mat_logp, na.rm = TRUE)),
                      c("#dc94ac", "#f9edf1","#ae3960",  "#c1afb5", "#c2cae5", "#536cb7"))


ht<-
  Heatmap(as.matrix(mat_logp),
          name = "-log10(p)",
          col = col_fun,
          #    row_order = common_order,
          cluster_columns = F,
          cluster_rows = F,   # disable row clustering, since we use our own order
          rect_gp = gpar(col = "white"),
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 7),
          column_names_rot = 75, show_column_dend = FALSE,
          column_title = "Pesticide adjustment",
          column_title_side = "bottom",
          column_title_gp = gpar(fontsize = 7),
          show_heatmap_legend = F)

lgd_discrete <- Legend(labels = c("No overlap", "Basic model", "After adjustment"),
                       labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 6),
                       title = "Loss of significance",  
                       border = "black", 
                       grid_height = unit(2, "mm"), grid_width = unit(2, "mm"),
                       legend_gp = gpar(fill = c("#c1afb5", "#dc94ac", "#f9edf1")))

# Continuous legend for significant values >= 1.3
lgd_cont <- Legend(title = "-log10(Pvalue)",
                   col_fun = colorRamp2(c(1.3, max(mat_logp, na.rm = TRUE)), c("#c2cae5", "#536cb7")),
                   at = c(1.3, round(max(mat_logp, na.rm = TRUE), 1)), direction = "horizontal",
                   labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 6),border = "black", 
                   grid_height =  unit(2, "mm"),
                   labels = c("1.3", round(max(mat_logp, na.rm = TRUE), 1)))

#pdf("Fig6B.pdf", width = 3.4 , height = 4)

draw(ht, show_heatmap_legend = FALSE, newpage = TRUE, padding = unit(c(0.5, 3, 0.5, -0.50), "mm"))
draw(lgd_discrete, x = unit(0.75, "npc"), y = unit(0.06, "npc"), just = c("left", "center"))
draw(lgd_cont, x = unit(0.75, "npc"), y = unit(0.175, "npc"), just = c("left", "center"))

#dev.off()


## Organize top 5 results


PestiSESadj_results <-read.csv("Results/SEER_pesticide/AllPestiAdjSESpesti.csv", header=T)

pes <- c("XATRAZINE", "XESFENVALERATE", "XGLYPHOSATE", "XNICOSULFURON", "XPICLORAM")

sub <- PestiSESadj_results[PestiSESadj_results$Pesticide1 %in% pes & PestiSESadj_results$Pesticide2 %in% pes, ]

df_plot <- sub %>%
  filter(SES_Adjusted != "NULL") %>%
  mutate(SES_Adjusted = factor(SES_Adjusted,
                               levels = c("NULL+Pesti2",  "Education.index", "Income", "home_value", "rent", "Poverty_150",  "Unemployed", "Working_class"))) %>%
  select(Pesticide1, Pesticide2, SES_Adjusted, Effect, P_value, Significant)


df_plot <- df_plot[complete.cases(df_plot$SES_Adjusted), ]

df_plot_mod <- df_plot %>%
  mutate(P_value = as.numeric(P_value)) %>%
  group_by(Pesticide1, Pesticide2) %>%
  mutate(
    # extract the NULL+Pesti2 p-value for this pesticide pair
    null_p = P_value[SES_Adjusted == "NULL+Pesti2"],
    
    # new P_value2 according to rules
    P_value2 = case_when(
      SES_Adjusted == "NULL+Pesti2" & null_p > 0.05 ~ -2,
      SES_Adjusted == "NULL+Pesti2" & null_p <= 0.05 ~ -4,
      null_p > 0.05 ~ -5,
      null_p <= 0.05 & P_value < 0.05 ~ -10,
      null_p <= 0.05 & P_value >= 0.05 ~ -15
    )
  ) %>%
  ungroup() %>%
  select(-null_p)

df_plot_mod$Pesticide1 <- str_remove(df_plot_mod$Pesticide1, "X")
df_plot_mod$Pesticide2 <- str_remove(df_plot_mod$Pesticide2, "X")

df_plot_mod$Pesticide1 <- str_to_title(df_plot_mod$Pesticide1)
df_plot_mod$Pesticide2 <- str_to_title(df_plot_mod$Pesticide2)

df_plot_mod$Pesticide1 <- gsub("Esfenvalerate", "Esfenval", df_plot_mod$Pesticide1)
df_plot_mod$Pesticide2 <- gsub("Esfenvalerate", "Esfenval", df_plot_mod$Pesticide2)

df_plot_mod$Pesticide1 <- gsub("Nicosulfuron", "Nicosulf", df_plot_mod$Pesticide1)
df_plot_mod$Pesticide2 <- gsub("Nicosulfuron", "Nicosulf", df_plot_mod$Pesticide2)

rename_map <- c(
  "NULL+Pesti2" = "Pesticide",
  "Education.index"= "Education index",
  "Income"= "Household income",
  "home_value" = "House value",     
  "rent"  =  "Rent" ,
  "Poverty_150" = "150% poverty (%)",    
  "Unemployed"   = "Unemployed (%)",
  "Working_class"= "Working class (%)"
)

df_plot_mod$SES_Adjusted <- unname(rename_map[df_plot_mod$SES_Adjusted])
df_plot_mod$Adjusted <- paste(df_plot_mod$Pesticide1, df_plot_mod$Pesticide2, sep= "+")

mat_logp <-  as.data.frame(df_plot_mod %>%
                             select(P_value2, Adjusted, SES_Adjusted) %>%
                             pivot_wider(names_from = Adjusted, values_from = P_value2))

rownames(mat_logp) <- mat_logp$Pesticide1
rownames(mat_logp) <- mat_logp$SES_Adjusted
mat_logp$SES_Adjusted <- NULL

text_list <- sub("\\+.*", "", colnames(mat_logp))


# Row annotation for pesticide1
ha <- rowAnnotation(Pesticide1 = anno_block(gp = gpar(fill= "#FEE6CE", col="white"),
                                            labels = unique(text_list),
                                            labels_gp = gpar(col = "black", fontsize=7)     ))

names(mat_logp) <- sub(".*\\+", "", names(mat_logp))
wide <- mat_logp
over <- data.matrix(wide)

col_fun <- colorRamp2(c(-2, -4, -5, -10, -15),
                      c("#dc94ac", "#536cb7","#efcfda",  "#7c8fc8", "#c2cae5"))

lgd_discrete <- Legend(labels = c("Pesticide n.s.",  "SES not tested", "Pesticide p<0.05", "Pesticide & SES p<0.05", "Pesticide & SES n.s."),
                       title = "Association eoCRC", 
                       title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5), 
                       grid_height = unit(2, "mm"), grid_width = unit(2, "mm"),
                       legend_gp = gpar(fill = c("#dc94ac", "#efcfda",  "#536cb7","#7c8fc8", "#c2cae5")))


ht<-Heatmap(data.matrix(t(over)), 
            show_row_names = T, 
            show_row_dend = FALSE, 
            show_column_dend = FALSE,
            show_column_names = T,
            show_parent_dend_line = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, # row_title=NULL, 
            col = col_fun,  rect_gp = gpar(col = "white", lwd = 0.5),
            show_heatmap_legend = F,
            left_annotation = ha ,
            row_split = text_list, 
            row_names_gp = gpar(fontsize = 7),  
            row_title = "Pesticide adjustment", 
            row_title_side = "right",
            row_title_gp = gpar(fontsize = 7),  
            column_gap=unit(.0125, "npc"),
            column_names_gp = gpar(fontsize = 7),
            column_names_rot = 75,
            column_title = "Adjustment pesticide & socioeconomic indicators", 
            column_title_side = "bottom",
            column_title_gp = gpar(fontsize = 7),)

#pdf("Fig6C.pdf", width = 3.4 , height = 4)

draw(ht, show_heatmap_legend = FALSE, newpage = TRUE, padding = unit(c(1, 1, 1, 5), "mm"))
draw(lgd_discrete, x = unit(0.7, "npc"), y = unit(0.13, "npc"), just = c("left", "center"))

#dev.off()