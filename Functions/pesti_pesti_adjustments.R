### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


# Assess whether significant pesticides lose their association when adjusted for each other 


# Function to test mutual adjustment of pesticides
library(lmerTest)
library(dplyr)
library(tidyr)
library(progress)

mutual_adjustment_analysis <- function(data,
                                       significant_pesticides,
                                       outcome = "value",
                                       time_var = "variable",
                                       group_var = "SC",
                                       filter_group = TRUE,
                                       group_min = 10,
                                       filter_top = TRUE,
                                       top_percent = 0.05,
                                       filter_bottom = FALSE,
                                       bottom_percent = 0.05) {
  
  results <- data.frame()
  
  pb <- progress_bar$new(
    total = (length(significant_pesticides) * (length(significant_pesticides) - 1)),
    format = " Adjusting [:bar] :percent eta: :eta"
  )
  
  for (i in seq_along(significant_pesticides)) {
    pest1 <- significant_pesticides[i]
    
    for (j in setdiff(significant_pesticides, pest1)) {
      pb$tick()
      pest2 <- j
      
      # Base selection for pesti1
      df_pest1 <- data %>%
        select(all_of(c(outcome, time_var, group_var, pest1))) %>%
        filter(complete.cases(.))
      
      # 1. Group filtering
      if (filter_group) {
        df_pest1 <- df_pest1 %>%
          group_by(.data[[group_var]]) %>%
          filter(n() >= group_min) %>%
          ungroup()
      }
      
      # 2. Top percentile filter
      if (filter_top) {
        df_pest1 <- df_pest1 %>%
          filter(.data[[pest1]] <= quantile(.data[[pest1]], 1 - top_percent, na.rm = TRUE))
      }
      
      # 3. Bottom percentile filter
      if (filter_bottom) {
        df_pest1 <- df_pest1 %>%
          filter(.data[[pest1]] >= quantile(.data[[pest1]], bottom_percent, na.rm = TRUE))
      }
      
      # Skip if not enough data
      if (nrow(df_pest1) < 50) next
      
      # Merge in pest2 from full data
      df_complete <- df_pest1 %>%
        left_join(data %>% select(all_of(c(group_var, time_var, pest2))), by = c(group_var, time_var)) %>%
        filter(complete.cases(.))
      
      if (nrow(df_complete) < 50) next
      
      # Mutual adjustment model
      model_adj <- lmerTest::lmer(
        as.formula(paste(outcome, "~", pest1, "+", pest2, "+", time_var, "+ (1|", group_var, ")")),
        data = df_complete
      )
      sum_adj <- summary(model_adj)$coefficients
      
      # Correlation between pesticides
      cor_test <- cor.test(df_complete[[pest1]], df_complete[[pest2]])
      pest_cor <- cor_test$estimate
      pest_cor_pval <- cor_test$p.value
      
      # Unadjusted model for pest1
      model_pest1 <- lmerTest::lmer(
        as.formula(paste(outcome, "~", pest1, "+", time_var, "+ (1|", group_var, ")")),
        data = df_complete
      )
      sum_pest1 <- summary(model_pest1)$coefficients
      
      # Unadjusted model for pest2
      
        model_pest2 <- lmerTest::lmer(
          as.formula(paste(outcome, "~", pest2, "+", time_var, "+ (1|", group_var, ")")),
          data = df_complete
        )
        sum_pest2 <- summary(model_pest2)$coefficients
        
        pest2_est <- sum_pest2[pest2, "Estimate"]
        pest2_pval <- sum_pest2[pest2, "Pr(>|t|)"]
        n_SC_pest2 <- n_distinct(df_complete[[group_var]])
      
      
      # Describe filters
      filter_desc <- paste0("GroupFilter=", filter_group,
                            ";MinGroup=", group_min,
                            ";Top=", filter_top,
                            ";TopPct=", top_percent,
                            ";Bottom=", filter_bottom,
                            ";BottomPct=", bottom_percent)
      
      # Collect results
      results <- bind_rows(results, data.frame(
        Pesticide_1 = pest1,
        Pesticide_2 = pest2,
        Estimate_Pest1_adj = sum_adj[pest1, "Estimate"],
        Pvalue_Pest1_adj = sum_adj[pest1, "Pr(>|t|)"],
        Estimate_Pest2_adj = sum_adj[pest2, "Estimate"],
        Pvalue_Pest2_adj = sum_adj[pest2, "Pr(>|t|)"],
        Correlation = pest_cor,
        Correlation_pvalue = pest_cor_pval,
        N_obs = nrow(df_complete),
        N_SC = n_distinct(df_complete[[group_var]]),
        Estimate_Pest1_unadj = sum_pest1[pest1, "Estimate"],
        Pvalue_Pest1_unadj = sum_pest1[pest1, "Pr(>|t|)"],
        N_SC_Pest1_unadj = n_distinct(df_pest1[[group_var]]),
        Estimate_Pest2_unadj = pest2_est,
        Pvalue_Pest2_unadj = pest2_pval,
        N_SC_Pest2_unadj = n_SC_pest2,
        Filters = filter_desc
      ))
    }
  }
  
  return(results)
}
