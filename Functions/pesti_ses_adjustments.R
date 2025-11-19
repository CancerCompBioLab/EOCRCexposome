### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


library(lmerTest)
library(dplyr)
library(tidyr)
library(progress)



# Do the SES variables have a significant influence on the association between the pesticides and eoCRC incidence? 
pesticide_ses_analysis <- function(data,
                                   pesticides,
                                   ses_vars = 3:10,
                                   outcome = "value",
                                   time_var = "variable",
                                   group_var = "SC",
                                   filter_group = TRUE,
                                   group_min = 10,
                                   filter_top = TRUE,
                                   top_percent = 0.05,
                                   filter_bottom = FALSE,
                                   bottom_percent = 0.05) {
  
  adjustment_results <- data.frame()
  SES_results <- data.frame()
  summary_table <- data.frame()
  skipped_dfs <- list()
  
  pb <- progress_bar$new(
    total = length(pesticides),
    format = " Processing [:bar] :percent eta: :eta"
  )
  
  for (pest in pesticides) {
    pb$tick()
    
    # Select and filter complete cases
    df_pest <- data %>%
      select(all_of(c(outcome, time_var, group_var, pest, names(data)[ses_vars]))) %>%
      filter(complete.cases(.))
    
    # Optional group filter
    if (filter_group) {
      df_pest <- df_pest %>%
        group_by(.data[[group_var]]) %>%
        filter(n() >= group_min) %>%
        ungroup()
    }
    
    # Optional top filter
    if (filter_top) {
      df_pest <- df_pest %>%
        filter(.data[[pest]] <= quantile(.data[[pest]], 1 - top_percent, na.rm = TRUE))
    }
    
    # Optional bottom filter
    if (filter_bottom) {
      df_pest <- df_pest %>%
        filter(.data[[pest]] >= quantile(.data[[pest]], bottom_percent, na.rm = TRUE))
    }
    
    # Skip if not enough data
    if (nrow(df_pest) < 50) {
      skipped_dfs[[paste0("pest_", pest)]] <- df_pest
      next
    }
    
    # Metadata summary
    total_obs <- nrow(df_pest)
    n_unique_group <- length(unique(df_pest[[group_var]]))
    n_unique_time <- length(unique(df_pest[[time_var]]))
    filter_desc <- paste0("GroupFilter=", filter_group,
                          ";MinGroup=", group_min,
                          ";Top=", filter_top,
                          ";TopPct=", top_percent,
                          ";Bottom=", filter_bottom,
                          ";BottomPct=", bottom_percent)
    
    # ---- Null Model ----
    null_model <- lmerTest::lmer(as.formula(paste(outcome, "~", pest, "+", time_var, "+ (1|", group_var, ")")), data = df_pest)
    null_summary <- summary(null_model)$coefficients
    null_effect <- null_summary[pest, "Estimate"]
    null_pval <- null_summary[pest, "Pr(>|t|)"]
    
    adjustment_results <- bind_rows(adjustment_results, data.frame(
      Pesticide = pest,
      SES_Adjusted = "NULL",
      Effect = null_effect,
      P_value = null_pval,
      Significant = null_pval < 0.05,
      Included_SES = "",
      N = total_obs,
      N_SC = n_unique_group,
      N_Time = n_unique_time,
      Filters = filter_desc
    ))
    
    ses_pval_list <- list()
    ses_sig_count <- 0
    ses_adj_sig_count <- 0
    
    # ---- SES Adjusted Models (individual) ----
    for (ses_var in names(data)[ses_vars]) {
      model <- lmerTest::lmer(as.formula(paste(outcome, "~", pest, "+", time_var, "+", ses_var, "+ (1|", group_var, ")")), data = df_pest)
      model_sum <- summary(model)$coefficients
      
      pest_effect <- model_sum[pest, "Estimate"]
      pest_pval <- model_sum[pest, "Pr(>|t|)"]
      ses_pval <- model_sum[ses_var, "Pr(>|t|)"]
      
      adjustment_results <- bind_rows(adjustment_results, data.frame(
        Pesticide = pest,
        SES_Adjusted = ses_var,
        Effect = pest_effect,
        P_value = pest_pval,
        Significant = pest_pval < 0.05,
        Included_SES = ses_var,
        N = total_obs,
        N_SC = n_unique_group,
        N_Time = n_unique_time,
        Filters = filter_desc
      ))
      
      SES_results <- bind_rows(SES_results, data.frame(
        Pesticide = pest,
        SES_Variable = ses_var,
        SES_P_value = ses_pval,
        Significant = ses_pval < 0.05,
        N = total_obs,
        N_SC = n_unique_group,
        N_Time = n_unique_time,
        Filters = filter_desc
      ))
      
      ses_pval_list[[ses_var]] <- ses_pval
      
      if (ses_pval < 0.05) ses_sig_count <- ses_sig_count + 1
      if (pest_pval < 0.05) ses_adj_sig_count <- ses_adj_sig_count + 1
    }
    
    # ---- Combined SES Model for significant SES ----
    significant_SES <- names(ses_pval_list)[unlist(ses_pval_list) < 0.05]
    
    if (length(significant_SES) > 0) {
      combined_formula <- paste(outcome, "~", pest, "+", time_var, "+", paste(significant_SES, collapse = " + "), "+ (1|", group_var, ")")
      combined_model <- lmerTest::lmer(as.formula(combined_formula), data = df_pest)
      combined_sum <- summary(combined_model)$coefficients
      
      pest_effect_comb <- combined_sum[pest, "Estimate"]
      pest_pval_comb <- combined_sum[pest, "Pr(>|t|)"]
      
      adjustment_results <- bind_rows(adjustment_results, data.frame(
        Pesticide = pest,
        SES_Adjusted = "Combined_SES",
        Effect = pest_effect_comb,
        P_value = pest_pval_comb,
        Significant = pest_pval_comb < 0.05,
        Included_SES = paste(significant_SES, collapse = ", "),
        N = total_obs,
        N_SC = n_unique_group,
        N_Time = n_unique_time,
        Filters = filter_desc
      ))
    }
    
    # ---- Summary Table ----
    summary_table <- bind_rows(summary_table, data.frame(
      Pesticide = pest,
      NULL_Model_Significant = null_pval < 0.05,
      Significant_SES_Count = ses_sig_count,
      Pesticide_Significant_After_SES_Adjustment = ses_adj_sig_count,
      N = total_obs,
      N_SC = n_unique_group,
      N_Time = n_unique_time,
      Filters = filter_desc
    ))
  }
  
  return(list(adjustment_results = adjustment_results,
              SES_results = SES_results,
              summary_table = summary_table,
              skipped_dfs = skipped_dfs
              ))
}
