### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


pesticide_interaction_adjustment <- function(data, pesticides, ses_vars = 3:10, outcome = "value",
                                               time_var = "variable", group_var = "SC",
                                               SC_filter = TRUE, SC_n_min = 10,
                                               top_filter = TRUE, top_percent = 0.05,
                                               bottom_filter = FALSE, bottom_percent = 0.05) {
  
  adjustment_results <- data.frame()
  SES_results <- data.frame()
  summary_table <- data.frame()
  
  pb <- progress_bar$new(
    total = length(pesticides),
    format = " Processing [:bar] :percent eta: :eta"
  )
  
  for (pesti1 in pesticides) {
    pb$tick()
    
    # Select relevant columns
    df_pest <- data %>%
      select(all_of(c(outcome, time_var, group_var, pesti1, names(data)[ses_vars]))) %>%
      filter(complete.cases(.))
    
    filter_description <- c()
    
    # Group filtering
    if (SC_filter) {
      df_pest <- df_pest %>%
        group_by(.data[[group_var]]) %>%
        filter(n() >= SC_n_min) %>%
        ungroup()
      filter_description <- c(filter_description, paste0("SC≥", SC_n_min))
    } else {
      filter_description <- c(filter_description, "No_SC_Filter")
    }
    
    
    # Top filtering
    if (top_filter) {
      q_high <- quantile(df_pest[[pesti1]], 1 - top_percent, na.rm = TRUE)
      df_pest <- df_pest %>% filter(.data[[pesti1]] <= q_high)
      filter_description <- c(filter_description, paste0("Top", top_percent * 100, "%"))
    }
    
    # Bottom filtering
    if (bottom_filter) {
      q_low <- quantile(df_pest[[pesti1]], bottom_percent, na.rm = TRUE)
      df_pest <- df_pest %>% filter(.data[[pesti1]] >= q_low)
      filter_description <- c(filter_description, paste0("Bottom", bottom_percent * 100, "%"))
    }    
    filter_desc <- paste(filter_description, collapse = "|")
    
    for (pesti2 in setdiff(pesticides, pesti1)) {
      
      # Add pesti2 to filtered df_pest
      df_complete <- df_pest %>%
        left_join(data %>% select(all_of(c(group_var, time_var, pesti2))), by = c(group_var, time_var)) %>%
        filter(complete.cases(.))
      
      if (nrow(df_complete) < 50) next
      
      # Null Model
      null_model <- lmerTest::lmer(as.formula(paste(outcome, "~", pesti1, "+", time_var, "+ (1|", group_var, ")")), data = df_complete)
      null_summary <- summary(null_model)$coefficients
      
      adjustment_results <- bind_rows(adjustment_results, data.frame(
        Pesticide1 = pesti1,
        Pesticide2 = pesti2,
        SES_Adjusted = "NULL",
        Effect = null_summary[pesti1, "Estimate"],
        P_value = null_summary[pesti1, "Pr(>|t|)"],
        Significant = null_summary[pesti1, "Pr(>|t|)"] < 0.05,
        Included_SES = "None",
        N_obs = nrow(df_complete),
        N_years = n_distinct(df_complete[[time_var]]),
        N_SC = n_distinct(df_complete[[group_var]]),
        Filter_Description = filter_desc
      ))
      
      # Null + Pesti2
      null_p2_model <- lmerTest::lmer(as.formula(paste(outcome, "~", pesti1, "+", pesti2, "+", time_var, "+ (1|", group_var, ")")), data = df_complete)
      null_p2_summary <- summary(null_p2_model)$coefficients
      
      adjustment_results <- bind_rows(adjustment_results, data.frame(
        Pesticide1 = pesti1,
        Pesticide2 = pesti2,
        SES_Adjusted = "NULL+Pesti2",
        Effect = null_p2_summary[pesti1, "Estimate"],
        P_value = null_p2_summary[pesti1, "Pr(>|t|)"],
        Significant = null_p2_summary[pesti1, "Pr(>|t|)"] < 0.05,
        Included_SES = "None",
        N_obs = nrow(df_complete),
        N_years = n_distinct(df_complete[[time_var]]),
        N_SC = n_distinct(df_complete[[group_var]]),
        Filter_Description = filter_desc
      ))
      
      # Individual SES
      ses_pval_list <- list()
      ses_sig_count <- 0
      ses_adj_sig_count <- 0
      
      for (ses_var in names(data)[ses_vars]) {
        model_formula <- paste(outcome, "~", pesti1, "+", pesti2, "+", time_var, "+", ses_var, "+ (1|", group_var, ")")
        model <- lmerTest::lmer(as.formula(model_formula), data = df_complete)
        model_sum <- summary(model)$coefficients
        
        ses_pval_list[[ses_var]] <- model_sum[ses_var, "Pr(>|t|)"]
        
        adjustment_results <- bind_rows(adjustment_results, data.frame(
          Pesticide1 = pesti1,
          Pesticide2 = pesti2,
          SES_Adjusted = ses_var,
          Effect = model_sum[pesti1, "Estimate"],
          P_value = model_sum[pesti1, "Pr(>|t|)"],
          Significant = model_sum[pesti1, "Pr(>|t|)"] < 0.05,
          Included_SES = ses_var,
          N_obs = nrow(df_complete),
          N_years = n_distinct(df_complete[[time_var]]),
          N_SC = n_distinct(df_complete[[group_var]]),
          Filter_Description = filter_desc
        ))
        
        SES_results <- bind_rows(SES_results, data.frame(
          Pesticide1 = pesti1,
          Pesticide2 = pesti2,
          SES_Variable = ses_var,
          SES_P_value = model_sum[ses_var, "Pr(>|t|)"],
          Significant = model_sum[ses_var, "Pr(>|t|)"] < 0.05
        ))
        
        if (model_sum[ses_var, "Pr(>|t|)"] < 0.05) ses_sig_count <- ses_sig_count + 1
        if (model_sum[pesti1, "Pr(>|t|)"] < 0.05) ses_adj_sig_count <- ses_adj_sig_count + 1
      }
      
      # Combined SES Model
      significant_SES <- names(ses_pval_list)[unlist(ses_pval_list) < 0.05]
      
      if (length(significant_SES) > 0) {
        combined_formula <- paste(outcome, "~", pesti1, "+", pesti2, "+", time_var, "+", paste(significant_SES, collapse = " + "), "+ (1|", group_var, ")")
        combined_model <- lmerTest::lmer(as.formula(combined_formula), data = df_complete)
        combined_sum <- summary(combined_model)$coefficients
        
        adjustment_results <- bind_rows(adjustment_results, data.frame(
          Pesticide1 = pesti1,
          Pesticide2 = pesti2,
          SES_Adjusted = "Combined_SES",
          Effect = combined_sum[pesti1, "Estimate"],
          P_value = combined_sum[pesti1, "Pr(>|t|)"],
          Significant = combined_sum[pesti1, "Pr(>|t|)"] < 0.05,
          Included_SES = paste(significant_SES, collapse = ", "),
          N_obs = nrow(df_complete),
          N_years = n_distinct(df_complete[[time_var]]),
          N_SC = n_distinct(df_complete[[group_var]]),
          Filter_Description = filter_desc
        ))
      }
      
      summary_table <- bind_rows(summary_table, data.frame(
        Pesticide1 = pesti1,
        Pesticide2 = pesti2,
        Significant_SES_Count = ses_sig_count,
        Pesticide1_Significant_After_SES_Adjustment = ses_adj_sig_count,
        N_obs = nrow(df_complete),
        N_years = n_distinct(df_complete[[time_var]]),
        N_SC = n_distinct(df_complete[[group_var]]),
        Filter_Description = filter_desc
      ))
    }
  }
  
  return(list(adjustment_results = adjustment_results, SES_results = SES_results, summary_table = summary_table))
}
