# ' Author: Qiaochu Chen
# ‘ Date: Jun 26th, 2025

library(parallel)
library(pbapply)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(diann)
library(sva)
library(harmony)
library(RUVIIIC)
library(WaveICA2.0)
library(JADE)
library(corpcor)
library(lme4)


## functions for correction. --------------------
ratio_by_d6 <- function(df_quant_wide, df_meta, reference_sample) {
  
  colnames(df_meta)[1] <- "run_id"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    match_d6 <- match_batch[which(grepl(reference_sample, samples[match_batch]))]
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(D6mean = apply(.[match_d6], 1, mean)) %>%
      mutate_at(match_batch, ~ (. - D6mean)) %>%
      dplyr::select(!D6mean)
    
  }
  
  return(df_corrected)
}

median_centering <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "run_id"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(all_median = apply(.[match_batch], 1, median, na.rm = TRUE)) %>%
      mutate_at(match_batch, ~ (. - all_median)) %>%
      dplyr::select(!all_median)
    
  }
  
  return(df_corrected)
}

mean_centering <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "run_id"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(all_mean = apply(.[match_batch], 1, mean, na.rm = TRUE)) %>%
      mutate_at(match_batch, ~ (. - all_mean)) %>%
      dplyr::select(!all_mean)
    
  }
  
  return(df_corrected)
}

combat_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ as.factor(sample), data = df_meta)
  
  df_corrected <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0, log(.))) %>%
    as.matrix() %>%
    ComBat(., df_meta$batch, data_mode) %>%
    as.data.frame() %>%
    mutate_all(~ na_if(., 0)) %>%
    rownames_to_column(feature_header)
  
  return(df_corrected)
}

sva_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  
  df_log_wide <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0.01, log(.))) %>%
    as.matrix()
  
  df_meta <- df_meta %>%
    column_to_rownames("run_id") %>%
    dplyr::select(sample, batch)
  
  data_mode <- model.matrix(~ as.factor(sample), data = df_meta)
  data_mode0 <- model.matrix(~ 1, data = df_meta)
  
  sv_num <- num.sv(df_log_wide, data_mode, method = "be")
  sv_obj <- sva(df_log_wide, data_mode, data_mode0, n.sv = sv_num)
  
  sv_mode <- cbind(data_mode, sv_obj$sv)
  hat <- solve(t(sv_mode) %*% sv_mode) %*% t(sv_mode)
  beta <- (hat %*% t(df_log_wide))
  
  num_type <- ncol(data_mode)
  sva_be <- t((sv_mode[, -c(1:num_type)]) %*% beta[-c(1:num_type), ])
  sva_corrected <- (df_log_wide - sva_be) %>%
    as.data.frame() %>%
    rownames_to_column(feature_header)
  
  sva_corrected[df_quant_wide == 0] <- NA
  
  return(sva_corrected)
  
}

harmony_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  
  df_log_wide <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0, log(.))) %>%
    as.matrix()
  
  harmony_corrected <- df_log_wide %>%
    HarmonyMatrix(., df_meta, vars_use = "batch", do_pca = FALSE) %>%
    t %>%
    as.data.frame() %>%
    rownames_to_column(feature_header)
  
  harmony_corrected[df_quant_wide == 0] <- NA
  
  return(harmony_corrected)
  
}

ruv3c_solve <- function(df_quant_wide, df_meta, ctrl_features, k=3) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ sample - 1, data = df_meta)
  
  colnames(df_quant_wide)[1] <- "feature"
  vary_features <- setdiff(df_quant_wide$feature, ctrl_features)
  
  df_corrected <- df_quant_wide %>%
    column_to_rownames("feature") %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, NA, log(.))) %>%
    t %>%
    RUVIII_C(k, ., data_mode,
             toCorrect = vary_features,
             controls = ctrl_features) %>%
    t %>%
    as.data.frame %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    rownames_to_column(feature_header)
  
  return(df_corrected)
}

waveica2_solve <- function(df_quant_wide, df_meta, alpha=0, cutoff=0.1, k=10) {
  
  colnames(df_quant_wide)[1] <- "feature"
  colnames(df_meta)[1] <- "run_id"
  
  df_quant_wide_t <- df_quant_wide %>%
    mutate_at(df_meta$run_id, ~ ifelse(. == 0, 0, log(.))) %>%
    reshape2::melt(., id = 1, variable.name = "run_id") %>%
    left_join(., df_meta, by = "run_id") %>%
    reshape2::dcast(., order + sample + batch ~ feature, value.var = "value") %>%
    mutate_at("batch", ~ as.numeric(factor(.)))
  
  df_corrected <- WaveICA_2.0(data = df_quant_wide_t[, -(1:3)], wf = "haar",
                              Injection_Order = df_quant_wide_t$order,
                              Cutoff = cutoff, alpha = alpha, K = k)
  
  df_corrected_final <- df_corrected$data_wave %>%
    as.data.frame() %>%
    tibble::rowid_to_column("order") %>%
    left_join(., df_meta, by = "order") %>%
    dplyr::select(run_id, 2:(nrow(df_quant_wide) + 1)) %>%
    reshape2::melt(., id = 1, variable.name = "feature") %>%
    mutate_at("feature", as.character) %>%
    mutate_at("feature", ~ ifelse(grepl("^Var", .), "", .)) %>%
    reshape2::dcast(., feature ~ run_id, value.var = "value")
  
  df_corrected_final[df_quant_wide == 0] <- NA
  
  return(df_corrected_final)
}

correct_by_loess <- function(y, x) {
  
  data_input <- data.frame(value = y, order = x)
  
  values <- data_input$value
  orders <- data_input$order
  
  if (length(values) > 2) {
    
    drift_loess <- loess(
      formula = value ~ order,
      data = data_input,
      control = loess.control(surface = "direct"),
      degree = 2,
      normalize = FALSE)
    
    drift <- predict(drift_loess, orders)
    values_final <- values - drift + median(values)
    
  } else {
    
    values_final <- values
    
  }
  
  return(values_final)
}

loess_solve <- function(df_quant_wide, df_meta) {
  
  colnames(df_quant_wide)[1] <- "feature"
  colnames(df_meta)[1] <- "run_id"
  feature_n <- nrow(df_quant_wide)
  
  batches <- unique(df_meta$batch)
  df_quant_wide <- df_quant_wide %>%
    mutate_at("feature", as.character) %>%
    mutate_if(is.numeric, ~ ifelse(. == 0, NA, log(.)))
  
  expr_quantile <- df_quant_wide %>%
    column_to_rownames("feature") %>%
    limma::normalizeQuantiles(., ties = T)
  
  loess_tables <- mclapply(1:feature_n, function(i) {
    
    loess_tables_i <- mclapply(batches, function(batch_id) {
      
      data_test_j <- df_quant_wide[i, ] %>%
        reshape2::melt(., id.vars = 1, verbose = FALSE) %>%
        as.data.frame %>%
        na.omit %>%
        dplyr::rename(run_id = variable) %>%
        left_join(., df_meta, by = "run_id") %>%
        filter(batch %in% batch_id) %>%
        mutate(value_loess = correct_by_loess(value, order))
      
      return(data_test_j)
    }, mc.cores = 2)
    
    data_test_i <- rbindlist(loess_tables_i)
    
    return(data_test_i)
  }, mc.cores = 2)
  
  data_loess <- data.table::rbindlist(loess_tables)
  
  expr_loess <- reshape2::dcast(data_loess, feature ~ run_id, value.var = "value_loess")
  
  return(expr_loess)
}

correct_by_x <- function(df_quant_wide, df_meta,
                         method = c("ratio_by_d6",
                                    "median_centering",
                                    "mean_centering",
                                    "combat",
                                    "sva",
                                    "harmony",
                                    "RUV-III-C",
                                    "WaveICA"),
                         ctrl_features = NULL,
                         reference_sample = "D6") {
  
  if (method == "ratio_by_d6") {
    df_corrected <- ratio_by_d6(df_quant_wide, df_meta, reference_sample)
  } else if (method == "median_centering") {
    df_corrected <- median_centering(df_quant_wide, df_meta)
  } else if (method == "mean_centering") {
    df_corrected <- mean_centering(df_quant_wide, df_meta)
  } else if (method == "combat") {
    df_corrected <- combat_solve(df_quant_wide, df_meta)
  } else if(method == "sva") {
    df_corrected <- sva_solve(df_quant_wide, df_meta)
  } else if(method == "harmony") {
    df_corrected <- harmony_solve(df_quant_wide, df_meta)
  } else if(method == "RUV-III-C") {
    df_corrected <- ruv3c_solve(df_quant_wide, df_meta, ctrl_features)
  } else if(method == "WaveICA") {
    df_corrected <- waveica2_solve(df_quant_wide, df_meta)
  }
  
  return(df_corrected)
}


## functions for quantification. --------------------
maxlfq_solve <- function(quantities,
                         peptides,
                         samples,
                         margin = -10.0001) {
  
  .Call('_diann_maxlfq_solve',
        PACKAGE = 'diann',
        quantities, peptides, samples, margin)
}

quantify_by_x <- function(df_quant, df_mapping,
                          group.header = "peptide",
                          subgroup.header = "precursor",
                          sample.id.header = "raw_file",
                          intensity.header = "intensity",
                          method = c("MaxLFQ", "iBAQ", "TopPep3"), 
                          mc.cores = 4) {
  
  df_quant_tmp <- df_quant %>%
    dplyr::rename(c("sample_id_header" = sample.id.header,
                    "group_header" = group.header,
                    "subgroup_header" = subgroup.header,
                    "intensity_header" = intensity.header))
  
  all_samples <- df_quant_tmp %>% pull(sample_id_header) %>% unique %>% sort
  all_groups <- df_quant_tmp %>% pull(group_header) %>% unique
  
  if (method == "MaxLFQ") {
    
    lfq_tables <- mclapply(all_groups, function(group_id) {
      
      df_quant_i <- df_quant_tmp %>%
        filter(group_header %in% group_id) %>%
        distinct(sample_id_header, subgroup_header, group_header, intensity_header)
      
      tmp_subgroups <- df_quant_i %>% pull(subgroup_header) %>% unique
      tmp_samples <- df_quant_i %>% pull(sample_id_header) %>% unique %>% sort
      
      lfq_i <- rep(NA, length(all_samples))
      names(lfq_i) <- all_samples
      
      if (length(tmp_samples) > 1) {
        
        quant_validated_j <- df_quant_i %>%
          reshape2::acast(., subgroup_header ~ sample_id_header,
                          value.var = "intensity_header") %>%
          log
        
        quant_validated_j[is.infinite(quant_validated_j)] <- NA
        
        lfq_values <- maxlfq_solve(as.vector(quant_validated_j),
                                   nrow(quant_validated_j),
                                   ncol(quant_validated_j))
        
        names(lfq_values) <- tmp_samples
        
        lfq_i[match(names(lfq_values), names(lfq_i))] <- lfq_values
        
      }
      
      quant_lfq_i <- data.frame(group_header = group_id,
                                sample_id_header = names(lfq_i),
                                lfq = lfq_i)
      
      return(quant_lfq_i)
      
    }, mc.cores = mc.cores)
    
    df_lfq <- data.table::rbindlist(lfq_tables)
    
    df_final <- df_lfq %>%
      plyr::rename(c("group_header" = group.header,
                     "sample_id_header" = sample.id.header,
                     "lfq" = intensity.header))
    
  }
  
  if (method == "iBAQ") {
    
    ibaq_tables <- mclapply(all_groups, function(group_id) {
      
      df_quant_i <- df_quant_tmp %>%
        filter(group_header %in% group_id) %>%
        distinct(sample_id_header, subgroup_header, group_header, intensity_header)
      
      quant_ibaq_i <- df_quant_i %>%
        group_by(sample_id_header, group_header) %>%
        dplyr::summarise(ibaq = mean(intensity_header, na.rm = TRUE), .groups = "drop") %>%
        mutate_at("ibaq", ~ ifelse(. == 0, NA, log(.)))
      
      return(quant_ibaq_i)
      
    }, mc.cores = mc.cores)
    
    df_ibaq <- data.table::rbindlist(ibaq_tables)
    
    df_final <- df_ibaq %>%
      plyr::rename(c("group_header" = group.header,
                     "sample_id_header" = sample.id.header,
                     "ibaq" = intensity.header))
  }
  
  if (method == "TopPep3") {
    
    top3_tables <- mclapply(all_groups, function(group_id) {
      
      df_quant_i <- df_quant_tmp %>%
        filter(group_header %in% group_id) %>%
        distinct(sample_id_header, subgroup_header, group_header, intensity_header)
      
      df_sum_i <- df_quant_i %>%
        group_by(subgroup_header) %>%
        dplyr::summarise(total = mean(intensity_header), .groups = "drop")
      
      if (nrow(df_sum_i) >= 3) {
        
        kept_subgroups <- df_sum_i %>%
          filter(total >= sort(total, decreasing = TRUE)[3]) %>%
          pull(subgroup_header)
      } else {
        
        kept_subgroups <- df_sum_i %>%
          pull(subgroup_header)
      }
      
      quant_top3_i <- df_quant_i %>%
        filter(subgroup_header %in% kept_subgroups) %>%
        group_by(sample_id_header, group_header) %>%
        dplyr::summarise(top3 = mean(intensity_header, na.rm = TRUE), .groups = "drop") %>%
        mutate_at("top3", ~ ifelse(. == 0, NA, log(.)))
      
      return(quant_top3_i)
      
    }, mc.cores = mc.cores)
    
    df_top3 <- data.table::rbindlist(top3_tables)
    
    df_final <- df_top3 %>%
      plyr::rename(c("group_header" = group.header,
                     "sample_id_header" = sample.id.header,
                     "top3" = intensity.header))
  }
  
  return(df_final)
}


