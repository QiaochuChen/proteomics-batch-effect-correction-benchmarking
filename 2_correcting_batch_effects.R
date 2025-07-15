# ' Title: Data cleaning
# ' Author: Qiaochu Chen
# ‘ Date: Jun 30th, 2025

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


## global. ------------------
if (!interactive()) pboptions(type = "none")
# options(mc.cores = min(24, max(1, parallel::detectCores() - 8)))
options(mc.cores = 4)
source("./utils/BEC.R")
set.seed(20250630)


## correct batch effects at precursor/peptide/protein levels. 约30秒钟 ---------
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("expdata_log", all_files)]

correct_methods <- c("ratio_by_d6", "median_centering", "combat", "harmony", "WaveICA", "RUV-III-C")
label_names <- c("ratio", "median", "combat", "harmony", "waveica", "ruv")

corrected_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "correcting ", file_id, "\n"))
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded)")
  file_name <- str_extract(file_id, "^.*(?=expdata)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  df_quant_wide_i <- fread(file_id)
  
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  
  if (dataset %in% "simulated") {
    df_features <- df_quant_wide_i %>%
      tibble::rowid_to_column("feature") %>%
      distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
    
    df_mapping <- fread("./data/expfiles/simulated/mapping.csv")
    
    ctrl_features <- df_mapping %>%
      filter(class %in% 3:4) %>%
      left_join(., df_features) %>%
      pull(feature) %>%
      unique
    
    df_quant_wide_i <- df_quant_wide_i %>%
      tibble::rowid_to_column("feature") %>%
      dplyr::select(feature, all_of(df_meta$run_id)) %>%
      mutate(across(all_of(df_meta$run_id), ~ ifelse(is.na(.)|. < 0, 0, exp(.))))
    
    reference_sample <- "group1"
    
  } else {
    df_anova_i <- fread(paste(file_name, "anova_twoway_test.csv", sep = ""))
    
    ctrl_features <- df_anova_i %>%
      reshape2::dcast(., feature ~ Effect, value.var = "p<.05") %>%
      filter(`sample` %in% "" & `sample:batch` %in% "") %>%
      pull(feature)
    
    rm(df_anova_i)
    gc()
    
    df_meta <- df_meta %>%
      mutate(batch = paste(mode, lab, sep = "_")) %>%
      mutate(order = 1:nrow(df_meta))
    
    df_features <- df_quant_wide_i %>%
      tibble::rowid_to_column("feature") %>%
      distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
    
    df_quant_wide_i <- df_quant_wide_i %>%
      tibble::rowid_to_column("feature") %>%
      dplyr::select(feature, all_of(df_meta$run_id)) %>%
      mutate(across(all_of(df_meta$run_id), ~ ifelse(is.na(.)|. < 0, 0, exp(.))))
    
    reference_sample <- "D6"
  }
  
  corrected_tables_i <- mclapply(1:5, function(j) {
    
    if (!correct_methods[j] %in% "RUV-III-C") ctrl_features <- NULL

    null_device <- if (.Platform$OS.type == "unix") "/dev/null" else "NUL"
    con <- file(null_device, open = "w")
    sink(con, type = "output")
    sink(con, type = "message")
    df_corrected_j <- correct_by_x(df_quant_wide_i, df_meta, correct_methods[j],
                                   ctrl_features, reference_sample)
    sink(type = "message")
    sink(type = "output")
    close(con)
    
    rm(df_quant_wide_i, df_meta)
    gc()
    
    df_corrected_j_final <- df_corrected_j %>%
      mutate_at("feature", as.integer) %>%
      left_join(., df_features, by = "feature") %>%
      dplyr::select(!feature) %>%
      dplyr::select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
    
    file_id2save <- paste(file_name, "expdata_", label_names[j], ".csv", sep = "")
    fwrite(df_corrected_j_final, file = file_id2save)
    
    gc()
    cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))

  }, mc.cores = 5)
  
  ## RUV-III-C单独运行
  df_corrected_j <- correct_by_x(df_quant_wide_i, df_meta, "RUV-III-C", ctrl_features)
  rm(df_quant_wide_i, df_meta)
  gc()
  
  df_corrected_j_final <- df_corrected_j %>%
    mutate_at("feature", as.integer) %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
    
  file_id2save <- paste(file_name, "expdata_ruv.csv", sep = "")
  fwrite(df_corrected_j_final, file = file_id2save)
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))

  cat(paste0("[", Sys.time(), "] ", "Successfully corrected.\n"))
})


