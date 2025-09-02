# ' Title: Data cleaning
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


if (!interactive()) pboptions(type = "none")
set.seed(20250630)


## set methods ------------
source("./utils/BEC.R")
quantification_methods <- c("MaxLFQ", "iBAQ", "TopPep3")
label_names <- c("maxlfq", "ibaq", "toppep3")

if (!interactive()) pboptions(type = "none")
options(mc.cores = min(24, max(1, parallel::detectCores() - 8)))


## aggregate precursors into peptides: about 25 min -----------------
gc()
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("precursor/expdata_log", all_files)]
pep_tables <- pblapply(queried_files, function(file_id) {
  
  print(file_id)
  file_path <- str_extract(file_id, ".+(?=precursor)")
  file_name <- str_extract(file_id, "expdata.+")
  df_quant_wide_i <- fread(file_id)
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  df_quant_pre2pep <- df_quant_wide_i %>%
    dplyr::select(precursor, peptide, protein, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1:3, variable.name = "run_id", na.rm = TRUE) %>%
    mutate(peptide2protein = paste(peptide, protein, sep = "_")) %>%
    mutate(intensity = exp(value))
  
  pep_tables <- mclapply(1:3, function(j) {
    
    df_final <- quantify_by_x(df_quant = df_quant_pre2pep,
                              group.header = "peptide2protein",
                              subgroup.header = "precursor",
                              sample.id.header = "run_id",
                              intensity.header = "intensity",
                              method = quantification_methods[j],
                              mc.cores = 8)
    
    df_final_complete <- df_final %>%
      reshape2::dcast(., peptide2protein ~ run_id, value.var = "intensity") %>%
      tidyr::separate(peptide2protein, c("peptide", "protein"), sep = "_")
    
    fwrite(df_final_complete,
           file = paste(file_path, "peptide/", label_names[j], "/", file_name, sep = ""))
  }, mc.cores = 3)
})


## aggregate peptides into proteins: about 5 min -----------------
rm(pep_tables)
gc()
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(paste("peptide", sep = ""), all_files)]
queried_files <- queried_files[grepl("expdata_log", queried_files)]
pro_tables <- pblapply(queried_files, function(file_id) {
  
  print(file_id)
  file_path <- str_extract(file_id, ".+(?=peptide)")
  file_name <- str_extract(file_id, "expdata.+")
  quant_method_label <- str_match(file_id, "\\/([^/]+)\\/[^/]+$")[, 2]
  df_quant_wide_i <- fread(file_id)
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  df_quant_pep2pro <- df_quant_wide_i %>%
    dplyr::select(peptide, protein, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1:2, variable.name = "run_id") %>%
    mutate(intensity = exp(value))
  
  method_j <- quantification_methods[match(quant_method_label,label_names)]
  
  df_final <- quantify_by_x(df_quant = df_quant_pep2pro,
                            group.header = "protein",
                            subgroup.header = "peptide",
                            sample.id.header = "run_id",
                            intensity.header = "intensity",
                            method = method_j,
                            mc.cores = 24)
  
  df_final_complete <- df_final %>%
    reshape2::dcast(., protein ~ run_id, value.var = "intensity")
  
  fwrite(df_final_complete,
         file = paste(file_path, "protein/", quant_method_label, "/", file_name, sep = ""))
  
})


## set negative-control features by two-way anova. 约40分钟 ---------
rm(list = ls())
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("expdata_log", all_files)]
queried_files <- queried_files[!grepl("simulated", queried_files)]

anova_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "processing ", file_id, "\n"))
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded)")
  file_name <- str_extract(file_id, "^.*(?=expdata)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  df_meta <- df_meta %>% mutate(batch = paste(mode, lab, sep = "_"))
  df_quant_wide_i <- fread(file_id)
  
  df_features <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_quant_long_i <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    select(feature, starts_with("D")) %>%
    mutate(across(starts_with("D"), ~ ifelse(is.na(.)|. < 0, 0, exp(.)))) %>%
    # filter(apply(., 1, function(x) sum(x == 0) < .2 * nrow(df_meta))) %>%
    filter(apply(., 1, function(x) sum(x == 0) == 0)) %>%
    reshape2::melt(., id = 1, variable.name = "run_id", value.name = "intensity") %>%
    left_join(., df_meta, by = "run_id")
  
  rm(df_quant_wide_i, df_meta)
  gc()
  all_features <- unique(df_quant_long_i$feature)
  
  anova_tables_i <- mclapply(all_features, function(feature_id) {
    
    df_anova_j <- df_quant_long_i %>%
      filter(feature == feature_id) %>%
      rstatix::anova_test(intensity ~ sample * batch)
    
    df_anova_j$feature <- feature_id
    
    return(df_anova_j)
    
  }, mc.cores = 16)
  
  df_anova_i <- rbindlist(anova_tables_i)
  rm(anova_tables_i)
  gc()
  
  df_anova_i_final <- df_anova_i %>%
    left_join(., df_features, by = "feature") %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
  
  file_id2save <- paste(file_name, "anova_twoway_test.csv", sep = "")
  fwrite(df_anova_i_final, file_id2save)
  
  rm(df_anova_i, df_anova_i_final, df_features, file_path, file_name)
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
})

