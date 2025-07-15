# ' Title: Batch effect diagnosis by PCA.
# ' Author: Qiaochu Chen
# ‘ Date: Jun 30th, 2025

library(parallel)
library(pbapply)
library(readxl)
library(dplyr)
library(tibble)
library(reshape2)
library(data.table)
library(stringr)
if (!"gPCA" %in% installed.packages()) {
  install.packages("./utils/gPCA_1.0.tar.gz", repos = NULL)
  library(gPCA)
} else {
  library(gPCA)
}


## Global ---------------------
if (!interactive()) pboptions(type = "none")
options(mc.cores = 4)
dir2work <- "./results"
if (!dir.exists(dir2work)) dir.create(dir2work)
dir2save1 <- paste(dir2work, "/figures", sep = "")
if (!dir.exists(dir2save1)) dir.create(dir2save1)
dir2save2 <- paste(dir2work, "/tables", sep = "")
if (!dir.exists(dir2save2)) dir.create(dir2save2)
set.seed(20250630)
source("./utils/PCA.R")


## input files------------
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("expdata_intensity", queried_files)]
queried_files <- queried_files[!grepl("1_expdata_gamma_log", queried_files)]


## calculate PVCA -----------------------
pvca_tables <- pblapply(queried_files, function(file_id) {
  
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded|simulated)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  if (dataset %in% c("quartet")) {
    df_meta <- df_meta %>%
      mutate(batch = paste(mode, lab, sep = "_")) %>%
      dplyr::select(run_id, sample, any_of(c("lab", "mode")))
  } else {
    df_meta <- df_meta %>%
      dplyr::select(run_id, batch, sample)
  }
  
  metadata <- df_meta %>%
    column_to_rownames("run_id")
  
  expr <- fread(file_id)
  
  exprdata_t <- expr %>%
    select(all_of(rownames(metadata))) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    na.omit %>%
    t
  
  pvca_results <- main_pvca(exprdata_t = exprdata_t, metadata = metadata, plot = FALSE)
  
  sub_pvca <- pvca_results$table %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(sub_pvca)
})
sub_pvca_final <- rbindlist(pvca_tables)
fwrite(sub_pvca_final, "./results/tables/pvca.csv")


## calculate gPCA -----------------------
gpca_tables <- pblapply(queried_files, function(file_id) {
  
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded|simulated)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  if (dataset %in% c("quartet")) {
    df_meta <- df_meta %>% mutate(batch = paste(mode, lab, sep = "_"))
  }
  
  metadata <- df_meta %>%
    column_to_rownames("run_id")
  
  expr <- fread(file_id)
  
  exprdata_t <- expr %>%
    select(all_of(rownames(metadata))) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    na.omit %>%
    t
  
  gpca_results <- gPCA.batchdetect(exprdata_t, metadata$batch, center=TRUE,
                                   nperm = 250, seed = 2025)
  
  sub_gpca <- data.frame(gpca_delta = gpca_results$delta,
                         gpca_p = gpca_results$p.val) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(sub_gpca)
})
sub_gpca_final <- rbindlist(gpca_tables)
fwrite(sub_gpca_final, "./results/tables/gpca.csv")


