# ' Title: Signal-to-noise ratio based on PCA.
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


## global ---------------------
if (!interactive()) pboptions(type = "none")
options(mc.cores = 4)
dir2work <- "./results"
if (!dir.exists(dir2work)) dir.create(dir2work)
dir2save1 <- paste(dir2work, "/figures", sep = "")
if (!dir.exists(dir2save1)) dir.create(dir2save1)
dir2save2 <- paste(dir2work, "/tables", sep = "")
if (!dir.exists(dir2save2)) dir.create(dir2save2)
source("./utils/PCA.R")
set.seed(20250630)


## input files------------
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("expdata_intensity", queried_files)]
queried_files <- queried_files[!grepl("1_expdata_gamma_log", queried_files)]


## calculate PCA and SNR --------------------
pca_objs <- pblapply(queried_files, function(file_id) {
  
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded|simulated)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  metadata <- df_meta %>%
    dplyr::select(run_id, sample, everything()) %>%
    column_to_rownames("run_id")
  
  expr <- fread(file_id)
  
  exprdata_t <- expr %>%
    dplyr::select(all_of(rownames(metadata))) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    na.omit %>%
    t
  
  pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                          center = TRUE, scale = TRUE, group = "sample",
                          biplot = FALSE, dictGroups = metadata$sample,
                          snr = TRUE, plot = FALSE)
  
  pca_objs_i <- c(pca_results,
                   list(dataset = dataset,
                        scenario = scenario,
                        data_level = data_level,
                        correct_level = correct_level,
                        correct_method = correct_method,
                        quant_method = quant_method))
  
  return(pca_objs_i)
})
names(pca_objs) <- queried_files
saveRDS(pca_objs, paste("./results/tables/pca_snr.rds", sep = ""))

