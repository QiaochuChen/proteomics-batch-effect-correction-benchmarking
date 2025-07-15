# ' Title: Coefficient of variation.
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
set.seed(20250630)


## input files------------
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("expdata_intensity", queried_files)]


## calculate CV ------------------------
cv_results <- pblapply(queried_files, function(file_id) {
  
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded|simulated)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  expr <- fread(file_id)
  
  df_features <- expr %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_expr <- expr %>%
    tibble::rowid_to_column("feature") %>%
    select(feature, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1, variable.name = "run_id", na.rm = TRUE) %>%
    left_join(., df_meta, by = "run_id")
  
  df_cv <- df_expr %>%
    mutate_at("value", exp) %>%
    group_by(feature, sample) %>%
    summarise(cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
              .groups = "drop") %>%
    ungroup() %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything()) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(df_cv)
})
names(cv_results) <- queried_files
saveRDS(cv_results, paste("./results/tables/cv.rds", sep = ""))


