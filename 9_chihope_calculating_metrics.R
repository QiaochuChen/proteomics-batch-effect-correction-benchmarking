# ' Title: ChiHOPE.
# ' Author: Qiaochu Chen
# ‘ Date: Jul 10th, 2025

library(parallel)
library(pbapply)
library(readxl)
library(openxlsx)
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
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("precursor/expdata_log.csv", queried_files)]
queried_files <- queried_files[grepl("protein", queried_files)]


## calculate CV ------------------------
cv_results <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
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
    filter(!sample %in% "Study sample") %>%
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
saveRDS(cv_results, paste("./results/tables/cv_chihope.rds", sep = ""))


## calculate PCA and SNR --------------------
source("./utils/PCA.R")
pca_objs <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  metadata <- df_meta %>%
    dplyr::select(run_id, sample, everything()) %>%
    filter(!sample %in% "Study sample") %>%
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
                          pct_threshold = NULL,
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
saveRDS(pca_objs, paste("./results/tables/pca_snr_chihope.rds", sep = ""))


## calculate PCA --------------------
source("./utils/PCA.R")
pca_objs <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
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
                          snr = FALSE, plot = FALSE)
  
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
saveRDS(pca_objs, paste("./results/tables/pca_chihope.rds", sep = ""))


## calculate UMAP --------------------
source("./utils/PCA.R")
umap_objs <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
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
  
  umap_results <- calculate_umap(exprdata_t = exprdata_t, metadata = metadata,
                                center = TRUE, scale = TRUE, group = "sample",
                                pct_threshold = .8)
  
  umap_objs_i <- umap_results %>%
    dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(umap_objs_i)
})
names(umap_objs) <- queried_files
saveRDS(umap_objs, paste("./results/tables/umap_chihope.rds", sep = ""))


## calculate PVCA -----------------------
pvca_tables <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  metadata <- df_meta %>%
    dplyr::select(run_id, sample, batch, Date) %>%
    filter(!sample %in% "Study sample") %>%
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
fwrite(sub_pvca_final, "./results/tables/pvca_chihope.csv")


## prepare for model -----------------------------
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("precursor/expdata_log.csv", queried_files)]
queried_files <- queried_files[grepl("protein", queried_files)]
# queried_files <- queried_files[grepl("log|loess", queried_files)]

halfwide_list <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  df_quant_wide_i <- fread(file_id)
  
  df_quant_long_i <- df_quant_wide_i %>%
    reshape2::melt(., id = 1, variable.name = "run_id", na.rm = TRUE) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(df_quant_long_i)
})

df_meta <- fread("./ChiHOPE/meta.csv")

df_quant_wide_t <- halfwide_list %>%
  rbindlist %>%
  reshape2::dcast(., correct_method + run_id ~ protein,
                  value.var = "value") %>%
  left_join(., df_meta, by = "run_id") %>%
  # filter(Week %in% "Baseline") %>%
  filter(!grepl("QC", sample))

unique(df_quant_wide_t$correct_method)
write.xlsx(df_quant_wide_t, "./results/tables/expdata_chihope_1431.xlsx")
# write.xlsx(df_quant_wide_t, "./results/tables/expdata_chihope_baseline_744.xlsx")  

## randomize true labels --------------------
set.seed(20250630)
df_quant_wide_t_tmp <- df_quant_wide_t %>%
  filter(correct_method %in% "log") %>%
  mutate_at("Sex", ~ sample(.)) %>%
  mutate_at("Age", ~ sample(.))

write.xlsx(df_quant_wide_t_tmp, "./results/tables/expdata_chihope_1431_negctrl.xlsx")

df_tmp <- read_excel("./results/tables/expdata_chihope_1431_negctrl.xlsx")
length(unique(df_tmp$Subject))
length(unique(df_tmp$Subject[df_tmp$Training.Validation == 1]))
length(unique(df_tmp$Subject[df_tmp$Training.Validation == 2]))


## calculate MCC & R square -----------------------------
df_predicted <- read_excel("./results/tables/expdata_chihope_1431_with_predictions.xlsx")
df_predicted_test <- df_predicted %>%
  filter(Training.Validation == 2) %>%
  select(correct_method, Training.Validation, Sex, Age, prediction_sex, prediction_age) %>%
  mutate(label = case_when(Sex %in% "F" & prediction_sex %in% "F" ~ "TP",
                           Sex %in% "F" & prediction_sex %in% "M" ~ "FN",
                           Sex %in% "M" & prediction_sex %in% "M" ~ "TN",
                           Sex %in% "M" & prediction_sex %in% "F" ~ "FP"))

df_predicted <- read_excel("./results/tables/expdata_chihope_1431_negctrl_with_predictions.xlsx")
df_predicted_negctrl <- df_predicted %>%
  filter(Training.Validation == 2) %>%
  mutate(correct_method = "negctrl") %>%
  select(correct_method, Training.Validation, Sex, Age, prediction_sex, prediction_age) %>%
  mutate(label = case_when(Sex %in% "F" & prediction_sex %in% "F" ~ "TP",
                           Sex %in% "F" & prediction_sex %in% "M" ~ "FN",
                           Sex %in% "M" & prediction_sex %in% "M" ~ "TN",
                           Sex %in% "M" & prediction_sex %in% "F" ~ "FP"))

df_mcc <- df_predicted_test %>%
  rbind(., df_predicted_negctrl) %>%
  reshape2::dcast(., correct_method + Training.Validation ~ label, fun.aggregate = length) %>%
  mutate_at(3:6, as.numeric) %>%
  mutate(mcc = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))) %>%
  mutate(recall = TP / (TP + FN),
         precision = TP / (TP + FP),
         f1 = 2 * precision * recall / (precision + recall))

df_rsquare <- df_predicted_test %>%
  rbind(., df_predicted_negctrl) %>%
  group_by(correct_method) %>%
  summarise(ss_res = sum((Age - prediction_age)^2),
            ss_tot = sum((Age - mean(Age)) ^ 2)) %>%
  mutate(r_square = 1 - (ss_res / ss_tot))

df_final <- full_join(df_mcc, df_rsquare)

fwrite(df_final, "./results/tables/mcc_rsquare_chihope.csv")

