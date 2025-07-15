# ' Title: ChiHOPE
# ' Author: Qiaochu Chen
# ‘ Date: Jul 10th, 2025

library(parallel)
library(pbapply)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(diann)


if (!interactive()) pboptions(type = "none")
options(mc.cores = min(24, max(1, parallel::detectCores() - 8)))
set.seed(20250630)


## create unexisted directory -------------------------
dir2work <- "./ChiHOPE"
if (!dir.exists(dir2work)) dir.create(dir2work)

dir2save_pre <- paste(dir2work, "/", "precursor", sep = "")
if (!dir.exists(dir2save_pre)) dir.create(dir2save_pre)

for (l in c("peptide", "protein")) {
  dir2save_p <- paste(dir2work, "/", l, sep = "")
  if (!dir.exists(dir2save_p)) dir.create(dir2save_p)
  
  for (m in c("maxlfq")) {
    dir2save_tmp <- paste(dir2save_p, "/", m, sep = "")
    if (!dir.exists(dir2save_tmp)) dir.create(dir2save_tmp)
  }
}
for (l in c("peptide", "protein")) {
  for (m in c("maxlfq")) {
    dir_tmp <- paste(dir2work, "/", l, "/", m, "/precursor_corrected/", sep = "")
    if (!dir.exists(dir_tmp)) dir.create(dir_tmp)
  }
}
for (m in c("maxlfq")) {
  dir_tmp <- paste(dir2work, "/protein/", m, "/peptide_corrected/", sep = "")
  if (!dir.exists(dir_tmp)) dir.create(dir_tmp)
}


## cleaning meta data into uniformed format -----------------
df_meta <- fread("./ChiHOPE/metadata.csv")
df_meta <- df_meta %>%
  rename(run_id = ID, batch = Cleaning, order = Order, sample = Type) %>%
  mutate_at("batch", ~ paste("Batch", as.numeric(factor(.)), sep = "")) %>%
  select(run_id, batch, order, sample, everything())

fwrite(df_meta, "./ChiHOPE/meta.csv")


## aggregate precursors into peptides: about  min -----------------
rm(list = ls())
gc()
source("./utils/BEC.R")
file_id <- "./ChiHOPE/precursor/expdata_log.csv"

cat(paste0("\n[", Sys.time(), "] ", "reading ", file_id, "\n"))
file_path <- str_extract(file_id, ".+(?=precursor)")
file_name <- str_extract(file_id, "expdata.+")
# df_quant_wide_i <- fread(file_id)
df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))

cat(paste0("\n[", Sys.time(), "] ", "missing values cutoff: 20%", "\n"))
# df_test <- df_quant_wide_i %>%
#   rowwise() %>%
#   mutate(na_sum = sum(is.na(c_across(where(is.numeric))))) %>%
#   ungroup() %>%
#   filter(na_sum <= .2 * nrow(df_meta)) %>%
#   dplyr::select(precursor, peptide, protein, all_of(df_meta$run_id))

# fwrite(df_test, "./ChiHOPE/precursor/expdata_log_cutoff_NAs.csv")
df_test <- fread("./ChiHOPE/precursor/expdata_log_cutoff_NAs.csv")

# rm(df_quant_wide_i)
gc()

# test_proteins <- unique(df_test$protein)[1:100] ## test only
# df_test <- df_test %>% filter(protein %in% test_proteins)

cat(paste0("\n[", Sys.time(), "] ", "aggregating ",
    nrow(df_test), " precursors into ",
    length(unique(df_test$peptide)), " peptides and ",
    length(unique(df_test$protein)), " proteins.", "\n"))

df_quant_pre2pep <- df_test %>%
  dplyr::select(precursor, peptide, protein, all_of(df_meta$run_id)) %>%
  reshape2::melt(., id = 1:3, variable.name = "run_id", na.rm = TRUE) %>%
  mutate(peptide2protein = paste(peptide, protein, sep = "_")) %>%
  mutate(intensity = exp(value))

df_final <- quantify_by_x(df_quant = df_quant_pre2pep,
                          group.header = "peptide2protein",
                          subgroup.header = "precursor",
                          sample.id.header = "run_id",
                          intensity.header = "intensity",
                          method = "MaxLFQ",
                          mc.cores = 2)

rm(df_quant_pre2pep)
gc()

df_final_complete <- df_final %>%
  reshape2::dcast(., peptide2protein ~ run_id, value.var = "intensity") %>%
  tidyr::separate(peptide2protein, c("peptide", "protein"), sep = "_")

rm(df_final)
gc()

cat(paste0("\n[", Sys.time(), "] ", nrow(df_final_complete), "rows aggregated", "\n"))
fwrite(df_final_complete,
       file = paste(file_path, "peptide/maxlfq/", file_name, sep = ""))


## aggregate peptides into proteins: about  min -----------------
rm(list = ls())
gc()
source("./utils/BEC.R")
file_id <- "./ChiHOPE/peptide/maxlfq/expdata_log.csv"

cat(paste0("\n[", Sys.time(), "] ", "reading ", file_id, "\n"))
file_path <- str_extract(file_id, ".+(?=peptide)")
file_name <- str_extract(file_id, "expdata.+")
df_quant_wide_i <- fread(file_id)
df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))

cat(paste0("\n[", Sys.time(), "] ", "processing ", nrow(df_quant_wide_i), "rows", "\n"))
df_quant_pep2pro <- df_quant_wide_i %>%
  dplyr::select(peptide, protein, all_of(df_meta$run_id)) %>%
  reshape2::melt(., id = 1:2, variable.name = "run_id") %>%
  mutate(intensity = exp(value))

rm(df_quant_wide_i)
gc()

df_final <- quantify_by_x(df_quant = df_quant_pep2pro,
                          group.header = "protein",
                          subgroup.header = "peptide",
                          sample.id.header = "run_id",
                          intensity.header = "intensity",
                          method = "MaxLFQ",
                          mc.cores = 2)

rm(df_quant_pep2pro)
gc()

df_final_complete <- df_final %>%
  reshape2::dcast(., protein ~ run_id, value.var = "intensity") %>%
  mutate_all(~ ifelse(. < 0, NA, .))

rm(df_final)
gc()

cat(paste0("\n[", Sys.time(), "] ", nrow(df_final_complete), "rows aggregated", "\n"))
fwrite(df_final_complete,
       file = paste(file_path, "protein/maxlfq/", file_name, sep = ""))


## loess corrected at precursor/peptide/protein levels. about  min -------
rm(list = ls())
gc()
source("./utils/BEC.R")
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("expdata_log", all_files)]
queried_files <- queried_files[!grepl("precursor/expdata_log.csv", queried_files)]

corrected_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "correcting ", file_id, "\n"))
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
  file_name <- str_extract(file_id, "^.*(?=expdata)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  df_quant_wide_i <- fread(file_id)
  
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  
  df_features <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_quant_wide_i <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    dplyr::select(feature, all_of(df_meta$run_id)) %>%
    mutate(across(all_of(df_meta$run_id), ~ ifelse(is.na(.)|. < 0, 0, exp(.))))
  
  ## Quantile + LOESS normalized
  df_loess <- loess_solve(df_quant_wide_i, df_meta)
  
  df_corrected_j_final <- df_loess %>%
    mutate_at("feature", as.integer) %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
  
  file_id2save <- paste(file_name, "expdata_loess.csv", sep = "")
  fwrite(df_corrected_j_final, file = file_id2save)
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
  
  cat(paste0("[", Sys.time(), "] ", "Successfully corrected.\n"))
})


## set negative-control features by two-way anova. about  min ---------
rm(list = ls())
gc()
source("./utils/BEC.R")
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("expdata_loess", all_files)]

anova_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "processing ", file_id, "\n"))
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
  file_name <- str_extract(file_id, "^.*(?=expdata)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  df_quant_wide_i <- fread(file_id)
  
  df_features <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_quant_long_i <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    select(feature, starts_with("Exp")) %>%
    mutate(across(starts_with("Exp"), ~ ifelse(is.na(.)|. < 0, 0, exp(.)))) %>%
    # filter(apply(., 1, function(x) sum(x == 0) < .2 * nrow(df_meta))) %>%
    filter(apply(., 1, function(x) sum(x == 0) == 0)) %>%
    reshape2::melt(., id = 1, variable.name = "run_id", value.name = "intensity") %>%
    left_join(., df_meta, by = "run_id") %>%
    filter(!sample %in% "Study sample")
  
  rm(df_quant_wide_i, df_meta)
  gc()
  all_features <- unique(df_quant_long_i$feature)
  
  anova_tables_i <- mclapply(all_features, function(feature_id) {
    
    df_anova_j <- df_quant_long_i %>%
      filter(feature == feature_id) %>%
      rstatix::anova_test(intensity ~ sample * batch)
    
    df_anova_j$feature <- feature_id
    
    return(df_anova_j)
    
  }, mc.cores = 8)
  
  df_anova_i <- rbindlist(anova_tables_i)
  rm(anova_tables_i)
  gc()
  
  df_anova_i_final <- df_anova_i %>%
    left_join(., df_features, by = "feature") %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
  
  file_id2save <- paste(file_name, "anova_twoway_test_loess.csv", sep = "")
  fwrite(df_anova_i_final, file_id2save)
  
  rm(df_anova_i, df_anova_i_final, df_features, file_path, file_name)
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
})


## correct batch effects at precursor/peptide/protein levels. about  min -------
rm(list = ls())
gc()
source("./utils/BEC.R")
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
# queried_files <- all_files[grepl("expdata_loess", all_files)]
queried_files <- all_files[grepl("expdata_log", all_files)]
queried_files <- queried_files[grepl("protein", queried_files)]

correct_methods <- c("ratio_by_d6", "median_centering", "combat", "harmony", "WaveICA", "RUV-III-C")
label_names <- c("ratio", "median", "combat", "harmony", "waveica", "ruv")

corrected_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "correcting ", file_id, "\n"))
  file_path <- str_extract(file_id, "^.*\\/(ChiHOPE)")
  file_name <- str_extract(file_id, "^.*(?=expdata)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  df_quant_wide_i <- fread(file_id)
  
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  
  # df_anova_i <- fread(paste(file_name, "anova_twoway_test_loess.csv", sep = ""))
  df_anova_i <- fread(paste(file_name, "anova_twoway_test.csv", sep = ""))
  
  ctrl_features <- df_anova_i %>%
    reshape2::dcast(., feature ~ Effect, value.var = "p<.05") %>%
    filter(`sample` %in% "" & `sample:batch` %in% "") %>%
    pull(feature)
  
  rm(df_anova_i)
  gc()
  
  df_meta2 <- df_meta %>%
    dplyr::rename(run_id_old = run_id) %>%
    mutate(run_id = paste(run_id_old, batch, sample, sep = "_")) %>%
    select(run_id, everything())
  dictColNames <- setNames(df_meta2$run_id_old, df_meta2$run_id)
  
  df_features <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_quant_wide_i <- df_quant_wide_i %>%
    tibble::rowid_to_column("feature") %>%
    dplyr::select(feature, all_of(df_meta2$run_id_old)) %>%
    dplyr::rename(dictColNames) %>%
    mutate(across(all_of(df_meta2$run_id), ~ ifelse(is.na(.)|. < 0, 0, exp(.))))
  
  corrected_tables_i <- mclapply(1:4, function(j) {
    
    if (!correct_methods[j] %in% "RUV-III-C") ctrl_features <- NULL
    
    null_device <- if (.Platform$OS.type == "unix") "/dev/null" else "NUL"
    con <- file(null_device, open = "w")
    sink(con, type = "output")
    sink(con, type = "message")
    df_corrected_j <- correct_by_x(df_quant_wide_i, df_meta2, correct_methods[j],
                                   ctrl_features, reference_sample = "PM")
    sink(type = "message")
    sink(type = "output")
    close(con)
    
    rm(df_quant_wide_i, df_meta)
    gc()
    
    dictColNames2 <- setNames(df_meta2$run_id, df_meta2$run_id_old)
    
    df_corrected_j_final <- df_corrected_j %>%
      mutate_at("feature", as.integer) %>%
      left_join(., df_features, by = "feature") %>%
      dplyr::select(!feature) %>%
      dplyr::rename(dictColNames2) %>%
      dplyr::select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
    
    # file_id2save <- paste(file_name, "expdata_", label_names[j], "_loess.csv", sep = "")
    file_id2save <- paste(file_name, "expdata_", label_names[j], ".csv", sep = "")
    fwrite(df_corrected_j_final, file = file_id2save)
    
    gc()
    cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
    
  }, mc.cores = 5)
  
  ## RUV-III-C单独运行
  df_corrected_j <- correct_by_x(df_quant_wide_i, df_meta, "RUV-III-C", ctrl_features)
  rm(df_quant_wide_i, df_meta)
  gc()
  
  dictColNames2 <- setNames(df_meta2$run_id, df_meta2$run_id_old)
  
  df_corrected_j_final <- df_corrected_j %>%
    mutate_at("feature", as.integer) %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    dplyr::rename(dictColNames2) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything())
  
  # file_id2save <- paste(file_name, "expdata_ruv_loess.csv", sep = "")
  file_id2save <- paste(file_name, "expdata_ruv.csv", sep = "")
  fwrite(df_corrected_j_final, file = file_id2save)
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
  
  cat(paste0("[", Sys.time(), "] ", "Successfully corrected.\n"))
})


## aggregate precursors into peptides: about  min -----------------
rm(list = ls())
gc()
source("./utils/BEC.R")
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("precursor/expdata", all_files)]
queried_files <- queried_files[!grepl("expdata_log", queried_files)]
pep_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "quantifying ", file_id, "\n"))
  file_path <- str_extract(file_id, ".+(?=precursor)")
  file_name <- str_extract(file_id, "expdata.+")
  df_quant_wide_i <- fread(file_id)
  df_meta <- fread(paste(file_path, "meta.csv", sep = ""))
  
  df_quant_pre2pep <- df_quant_wide_i %>%
    select(precursor, peptide, protein, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1:3, variable.name = "run_id", na.rm = TRUE) %>%
    mutate(peptide2protein = paste(peptide, protein, sep = "_")) %>%
    mutate(intensity = exp(value))
  
  df_final <- quantify_by_x(df_quant = df_quant_pre2pep,
                            group.header = "peptide2protein",
                            subgroup.header = "precursor",
                            sample.id.header = "run_id",
                            intensity.header = "intensity",
                            method = "MaxLFQ",
                            mc.cores = 2)
  
  df_final_complete <- df_final %>%
    reshape2::dcast(., peptide2protein ~ run_id, value.var = "intensity") %>%
    tidyr::separate(peptide2protein, c("peptide", "protein"), sep = "_")
  
  file_id2save <- paste(file_path, "peptide/maxlfq/precursor_corrected/", file_name, sep = "")
  fwrite(df_final_complete, file = file_id2save)
  
  gc()
  cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))
  
  cat(paste0("[", Sys.time(), "] ", "Successfully quantified into peptides.\n"))
})


## aggregate peptides into proteins: about  min -----------------
rm(list = ls())
gc()
source("./utils/BEC.R")
quantification_methods <- c("MaxLFQ", "iBAQ", "TopPep3")
label_names <- c("maxlfq", "ibaq", "toppep3")

all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("peptide/", all_files)]
queried_files <- queried_files[grepl("expdata", queried_files)]
queried_files <- queried_files[!grepl("expdata_log", queried_files)]
pro_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "quantifying ", file_id, "\n"))
  file_path <- str_extract(file_id, ".+(?=peptide)")
  file_name <- str_extract(file_id, "expdata.+")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  if (grepl("precursor_corrected", file_id)) {
    quant_method_label <- str_match(file_id, "\\/([^/]+)\\/([^/]+)\\/[^/]+$")[, 2]
    correct_level_label <- "precursor_corrected"
  } else {
    quant_method_label <- str_match(file_id, "\\/([^/]+)\\/[^/]+$")[, 2]
    correct_level_label <- "peptide_corrected"
  }
  
  df_quant_wide_i <- fread(file_id)
  
  df_quant_pep2pro <- df_quant_wide_i %>%
    select(peptide, protein, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1:2, variable.name = "run_id") %>%
    mutate(intensity = exp(value))
  
  method_j <- quantification_methods[match(quant_method_label,label_names)]
  
  df_final <- quantify_by_x(df_quant = df_quant_pep2pro,
                            group.header = "protein",
                            subgroup.header = "peptide",
                            sample.id.header = "run_id",
                            intensity.header = "intensity",
                            method = method_j,
                            mc.cores = 2)
  
  df_final_complete <- df_final %>%
    reshape2::dcast(., protein ~ run_id, value.var = "intensity")
  
  file_id2save <- paste(file_path, "protein/", quant_method_label, "/",
                        correct_level_label, "/", file_name, sep = "")
  fwrite(df_final_complete, file = file_id2save)
  
  cat(paste0("[", Sys.time(), "] ", "Successfully quantified into proteins.\n"))
})

df_test <- fread("./ChiHOPE/peptide/maxlfq/expdata_log.csv")
nrow(df_test)

df_test <- fread("./ChiHOPE/protein/maxlfq/expdata_log.csv")
nrow(df_test)
