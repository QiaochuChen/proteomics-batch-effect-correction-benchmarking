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


## global. ------------------
if (!interactive()) pboptions(type = "none")
options(mc.cores = min(48, max(1, parallel::detectCores() - 8)))

for (d in c("quartet", "simulated")) {
  for (s in c("balanced", "confounded")) {
    for (l in c("peptide", "protein")) {
      for (m in c("maxlfq", "toppep3", "ibaq")) {
        dir_tmp <- paste("./data/expfiles/", d, "/", s, "/", l, "/", m, "/precursor_corrected/", sep = "")
        if (!dir.exists(dir_tmp)) dir.create(dir_tmp)
      }
    }
    for (m in c("maxlfq", "toppep3", "ibaq")) {
      dir_tmp <- paste("./data/expfiles/", d, "/", s, "/protein/", m, "/peptide_corrected/", sep = "")
      if (!dir.exists(dir_tmp)) dir.create(dir_tmp)
    }
  }
}
source("./utils/BEC.R")
set.seed(20250630)


## set methods ------------
quantification_methods <- c("MaxLFQ", "iBAQ", "TopPep3")
label_names <- c("maxlfq", "ibaq", "toppep3")


## aggregate precursors into peptides: about 90 min -----------------
gc()
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("precursor/expdata", all_files)]
queried_files <- queried_files[!grepl("expdata_log", queried_files)]
queried_files <- queried_files[!grepl("expdata_intensity", queried_files)]
pep_tables <- pblapply(queried_files, function(file_id) {
  
  cat(paste0("\n[", Sys.time(), "] ", "quantifying ", file_id, "\n"))
  file_path <- str_extract(file_id, ".+(?=precursor)")
  file_name <- str_extract(file_id, "expdata.+")
  df_quant_wide_i <- fread(file_id)
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  df_quant_pre2pep <- df_quant_wide_i %>%
    select(precursor, peptide, protein, all_of(df_meta$run_id)) %>%
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
    
    file_id2save <- paste(file_path, "peptide/", label_names[j], "/precursor_corrected/", file_name, sep = "")
    fwrite(df_final_complete, file = file_id2save)
    
    gc()
    cat(paste0("[", Sys.time(), "] ", "Successfully saved: ", file_id2save, "\n"))

  }, mc.cores = 3)

  cat(paste0("[", Sys.time(), "] ", "Successfully quantified into peptides.\n"))
})


## aggregate peptides into proteins: about 5 min -----------------
rm(pep_tables)
gc()
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
all_files <- list.files("./data/expfiles/simulated/confounded", recursive = TRUE, full.names = TRUE)
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
                            mc.cores = 24)
  
  df_final_complete <- df_final %>%
    reshape2::dcast(., protein ~ run_id, value.var = "intensity")
  
  file_id2save <- paste(file_path, "protein/", quant_method_label, "/",
                        correct_level_label, "/", file_name, sep = "")
  fwrite(df_final_complete, file = file_id2save)
  
  cat(paste0("[", Sys.time(), "] ", "Successfully quantified into proteins.\n"))
})

