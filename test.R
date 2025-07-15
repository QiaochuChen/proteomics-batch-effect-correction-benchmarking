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
library(RUVIIIC)


rm(list = ls())
if (!interactive()) {
  pboptions(type = "none")
  print("Non-interactive session detected. Disabling progress bars.")
}
if (!dir.exists("./data/")) dir.create("./data/")
if (!dir.exists("./data/expfiles")) dir.create("./data/expfiles")
setwd("./data/")


## Balanced design--------------------------
## get file paths
all_files <- list.files("./rawfiles/MaxQuant", full.names = T)
all_dirs <-  all_files[!grepl("\\/raw$", all_files)]

## get precursor/peptide/protein-level data table
combined_tables <- pblapply(all_dirs, function(dir) {
  
  evidence_file_path <- paste(dir, "/evidence.txt", sep = "")
  evidence_txt <- fread(evidence_file_path, showProgress = FALSE)
  
  df_evidence <- evidence_txt %>%
    rename_with(., ~ tolower(gsub(" ", "_", .x, fixed = TRUE))) %>%
    filter(intensity > 0, !grepl(";", protein_names), !protein_names %in% "")
  
  df_combined <- df_evidence %>%
    mutate(unique_peptide = paste(modified_sequence, charge, sep = "+")) %>%
    group_by(raw_file, unique_peptide) %>%
    mutate(label = 1:length(unique_peptide)) %>%
    mutate(precursor =  paste(unique_peptide, label, sep = "_")) %>%
    ungroup() %>%
    dplyr::rename(peptide = sequence, protein = protein_names) %>%
    distinct(raw_file, precursor, peptide, protein, `m/z`, retention_time, intensity) %>%
    mutate(lab = str_extract(dir, "(?<=MaxQuant\\/).+(?=_)"),
           mode = str_extract(dir, "(?<=_).+"))
  
  return(df_combined)
  
})
names(combined_tables) <- sapply(all_dirs, function(x) str_extract(x, "(?<=MaxQuant\\/).+"))

## extract metadata
df_meta <- combined_tables %>%
  rbindlist %>%
  distinct(raw_file, lab, mode) %>%
  mutate(name_tmp = ifelse(grepl("20200706-[1-9]$", raw_file),
                           gsub("20200706-", "20200706-0", raw_file),
                           raw_file)) %>%
  arrange(name_tmp, raw_file, mode, lab) %>%
  mutate(order = c(rep(1:12, 5), c(1,5,9,2,6,10,3,7,11,4,8,12))) %>%
  arrange(mode, lab, order) %>%
  mutate(sample = rep(c("D5", "D6", "F7", "M8"), 18)) %>%
  mutate(tube = rep(rep(1:3, each = 4), 6)) %>%
  mutate(run_id = paste(mode, lab, sample, tube, sep = "_")) %>%
  select(run_id, mode, lab, sample, tube, order, raw_file)

## transform list into matrix
df_wide <- combined_tables %>%
  rbindlist %>%
  filter(intensity != 0) %>%
  group_by(precursor, peptide, protein, `m/z`) %>%
  mutate(mz = `m/z`, rt = mean(retention_time)) %>%
  ungroup() %>%
  reshape2::dcast(., precursor + peptide + protein + mz + rt ~ raw_file,
                  value.var = "intensity")

dictColumnNames <- df_meta$run_id
names(dictColumnNames) <- df_meta$raw_file

df_wide_final <- df_wide %>% plyr::rename(dictColumnNames)

## transform intensity into log-scale
df_log_final <- df_wide_final %>%
  mutate_if(is.numeric, ~ ifelse(. == 0, NA, log(.)))

## save
fwrite(df_meta, "./expfiles/meta_balanced.csv")

