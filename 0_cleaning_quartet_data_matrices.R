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


## create unexisted directory -------------------------
dir2work <- "./data/expfiles"
if (!dir.exists(dir2work)) dir.create(dir2work)

for (d in c("quartet", "simulated")) {
  dir_tmp1 <- paste(dir2work, "/", d, sep = "")
  if (!dir.exists(dir_tmp1)) dir.create(dir_tmp1)
  
  for (s in c("balanced", "confounded")) {
    dir_tmp <- paste(dir_tmp1, "/", s, sep = "")
    if (!dir.exists(dir_tmp)) dir.create(dir_tmp)
    
    dir2save_pre <- paste(dir_tmp, "/", "precursor", sep = "")
    if (!dir.exists(dir2save_pre)) dir.create(dir2save_pre)
    
    for (l in c("peptide", "protein")) {
      dir2save_p <- paste(dir_tmp, "/", l, sep = "")
      if (!dir.exists(dir2save_p)) dir.create(dir2save_p)
      
      for (m in c("maxlfq", "toppep3", "ibaq")) {
        dir2save_tmp <- paste(dir2save_p, "/", m, sep = "")
        if (!dir.exists(dir2save_tmp)) dir.create(dir2save_tmp)
      }
    }
  }
}


## Balanced design--------------------------
## get file paths
all_files <- list.files("./data/rawfiles/MaxQuant", full.names = T)
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
  mutate(batch = paste(mode, lab, sep = "_")) %>%
  mutate(run_id = paste(mode, lab, sample, tube, sep = "_")) %>%
  select(run_id, mode, lab, batch, sample, tube, order, raw_file)

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
fwrite(df_meta, "./data/expfiles/quartet/balanced/meta.csv")
fwrite(df_wide_final, "./data/expfiles/quartet/balanced/precursor/expdata_intensity.csv")
fwrite(df_log_final, "./data/expfiles/quartet/balanced/precursor/expdata_log.csv")


## Confounded design--------------------------
## set confounded design (2 sample types within each group)
df_meta_c <- df_meta %>%
  filter(grepl("D6|(DDA_APT_D5)|(DDA_FDU_F7)|(DDA_NVG_M8)|(DIA_APT_D5)|(DIA_FDU_F7)|(DIA_BGI_M8)", run_id))

## transform list into matrix
df_wide_c <- df_wide_final %>%
  select(precursor, peptide, protein, mz, rt, all_of(df_meta_c$run_id)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 36))

## transform intensity into log-scale
df_log_c <- df_log_final %>%
  select(precursor, peptide, protein, mz, rt, all_of(df_meta_c$run_id)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 36))

## save
fwrite(df_meta_c, "./data/expfiles/quartet/confounded/meta.csv")
fwrite(df_wide_c, "./data/expfiles/quartet/confounded/precursor/expdata_intensity.csv")
fwrite(df_log_c, "./data/expfiles/quartet/confounded/precursor/expdata_log.csv")


