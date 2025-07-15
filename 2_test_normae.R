# ' Title: Prepare for NormAE
# ' Author: Qiaochu Chen
# ‘ Date: Jun 26th, 2025

library(parallel)
library(pbapply)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)


if (!interactive()) pboptions(type = "none")

df_data <- fread("./data/expfiles/balanced/precursor/expdata_intensity.csv")
df_meta <- fread("./data/expfiles/balanced/meta.csv")

df_data_test <- df_data %>%
  mutate(feature = paste(precursor, peptide, protein, sep = "=")) %>%
  column_to_rownames("feature") %>%
  select(mz, rt, starts_with("D")) %>%
  mutate_at(3:74, ~ ifelse(is.na(.), 0, .))

df_meta_test <- df_meta %>%
  column_to_rownames("run_id") %>%
  mutate(injection.order = 1:72) %>%
  mutate(batch = paste(mode, lab, sep = "_")) %>%
  mutate(class = ifelse(sample %in% "D6", "QC", "Subject")) %>%
  select(class, injection.order, batch)

write.csv(df_data_test, "./data/expfiles/balanced/precursor/test_normae_data.csv")
write.csv(df_meta_test, "./data/expfiles/balanced/precursor/test_normae_info.csv")

## after NormAE correction
df_data <- read.csv("./data/expfiles/balanced/precursor/X_clean.csv", row.names = 1)

df_data_final <- df_data %>%
  rownames_to_column("feature") %>%
  tidyr::separate(feature, c("precursor", "peptide", "protein"), sep = "=")

fwrite(df_data_final, "./data/expfiles/balanced/precursor/expdata_normae.csv")

