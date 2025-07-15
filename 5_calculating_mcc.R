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
source("./utils/DEP.R")


## input files------------
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl(".+expdata", all_files)]
queried_files <- queried_files[!grepl("1_expdata_gamma_log", queried_files)]
queried_files <- queried_files[!grepl("expdata_intensity", queried_files)]
queried_files <- queried_files[grepl("simulated", queried_files)]


## DEP analysis ------------------
dep_tables <- pblapply(queried_files, function(file_id) {
  
  print(file_id)
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
    tibble::rowid_to_column("feature")
  
  sample_pairs <- list(c("group2", "group1"), c("group3", "group1"))
  
  dep_list <- mclapply(sample_pairs, function(sample_pair) {
    
    metadata <- df_meta %>%
      filter(sample %in% sample_pair)
    
    exprdata <- df_expr %>%
      column_to_rownames("feature") %>%
      select(all_of(metadata$run_id)) %>%
      mutate_all(~ ifelse(is.na(.)|. == Inf, 0, exp(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0))
    
    group <- factor(metadata$sample, levels = sample_pair, ordered = TRUE)
    
    df_dep_j <- calculate_p(exprdata, group, na_threshold = 0,
                            test_method = "t", p_adjust = "BH",
                            transform_to_log2 = TRUE)
    
    return(df_dep_j)
    
  }, mc.cores = 3)
  
  df_dep <- dep_list %>%
    rbindlist(.) %>%
    mutate_at("feature", as.integer) %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything()) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(df_dep)
  
})
names(dep_tables) <- queried_files
saveRDS(dep_tables, "./results/tables/dep.rds")


## calculate mcc score ------------------
dep_tables <- readRDS("./results/tables/dep.rds")
df_mapping <- fread("./data/expfiles/simulated/mapping.csv")
df_truth <- df_mapping %>% reshape2::melt(., id = 1:4, variable.name = "sample_pair")
mcc_objs <- pblapply(dep_tables, function(tmp_table) {
  
  df_combined_i <- tmp_table %>%
    # mutate_at("estimate", ~ ifelse(p > .05, 0, .)) %>%
    mutate_at("estimate", ~ ifelse(abs(.) < log(1.2), 0, .)) %>%
    mutate(sample_pair = paste(group1, group2, sep = "/")) %>%
    left_join(., df_truth) %>%
    select(protein, class, value, estimate, everything()) %>%
    mutate(label = apply(., 1, function(x) {
      if (as.numeric(x[3]) == 0 & as.numeric(x[4]) == 0) {
        a <- "TN"
      } else if (as.numeric(x[3]) == 0 & as.numeric(x[4]) != 0) {
        a <- "FP"
      } else if (as.numeric(x[3]) != 0 & as.numeric(x[4]) == 0) {
        a <- "FN"
      } else if (as.numeric(x[3]) * as.numeric(x[4]) > 0) {
        a <- "TP"
      } else if (as.numeric(x[3]) * as.numeric(x[4]) < 0) {
        a <- "FP"
      }
      return(a)
    })) 
  
  df_mcc_i <- df_combined_i %>%
    reshape2::dcast(., dataset + scenario + data_level + correct_level + correct_method + quant_method ~ label,
                    fun.aggregate = length) %>%
    mutate_if(is.integer, as.numeric) %>%
    mutate(mcc = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))) %>%
    mutate(recall = TP / (TP + FN),
           precision = TP / (TP + FP),
           f1 = 2 * precision * recall / (precision + recall))
  
  mcc_results_i <- list(df = df_combined_i, mcc = df_mcc_i)
  
  return(mcc_results_i)
})
saveRDS(mcc_objs, "./results/tables/mcc.rds")

