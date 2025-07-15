# ' Title: Data simulation. 参考文献http://dx.doi.org/10.1093/bib/bbaf168
# ' Author: Qiaochu Chen

library(pbapply)
library(dplyr)
library(data.table)
library(stringr)
library(fitdistrplus)

rm(list = ls())


## 0-1 基于真实观测数据的two-way anova结果确定各个类别特征比例 ------------- 
all_files <- list.files("./data/expfiles/quartet", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("anova_twoway_test", all_files)]
queried_files <- queried_files[grepl("balanced/protein", queried_files)]
prop_tables <- pblapply(queried_files, function(file_id) {
  
  df_anova_i <- fread(file_id)
  df_anova_test <- df_anova_i %>%
    reshape2::dcast(., feature ~ Effect, value.var = "p<.05")
  
  ## 总特征数
  n_total <- df_anova_test %>%
    pull(feature) %>%
    unique %>%
    length
  
  ## 类别1-纯生物学效应特征数
  n_biol <- df_anova_test %>%
    filter(sample %in% "*" & batch %in% "" & `sample:batch` %in% "") %>%
    pull(feature) %>%
    unique %>%
    length
  
  ## 类别2-生物学+批次效应特征数
  n_biol_tech <- df_anova_test %>%
    filter((sample %in% "*" & batch %in% "*")|`sample:batch` %in% "*") %>%
    pull(feature) %>%
    unique %>%
    length
  
  ## 类别3-纯批次效应特征数
  n_tech <- df_anova_test %>%
    filter(batch %in% "*" & sample %in% "" & `sample:batch` %in% "") %>%
    pull(feature) %>%
    unique %>%
    length
  
  ## 类别4-恒定无关特征数
  n_housekeep <- df_anova_test %>%
    filter(sample %in% "" & batch %in% "" & `sample:batch` %in% "") %>%
    pull(feature) %>%
    unique %>%
    length
  
  ## 各类别占比
  df_prop_i <- data.frame(
    type = "protein",
    quantification_method = str_extract(file_id, "(?<=balanced\\/).+(?=\\/anova_twoway_test)"),
    n_biol = n_biol,
    n_biol_tech = n_biol_tech,
    n_tech = n_tech,
    n_housekeep = n_housekeep,
    n_total = n_total
  )
  
  return(df_prop_i)
})
df_prop <- rbindlist(prop_tables)
df_prop_final <- df_prop %>%
  mutate_at(3:7, ~ round(. / n_total, digits = 4))
mean_props <- colMeans(df_prop_final[, 3:7])


## 0-2 基于真实观测数据模拟前体-肽段-蛋白质映射关系 ------------
df_quant_wide <- fread("./data/expfiles/quartet/balanced/precursor/expdata_intensity.csv")
df_pre2pep2pro <- df_quant_wide %>% distinct(precursor, peptide, protein)
df_stat <- df_pre2pep2pro %>%
  filter(!protein %in% "") %>%
  filter(!grepl(";", protein)) %>%
  group_by(protein) %>%
  summarise(n_pep = length(unique(peptide)),
            n_pre = length(unique(precursor)))
mean(df_stat$n_pep)
mean(df_stat$n_pre)

df_stat <- df_pre2pep2pro %>%
  group_by(peptide) %>%
  summarise(n_pre = length(unique(precursor)))
mean(df_stat$n_pre)


## 0-3 基于真实观测数据确定分布参数 ----------------
df_quant_wide <- fread("./data/expfiles/quartet/balanced/precursor/expdata_intensity.csv")
feature_means <- df_quant_wide %>%
  dplyr::select(starts_with("D")) %>%
  mutate_all(log) %>%
  rowMeans(., na.rm = TRUE)

## 估计 gamma 分布参数
fit1 <- fitdist(feature_means, distr = "gamma", method = "mle")
summary(fit1) ## 得到结果shape = 150.2, rate = 8.4, AIC = 390608.3, BIC = 390627.5 
par(mar = c(4, 4, 2, 1))
plot(fit1)

## 估计 norm 分布参数
fit2 <- fitdist(feature_means, distr = "norm", method = "mle")
summary(fit2) ## 得到结果mean = 17.9, sd = 1.5, AIC = 392005.8, BIC = 392025 
plot(fit2)


## 1. 设置模拟参数 ---------------
set.seed(20250520)
## 样本设置: 生物学分组、批次分组、重复测量次数
n_biol_group <- 3
n_batch <- 3
n_rep <- 3

## 特征设置1: 各类别蛋白质数目
n_total_protein <- 300
n_class1_protein <- 0 ## 类别1-纯生物学效应特征数
n_class2_protein <- 200 ## 类别2-生物学+批次效应特征数
n_class3_protein <- 96 ## 类别3-纯批次效应特征数
n_class4_protein <- 4 ## 类别4-恒定无关特征数

## 特征设置2: 前体、肽段数目
n_total_peptide <- n_total_protein * 10
n_total_precursor <- n_total_peptide * 1.65

## 基础表达量设置(真实数据更符合gamma分布)
gamma_shape <- 150.2
gamma_rate <- 8.4

## 生物学效应设置
gamma_shape_logfc <- 10
gamma_rate_logfc <- 6

## 随机误差设置
epsilon_noise <- 0.5

## 样本间缩放因子设置
epsilon_sample <- 0.2

## 批次效应设置
sigma_batch_mean <- 2 ## 批次均值效应
sigma_batch_spec <- 0.5 ## 批次特异效应


## 2. 模拟前体-肽段-蛋白质映射关系 --------------
## 确保每个蛋白质至少分配了1个肽段
df_pep2pro1 <- data.frame(peptide = 1:n_total_protein,
                          protein = 1:n_total_protein)
## 随机分配剩余肽段
df_pep2pro2 <- data.frame(
  peptide = (1+n_total_protein):n_total_peptide,
  protein = sample(1:n_total_protein,
                      n_total_peptide - n_total_protein, replace = TRUE))

## 完成模拟肽段-蛋白质映射关系
df_pep2pro <- rbind(df_pep2pro1, df_pep2pro2)

## 确保每个肽段至少分配了1个前体
df_pre2pep1 <- data.frame(precursor = 1:n_total_peptide,
                          peptide = 1:n_total_peptide)

## 随机分配剩余前体
df_pre2pep2 <- data.frame(
  precursor = (1+n_total_peptide):n_total_precursor,
  peptide = sample(1:n_total_peptide,
                   n_total_precursor - n_total_peptide, replace = TRUE))

## 完成模拟前体-肽段映射关系
df_pre2pep <- rbind(df_pre2pep1, df_pre2pep2)

## 完成模拟前体-肽段-蛋白质映射关系
df_pre2pep2pro <- df_pre2pep %>%
  full_join(., df_pep2pro, by = "peptide") %>%
  mutate(feature = paste("pre", precursor, "_pep", peptide, "_pro", protein, sep = "")) %>%
  mutate(class = sapply(protein, function(x) {
    if (x %in% 0: n_class1_protein) {
      a <- 1
    } else if (x %in% (n_class1_protein + 1):(n_class1_protein + n_class2_protein)) {
      a <- 2
    } else if (x %in% (n_class1_protein + n_class2_protein + 1):(n_class1_protein + n_class2_protein + n_class3_protein)) {
      a <- 3
    } else if (x %in% (n_class3_protein + 1):(n_total_protein)) {
      a <- 4
    }
  }))


## 3. 模拟基础表达量 --------------
## 设计矩阵结构
df_meta <- data.frame(sample = paste("group", rep(1:n_biol_group, n_batch * n_rep), sep = ""),
                      batch = paste("batch", rep(1:n_batch, each = n_biol_group * n_rep), sep = ""),
                      rep = paste("rep", rep(rep(1:n_rep, each = n_biol_group), n_batch), sep = ""),
                      order = 1:(n_biol_group * n_batch * n_rep)) %>%
  mutate(run_id = paste(batch, sample, rep, sep = "_")) %>%
  dplyr::select(run_id, everything())
fwrite(df_meta, "./data/expfiles/simulated/balanced/meta.csv")

data_matrix <- matrix(0, nrow = nrow(df_pre2pep2pro), ncol = nrow(df_meta),
                      dimnames = list(df_pre2pep2pro$feature, df_meta$run_id))

## 模拟基础表达量
log_psi <- rgamma(nrow(data_matrix), shape = gamma_shape, rate = gamma_rate)

## 将log_psi赋给data_matrix
for (i in 1:ncol(data_matrix)) data_matrix[, i] <- log_psi

## 储存中间结果
data_df <- data_matrix %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  full_join(df_pre2pep2pro, ., by = "feature") %>%
  dplyr::select(precursor, peptide, protein, starts_with("batch"))
fwrite(data_df, "./data/expfiles/simulated/balanced/precursor/1_expdata_gamma_log.csv")


## 4. 模拟含生物学效应的特征表达量并添加随机测量误差 --------------
## 设置各生物分组下的特征矩阵结构
log_rho <- matrix(0, nrow = n_total_protein, ncol = n_biol_group,
                  dimnames = list(1:n_total_protein,
                                  c("group1", "group2/group1", "group3/group1")))

## 以Group1为基准, 3个生物学分组对应8种表达模式, 每种模式分配15-16个蛋白质
log_rho <- log_rho %>%
  as.data.frame %>%
  tibble::rowid_to_column("protein") %>%
  mutate(biol_mode = sapply(protein, function(x) {
    if (x %in% 1:15) {
      a <- 1
    } else if (x %in% 16:30) {
      a <- 2
    } else if (x %in% 31:36) {
      a <- 3
    } else if (x %in% 37:51) {
      a <- 4
    } else if (x %in% 52:66) {
      a <- 5
    } else if (x %in% 67:82) {
      a <- 6
    } else if (x %in% 83:98) {
      a <- 7
    } else if (x %in% 99:123) {
      a <- 8
    } else a <- 0
  }))

## 遍历所有蛋白质
for (i in 1:nrow(log_rho)) {
  if (log_rho$biol_mode[i] == 0) {
    ## 模式0-无生物学效应
    log_fc_group2 <- 0
    log_fc_group3 <- 0
    
  } else if (log_rho$biol_mode[i] == 1) {
    ## 模式1-Group2上调, Group3不变; 15个蛋白质
    log_rho$`group2/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- 0
    
  } else if (log_rho$biol_mode[i] == 2) {
    ## 模式2-Group2下调, Group3不变; 15个蛋白质
    log_rho$`group2/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- 0
    
  } else if (log_rho$biol_mode[i] == 3) {
    ## 模式3-Group2不变, Group3上调; 16个蛋白质
    log_rho$`group2/group1`[i] <- 0
    log_rho$`group3/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  } else if (log_rho$biol_mode[i] == 4) {
    ## 模式4-Group2不变, Group3下调; 15个蛋白质
    log_rho$`group2/group1`[i] <- 0
    log_rho$`group3/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  } else if (log_rho$biol_mode[i] == 5) {
    ## 模式5-Group2上调, Group3上调; 15个蛋白质
    log_rho$`group2/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  } else if (log_rho$biol_mode[i] == 6) {
    ## 模式6-Group2上调, Group3下调; 16个蛋白质
    log_rho$`group2/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  } else if (log_rho$biol_mode[i] == 7) {
    ## 模式7-Group2下调, Group3上调; 16个蛋白质
    log_rho$`group2/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  } else if (log_rho$biol_mode[i] == 8) {
    ## 模式8-Group2下调, Group3下调; 15个蛋白质
    log_rho$`group2/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    log_rho$`group3/group1`[i] <- -rgamma(1, shape = gamma_shape_logfc, rate = gamma_rate_logfc)
    
  }
  
}

## 合并蛋白质-肽段-前体映射关系
log_rho_all <- log_rho %>%
  full_join(., df_pre2pep2pro, by = "protein") %>%
  dplyr::select(feature, `group2/group1`, `group3/group1`)

## 将log_rho赋给data_matrix
data_df_tmp <- data_matrix %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  reshape2::melt(., id = 1, variable.name = "run_id") %>%
  left_join(., df_meta, by = "run_id") %>%
  left_join(log_rho_all, ., by = "feature") %>%
  mutate(value_new = case_when(sample %in% "group1"~value,
                               sample %in% "group2"~value + `group2/group1`,
                               sample %in% "group3"~value + `group3/group1`)) %>%
  mutate(value_new2 = rnorm(n(), mean = value_new, sd = epsilon_noise))

data_matrix_withbiol <- data_df_tmp %>%
  reshape2::dcast(., feature ~ run_id, value.var = "value_new2") %>%
  dplyr::select(feature, all_of(colnames(data_matrix))) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix

## 储存中间结果
data_df <- data_matrix_withbiol %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  full_join(df_pre2pep2pro, ., by = "feature") %>%
  dplyr::select(precursor, peptide, protein, starts_with("batch"))
fwrite(data_df, "./data/expfiles/simulated/balanced/precursor/2_expdata_withbiol_log.csv")


## 5. 模拟样本间效应的特征表达量 --------------
log_alpha <- rnorm(n_biol_group * n_batch * n_rep, 0, epsilon_sample)
data_matrix_withscale <- sweep(data_matrix_withbiol, 2, log_alpha, `+`)

## 储存中间结果
data_df <- data_matrix_withscale %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  full_join(df_pre2pep2pro, ., by = "feature") %>%
  dplyr::select(precursor, peptide, protein, starts_with("batch"))
fwrite(data_df, "./data/expfiles/simulated/balanced/precursor/3_expdata_withscal_log.csv")


## 6. 模拟含批次效应的特征表达量 --------------
## 合并蛋白质-肽段-前体映射关系
data_df_withscale <- data_matrix_withscale %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  reshape2::melt(., id = 1, variable.name = "run_id") %>%
  left_join(., df_pre2pep2pro, by = "feature") %>%
  left_join(., df_meta, by = "run_id")

## 遍历所有前体模拟批次均值效应
log_beta_all <- data_df_withscale %>%
  distinct(feature, class, batch) %>%
  mutate(log_beta = 0)
for (i in 1:nrow(log_beta_all)) {
  if (log_beta_all$class[i] %in% c(2, 3)) {## 类别2/类别3-有批次效应特征数
    log_beta_all$log_beta[i] <- rnorm(1, 0, sigma_batch_mean)
  }
}

## 遍历所有前体模拟批次特异效应
log_omega_all <- data_df_withscale %>%
  left_join(., log_beta_all, by = c("feature", "class", "batch")) %>%
  mutate(log_omega = 0)
for (i in 1:nrow(log_omega_all)) {
  if (log_omega_all$class[i] %in% c(2, 3)) {## 类别2/类别3-有批次效应特征数
    mu_batch_spec <- log_omega_all$log_beta[i]
    log_omega_all$log_omega[i] <- rnorm(1, mu_batch_spec, sigma_batch_spec)
  }
}

## 将log_omega赋给data_matrix
data_matrix_withbatch <- log_omega_all %>%
  mutate(value_new = value + log_omega) %>%
  reshape2::acast(., feature ~ run_id, value.var = "value_new")

## 储存最终结果
data_df <- data_matrix_withbatch %>%
  as.data.frame %>%
  tibble::rownames_to_column("feature") %>%
  full_join(df_pre2pep2pro, ., by = "feature") %>%
  dplyr::select(precursor, peptide, protein, starts_with("batch"))
fwrite(data_df, "./data/expfiles/simulated/balanced/precursor/4_expdata_withbatch_log.csv")

## balanced design
fwrite(data_df, "./data/expfiles/simulated/balanced/precursor/expdata_log.csv")

## confounded design
df_meta <- fread("./data/expfiles/simulated/balanced/meta.csv")
df_meta_c <- df_meta %>%
  filter(grepl("group1|(batch1_group3)|(batch2_group2)|(batch3)", run_id)) %>%
  mutate(order = 1:21)
fwrite(df_meta_c, "./data/expfiles/simulated/confounded/meta.csv")

data_df <- fread("./data/expfiles/simulated/balanced/precursor/expdata_log.csv")
data_df_c <- data_df %>% select(precursor, peptide, protein, all_of(df_meta_c$run_id))
fwrite(data_df_c, "./data/expfiles/simulated/confounded/precursor/expdata_log.csv")

## 储存真实答案
df_mapping <- df_pre2pep2pro %>%
  left_join(., log_rho_all, by = "feature") %>%
  dplyr::select(!feature)

fwrite(df_mapping, "./data/expfiles/simulated/mapping.csv")
