# ' Title: Figure3-quality assessment-CV/MCC.
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
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

rm(list = ls())


## global ---------------------
if (!interactive()) pboptions(type = "none")
options(mc.cores = 4)
# dir2work <- "/app/results"
dir2work <- "./results"
if (!dir.exists(dir2work)) dir.create(dir2work)
dir2save1 <- paste(dir2work, "/figures", sep = "")
if (!dir.exists(dir2save1)) dir.create(dir2save1)
dir2save2 <- paste(dir2work, "/tables", sep = "")
if (!dir.exists(dir2save2)) dir.create(dir2save2)
set.seed(20250630)


## Set plot parameters ---------------------
dictLabelsScenario <- c("Quartet-B", "Quartet-C", "Simulated-B", "Simulated-C")
names(dictLabelsScenario) <- c("quartet_balanced", "quartet_confounded",
                               "simulated_balanced", "simulated_confounded")

dictColorsLevel <- c("#E41A1C", "#984EA3", "#377EB8", "#006D2C")
names(dictColorsLevel) <- c("Uncorrected", "Precursor-corrected", "Peptide-corrected", "Protein-corrected")

dictColorsSample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#999", "#2171B5", "#CB181D")
names(dictColorsSample) <- c("D5", "D6", "F7", "M8", "Group1", "Group2", "Group3")

dictColorsMethod <- brewer.pal(7, "Accent")
names(dictColorsMethod) <- c("Ratio", "Med-C", "Combat", "RUV-III-C", "Harmony", "WaveICA2", "NormAE")

dictColorsQuantMethods <- c("#66C2A5", "#8DA0CB","#FC8D62")
names(dictColorsQuantMethods) <- c("maxlfq", "toppep3", "ibaq")

dictShapesBatch <- c(12:17, 1:3)
names(dictShapesBatch) <- c("L1-1", "L2-1", "L3", "L1-2", "L2-2", "L4",
                            "Batch1", "Batch2", "Batch3")

dictLabelsBatch <- c("L1-1", "L2-1", "L3", "L1-2", "L2-2", "L4", "Batch1", "Batch2", "Batch3")
names(dictLabelsBatch) <- c("DDA_APT", "DDA_FDU", "DDA_NVG", "DIA_APT", "DIA_BGI", "DIA_FDU",
                            "batch1", "batch2", "batch3")

dictLabelsCorrectMethods <- c("Ratio", "Med-C", "Combat", "RUV-III-C", "Harmony", "WaveICA2", "NormAE")
names(dictLabelsCorrectMethods) <- c("ratio", "median", "combat", "ruv", "harmony", "waveica", "normae")


## figure 3a ------------------------
cv_objs <- readRDS("./results/tables/cv.rds")
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  tmp_table <- tmp_table %>%
    select(cv, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
    mutate(label = ifelse(correct_method %in% "log",
                          "Uncorrected",
                          paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")))
})
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  # filter(cv < 1) %>%
  group_by(label, scenario) %>%
  summarise(cv_mean = mean(cv),cv_sd = sd(cv), cv_median = median(cv))

p_cv_box <- ggplot(sub_cv_tmp, aes(x = label, y = cv)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 10, name = "CV") +
  facet_wrap( ~ scenario, ncol = 4)

figure3a <- p_cv_box



## supplementary figure 5a ------------------------
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  sub_tmp <- tmp_table %>%
    select(cv, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
    mutate(label = ifelse(correct_method %in% "log",
                          "Uncorrected",
                          paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")))
  return(sub_tmp)
})
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  # filter(cv < 1) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "ratio_median_combat_ruv_harmony_waveica_normae", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_") %>%
  # filter(!correct_method %in% "normae") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel))) %>%
  mutate(x_axis = as.numeric(label))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!is.na(cv) | is.infinite(cv)) %>%
  # filter(cv < 1) %>%
  group_by(label, correct_method) %>%
  summarise(cv_mean = mean(cv),cv_sd = sd(cv), cv_median = median(cv))

p_cv_bar <- ggplot(sub_cv_tmp, aes(x = label, y = cv)) +
  geom_boxplot(aes(x = label, y = cv, fill = label), width = .7,
               position = position_dodge(width = .8)) +
  # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Uncorrected", "Protein-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = TRUE,
              # y_position = 1,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .9), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 10, name = "CV") +
  facet_wrap( ~ correct_method,
              labeller = as_labeller(dictLabelsCorrectMethods), 
              nrow = 1);p_cv_bar

suppl5a <- p_cv_bar



## supplementary figure 5b ------------------
mcc_objs <- readRDS("./results/tables/mcc.rds")
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) tmp_obj$mcc)
sub_mcc_tmp <- mcc_objs_tmp %>%
  rbindlist %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  group_by(quant_method, scenario) %>%
  mutate(cut_off = mcc[correct_method %in% "log"]) %>%
  ungroup %>%
  filter(!correct_method %in% "log") %>%
  mutate(label = paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_mcc_stat <- sub_mcc_tmp %>%
  group_by(label, scenario, quant_method) %>%
  summarise(mcc_mean = mean(mcc), mcc_sd = sd(mcc))

p_mcc_bar <- ggplot(sub_mcc_tmp, aes(x = label, y = mcc)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
  geom_signif(comparisons = list(c("Precursor-corrected", "Peptide-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = TRUE,
              y_position = .75,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 8, name = "MCC") +
  facet_grid(cols = vars(quant_method), rows = vars(scenario))

suppl5b <- p_mcc_bar


## figure 3b ------------------
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) tmp_obj$mcc)
sub_mcc_tmp <- mcc_objs_tmp %>%
  rbindlist %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  group_by(quant_method, scenario) %>%
  mutate(cut_off = mcc[correct_method %in% "log"]) %>%
  ungroup %>%
  mutate(label = ifelse(correct_method %in% "log", "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "ratio_combat_ruv_harmony_median_waveica", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_") %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_mcc_stat <- sub_mcc_tmp %>%
  group_by(label, scenario) %>%
  summarise(mcc_mean = mean(mcc), mcc_sd = sd(mcc))

p_mcc_bar <- ggplot(sub_mcc_tmp, aes(x = label, y = mcc)) +
  geom_col(aes(fill = label), alpha = .7, width = .05, position = position_dodge(width = .9)) +
  geom_point(aes(color = label), size = 5, position = position_dodge(width = .9)) +
  geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsLevel) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 8, name = "MCC") +
  facet_grid(cols = vars(correct_method), rows = vars(scenario))

figure3b <- p_mcc_bar;figure3b


## figure 3c -------------------
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) {
  df_tmp <- tmp_obj$df %>%
    distinct(class, value, estimate, dataset, scenario, data_level,
             correct_level, correct_method, quant_method, label) %>%
    mutate(label2 = ifelse(correct_method %in% "log", "Uncorrected",
                           paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")))
  
})
sub_ref_tmp <- mcc_objs_tmp %>%
  rbindlist %>%
  filter(scenario %in% "confounded") %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "log_ratio_median_combat_ruv_harmony_waveica_normae", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_") %>%
  filter(correct_method %in% c("ruv", "combat", "ratio")) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label2", ~ factor(., levels = names(dictColorsLevel)))

p_tmp <- ggplot(sub_ref_tmp, aes(x = value, y = estimate)) +
  geom_point(aes(color = label2), alpha = .5, shape = 16) +
  geom_smooth(aes(color = label2), method = "lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), size = 6) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_color_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 5, name = "Test") +
  scale_x_continuous(limits = c(-3.6, 3.6), name = "Reference") +
  facet_grid(cols = vars(label2), rows = vars(correct_method))

figure3c <- p_tmp;figure3c



## combind supplementary figure5 -----------------------
supp5 <- plot_grid(suppl5a, suppl5b,
                   nrow = 2, rel_heights = c(.6, 1),
                   labels = c("a", "b"), label_size = 24)

ggsave("./results/figures/extended_figure5.pdf",
       supp5, width = 12, height = 14)


## combind figure3 -----------------------
figure3 <- plot_grid(figure3a + theme(panel.spacing = unit(.5, "cm"),
                                      plot.margin = unit(c(0, 1.7, .5, 1.1), "cm")) +
                       labs(fill = "Correction level"),
                     figure3b + theme(plot.margin = unit(c(.5, .5, .5, 1.1), "cm")),
                     figure3c + theme(panel.spacing = unit(.5, "cm"),
                                      plot.margin = unit(c(.5, .5, .5, .7), "cm")),
                     nrow = 3, rel_heights = c(.8, 1.2, 2),
                     labels = c("a", "b", "c"), label_size = 24)

ggsave(paste("./results/figures/figure3.pdf", sep = ""),
       figure3, width = 14, height = 20)


## 弃CV-Rank trend ------------------------
## get CVs
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  tmp_table <- tmp_table %>%
    filter(data_level %in% "protein") %>%
    select(protein, cv, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
    mutate(label = ifelse(correct_method %in% "log",
                          "Uncorrected",
                          paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")))
})
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!is.na(cv)) %>%
  filter(!is.infinite(cv)) %>%
  # filter(cv < 1) %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "ratio_median_combat_ruv_harmony_waveica_normae", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_")

## get rank of intensities
all_files <- list.files("./data/expfiles", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("expdata_log", all_files)]
queried_files <- queried_files[grepl("protein/", queried_files)]
rank_tables <- pblapply(queried_files, function(file_id) {
  dataset <- str_extract(file_id,"(?<=.expfiles/).+?(?=\\/.)")
  scenario <- str_extract(file_id,"(?<=.(simulated|quartet)/).+?(?=\\/.)")
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  
  file_path <- str_extract(file_id, "^.*\\/(balanced|confounded|simulated)")
  df_meta <- fread(paste(file_path, "/meta.csv", sep = ""))
  
  df_cv <- sub_cv_tmp %>%
    filter(dataset %in% dataset) %>%
    filter(scenario %in% scenario) %>%
    filter(data_level %in% data_level) %>%
    filter(quant_method %in% quant_method) %>%
    filter(correct_level %in% correct_level) %>%
    filter(correct_method %in% correct_method) %>%
    distinct(protein)
  
  expr <- fread(file_id)
  
  df_features <- expr %>%
    tibble::rowid_to_column("feature") %>%
    distinct(feature, across(any_of(c("precursor", "peptide", "protein", "mz", "rt"))))
  
  df_expr <- expr %>%
    filter(protein %in% df_cv$protein) %>%
    tibble::rowid_to_column("feature") %>%
    select(feature, all_of(df_meta$run_id)) %>%
    reshape2::melt(., id = 1, variable.name = "run_id", na.rm = TRUE) %>%
    left_join(., df_meta, by = "run_id")
  
  df_rank <- df_expr %>%
    mutate_at("value", exp) %>%
    filter(value != 0, !is.na(value)) %>%
    group_by(feature) %>%
    summarise(mean_value = mean(value), .groups = "drop") %>%
    ungroup() %>%
    mutate(rank = rank(-mean_value)) %>%
    left_join(., df_features, by = "feature") %>%
    select(!feature) %>%
    select(any_of(c("precursor", "peptide", "protein", "mz", "rt")), everything()) %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(df_rank)
})
sub_rank_tmp <- rank_tables %>%
  rbindlist(.) %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "log_ratio_median_combat_ruv_harmony_waveica_normae", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_") %>%
  select(!correct_method) %>%
  select(!correct_level)

sub_cv_rank <- sub_cv_tmp %>%
  left_join(., sub_rank_tmp, relationship = "many-to-many") %>%
  filter(correct_method %in% "combat") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel))) %>%
  na.omit

dictAlphaTmp <- c(0.01, 0.01, 0.1, 0.1)
names(dictAlphaTmp) <- unique(sub_cv_rank$scenario)

p_tmp <- ggplot(sub_cv_rank, aes(x = rank, y = cv)) +
  # geom_point(aes(color = label, alpha = scenario), shape = 16) +
  geom_smooth(aes(color = label), method = "gam", formula = y ~ s(x, bs = "cs")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsLevel) +
  # scale_alpha_manual(values = dictAlphaTmp) +
  scale_y_continuous(n.breaks = 6, name = "CV (%)") +
  scale_x_continuous(name = "Rank") +
  facet_wrap(~ scenario, nrow = 1, scales = "free")

figure3b <- p_tmp


