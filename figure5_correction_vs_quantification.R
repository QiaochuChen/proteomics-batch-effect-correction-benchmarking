# ' Title: Figure5-correction vs quantification
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


## figure 5a -----------------------
cv_objs <- readRDS("./results/tables/cv.rds")
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  cv_table <- tmp_table %>%
    select(dataset, scenario, data_level, correct_level, correct_method, quant_method, cv)
  return(cv_table)
})

sub_cv <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  filter(!grepl("_log", correct_method)) %>%
  group_by(dataset, scenario, correct_method, quant_method, data_level, correct_level) %>%
  dplyr::summarise_at("cv", mean) %>%
  ungroup() %>%
  mutate_at("scenario", Hmisc::capitalize) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate(label = apply(., 1, function(x) {
    if (as.character(x[5]) %in% as.character(x[6])) {
      a <- "No-further-aggregation"
    } else if (as.character(x[5]) %in% "peptide" & as.character(x[6]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[5]) %in% "protein" & as.character(x[6]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[5]) %in% "protein" & as.character(x[6]) %in% "peptide") {
      a <- "Aggregated"
    }
  })) %>%
  mutate_at("label", ~ ifelse(correct_level %in% "protein", "No-further-aggregation_Aggregated", .)) %>%
  tidyr::separate_rows(label, sep = "_")

sub_cv_tmp <- sub_cv %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  filter(!correct_method %in% "log") %>%
  mutate(x_axis = as.numeric(correct_level))

p_cv_scatter <- ggplot(sub_cv_tmp, aes(x = x_axis, y = cv)) +
  # stat_summary(aes(color = label), fun = mean, geom = "line", size = 1) +
  geom_smooth(aes(group = label, color = label), na.rm = TRUE, alpha = .2, method = "loess") +
  stat_summary(aes(color = label), fun = mean, geom = "point", size = 5) +
  # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        panel.spacing = unit(1.5, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, 1.2), "cm")) +
  scale_color_manual(values = c("#7FBC41", "#C51B7D"), name = "Group") +
  scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                     labels = c("Precursor", "Peptide", "Protein")) +
  scale_y_continuous(n.breaks = 5, name = "CV") +
  facet_grid(cols = vars(scenario))

figure5a <- p_cv_scatter;figure5a


## supplementary figure 10a -----------------------
sub_cv_tmp <- sub_cv %>%
  filter(label %in% "No-further-aggregation") %>%
  mutate(label2 = ifelse(is.na(correct_method), "Uncorrected", "Corrected")) %>%
  mutate(x_axis = as.numeric(data_level))

p_cv_scatter <- ggplot(sub_cv_tmp, aes(x = x_axis, y = cv)) +
  # stat_summary(aes(color = label2), fun = mean, geom = "line", size = 1) +
  geom_smooth(aes(group = label2, color = label2), na.rm = TRUE, alpha = .2, method = "loess") +
  stat_summary(aes(color = label2), fun = mean, geom = "point", size = 5) +
  # stat_summary(aes(fill = label2), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        panel.spacing = unit(1.5, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, 1.2), "cm")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02"), name = "Group") +
  scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                     labels = c("Precursor", "Peptide", "Protein")) +
  scale_y_continuous(n.breaks = 5, name = "CV") +
  facet_grid(cols = vars(scenario))

supp10a <- p_cv_scatter;supp10a


## figure 5b -----------------------
pca_objs <- readRDS(paste("./results/tables/pca_snr.rds", sep = ""))
pcs_tables <- pblapply(pca_objs, function(pca_obj) {
  sub_tmp <- pca_obj$pcs_values
  if (pca_obj$dataset %in% c("quartet")) {
    sub_tmp <- sub_tmp %>% mutate(batch = paste(mode, lab, sep = "_"))
  }
  sub_pcs <- sub_tmp %>%
    dplyr::select(library, batch, sample, PC1, PC2) %>%
    mutate(snr = pca_obj$snr_results$snr,
           dataset = pca_obj$dataset,
           scenario = pca_obj$scenario,
           data_level = pca_obj$data_level,
           correct_level = pca_obj$correct_level,
           correct_method = pca_obj$correct_method,
           quant_method = pca_obj$quant_method)
  
  return(sub_pcs)
})
sub_pca_final <- rbindlist(pcs_tables)

sub_snr <- sub_pca_final %>%
  distinct(snr, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate_at("scenario", Hmisc::capitalize) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate(label = apply(., 1, function(x) {
    if (as.character(x[4]) %in% as.character(x[5])) {
      a <- "No-further-aggregation"
    } else if (as.character(x[4]) %in% "peptide" & as.character(x[5]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[4]) %in% "protein" & as.character(x[5]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[4]) %in% "protein" & as.character(x[5]) %in% "peptide") {
      a <- "Aggregated"
    }
  })) %>%
  mutate_at("label", ~ ifelse(correct_level %in% "protein", "No-further-aggregation_Aggregated", .)) %>%
  tidyr::separate_rows(label, sep = "_")

sub_snr_tmp <- sub_snr %>%
  filter(!correct_method %in% "log") %>%
  mutate(x_axis = as.numeric(correct_level))

p_snr_scatter <- ggplot(sub_snr_tmp, aes(x = x_axis, y = snr)) +
  # stat_summary(aes(color = label), fun = mean, geom = "line", size = 1) +
  geom_smooth(aes(group = label, color = label), na.rm = TRUE, alpha = .2, method = "loess") +
  stat_summary(aes(color = label), fun = mean, geom = "point", size = 5) +
  # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        panel.spacing = unit(1.5, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, 1.3), "cm")) +
  scale_color_manual(values = c("#7FBC41", "#C51B7D"), name = "Group") +
  scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                     labels = c("Precursor", "Peptide", "Protein")) +
  scale_y_continuous(n.breaks = 5, name = "SNR") +
  facet_grid(cols = vars(scenario))

figure5b <- p_snr_scatter;figure5b


## supplementary figure 10b -----------------------
sub_snr_tmp <- sub_snr %>%
  filter(label %in% "No-further-aggregation") %>%
  mutate(label2 = ifelse(is.na(correct_method), "Uncorrected", "Corrected")) %>%
  mutate(x_axis = as.numeric(data_level))

p_snr_scatter <- ggplot(sub_snr_tmp, aes(x = x_axis, y = snr)) +
  # stat_summary(aes(color = label2), fun = mean, geom = "line", size = 1) +
  geom_smooth(aes(group = label2, color = label2), na.rm = TRUE, alpha = .2, method = "loess") +
  stat_summary(aes(color = label2), fun = mean, geom = "point", size = 5) +
  # stat_summary(aes(fill = label2), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        panel.spacing = unit(1.5, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02"), name = "Group") +
  scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                     labels = c("Precursor", "Peptide", "Protein")) +
  scale_y_continuous(n.breaks = 5, name = "SNR") +
  facet_grid(cols = vars(scenario))

supp10b <- p_snr_scatter;supp10b


## figure 5c -----------------------
mcc_objs <- readRDS("./results/tables/mcc.rds")
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) tmp_obj$mcc)

sub_mcc <- mcc_objs_tmp %>%
  rbindlist %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate_at("scenario", Hmisc::capitalize) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate(label = apply(., 1, function(x) {
    if (as.character(x[3]) %in% as.character(x[4])) {
      a <- "No-further-aggregation"
    } else if (as.character(x[3]) %in% "peptide" & as.character(x[4]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[3]) %in% "protein" & as.character(x[4]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[3]) %in% "protein" & as.character(x[4]) %in% "peptide") {
      a <- "Aggregated"
    }
  })) %>%
  mutate_at("label", ~ ifelse(correct_level %in% "protein", "No-further-aggregation_Aggregated", .)) %>%
  tidyr::separate_rows(label, sep = "_")

sub_mcc_tmp <- sub_mcc %>%
  filter(!correct_method %in% "log") %>%
  mutate(x_axis = as.numeric(correct_level)) %>%
  select(x_axis, label, scenario, FN, FP, TN, TP) %>%
  reshape2::melt(., id = 1:3, variable.name = "metric")

p_list <- pblapply(unique(sub_mcc_tmp$scenario), function(scenario_id) {
  sub_mcc_tmp_i <- sub_mcc_tmp %>% filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_mcc_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 20),
          plot.margin = unit(c(.5, .5, 0, 2.8), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_mcc_scatter <- ggplot(sub_mcc_tmp_i, aes(x = x_axis, y = value)) +
    # stat_summary(aes(color = label), fun = mean, geom = "line", size = 1) +
    geom_smooth(aes(group = label, color = label), na.rm = TRUE, alpha = .2, method = "loess") +
    stat_summary(aes(color = label), fun = mean, geom = "point", size = 5) +
    # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
    theme_bw() +
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(0, .5, .5, .8), "cm")) +
    scale_color_manual(values = c("#7FBC41", "#C51B7D"),name = "Group") +
    scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                       labels = c("Precursor", "Peptide", "Protein")) +
    scale_y_continuous(n.breaks = 5,
                       labels = scales::label_number(scale = 1/1000, suffix = "k", accuracy = .1),
                       name = "Number of features") +
    facet_wrap(~ metric, ncol = 2, scales = "free_y")
  
  p_i <- plot_grid(p_title, p_mcc_scatter, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_i)
})

figure5c <- plot_grid(plotlist = p_list, ncol = 2);figure5c


## supplementary figure 10c -----------------------
sub_mcc_tmp <- sub_mcc %>%
  filter(label %in% "No-further-aggregation") %>%
  mutate(label2 = ifelse(is.na(correct_method), "Uncorrected", "Corrected")) %>%
  mutate(x_axis = as.numeric(data_level)) %>%
  select(x_axis, label2, scenario, FN, FP, TN, TP) %>%
  reshape2::melt(., id = 1:3, variable.name = "metric")

p_list <- pblapply(unique(sub_mcc_tmp$scenario), function(scenario_id) {
  sub_mcc_tmp_i <- sub_mcc_tmp %>% filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_mcc_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 20),
          plot.margin = unit(c(.5, .5, 0, 2.8), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_mcc_scatter <- ggplot(sub_mcc_tmp_i, aes(x = x_axis, y = value)) +
    # stat_summary(aes(color = label2), fun = mean, geom = "line", size = 1) +
    geom_smooth(aes(group = label2, color = label2), na.rm = TRUE, alpha = .2, method = "loess") +
    stat_summary(aes(color = label2), fun = mean, geom = "point", size = 5) +
    # stat_summary(aes(fill = label2), fun = mean, geom = "bar", width = .7) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
    theme_bw() +
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(0, .5, .5, .8), "cm")) +
    scale_color_manual(values = c("#1B9E77", "#D95F02"),name = "Group") +
    scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                       labels = c("Precursor", "Peptide", "Protein")) +
    scale_y_continuous(n.breaks = 5, limits = c(0, NA),
                       labels = scales::label_number(scale = 1/1000, suffix = "k", accuracy = .1),
                       name = "Number of features") +
    facet_wrap(~ metric, ncol = 2, scales = "free_y")
  
  p_i <- plot_grid(p_title, p_mcc_scatter, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_i)
})

supp10c <- plot_grid(plotlist = p_list, ncol = 2);supp10c


## figure 5d -----------------------
sub_pvca_final <- data.table::fread("./results/tables/pvca.csv")

sub_pvca <- sub_pvca_final %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate(label2 = apply(., 1, function(x) {
    if (as.character(x[1] %in% "sample")) {
      return("Biological")
    } else if (as.character(x[1]) %in% "Residual") {
      "Residual"
    } else if (grepl(":sample", x[1])|grepl("sample:", x[1])) {
      "Biological:Technical"
    } else {
      "Technical"
    }
  })) %>%
  group_by(dataset, scenario, correct_method, quant_method, label2, data_level, correct_level) %>%
  dplyr::summarise(proportion = sum(Proportion)) %>%
  ungroup() %>%
  mutate_at("label2", ~ factor(., levels = c("Biological", "Technical", "Biological:Technical", "Residual"))) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("scenario", Hmisc::capitalize) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate(label = apply(., 1, function(x) {
    if (as.character(x[6]) %in% as.character(x[7])) {
      a <- "No-further-aggregation"
    } else if (as.character(x[6]) %in% "peptide" & as.character(x[7]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[6]) %in% "protein" & as.character(x[7]) %in% "precursor") {
      a <- "Aggregated"
    } else if (as.character(x[6]) %in% "protein" & as.character(x[7]) %in% "peptide") {
      a <- "Aggregated"
    }
  })) %>%
  mutate_at("label", ~ ifelse(correct_level %in% "protein", "No-further-aggregation_Aggregated", .)) %>%
  tidyr::separate_rows(label, sep = "_")

sub_pvca_tmp <- sub_pvca %>%
  filter(!correct_method %in% "log") %>%
  mutate(x_axis = as.numeric(correct_level))

p_list <- pblapply(unique(sub_pvca_tmp$scenario), function(scenario_id) {
  sub_pvca_tmp_i <- sub_pvca_tmp %>% filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_pvca_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 20),
          plot.margin = unit(c(.5, .5, 0, 2.8), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_pvca_scatter <- ggplot(sub_pvca_tmp_i, aes(x = x_axis, y = proportion)) +
    # stat_summary(aes(color = label), fun = mean, geom = "line", size = 1) +
    geom_smooth(aes(group = label, color = label), na.rm = TRUE, alpha = .2) +
    stat_summary(aes(color = label), fun = mean, geom = "point", size = 5) +
    # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(0, .5, .5, 1), "cm")) +
    scale_color_manual(values = c("#7FBC41", "#C51B7D"), name = "Group") +
    scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                       labels = c("Precursor", "Peptide", "Protein")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25),
                       labels = ~ . * 100, name = "PVCA (%)") +
    facet_wrap(~ label2, ncol = 2, scales = "free_y")
  
  p_i <- plot_grid(p_title, p_pvca_scatter, nrow = 2, rel_heights = c(.1, 1))

  return(p_i)
  
})

figure5d <- plot_grid(plotlist = p_list, ncol = 2);figure5d


## supplementary figure 10d -----------------------
sub_pvca_tmp <- sub_pvca %>%
  filter(label %in% "No-further-aggregation") %>%
  mutate(label3 = ifelse(is.na(correct_method), "Uncorrected", "Corrected")) %>%
  mutate(x_axis = as.numeric(data_level))

p_list <- pblapply(unique(sub_pvca_tmp$scenario), function(scenario_id) {
  sub_pvca_tmp_i <- sub_pvca_tmp %>% filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_pvca_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 20),
          plot.margin = unit(c(.5, .5, 0, 2.8), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_pvca_scatter <- ggplot(sub_pvca_tmp_i, aes(x = x_axis, y = proportion)) +
    # stat_summary(aes(color = label3), fun = mean, geom = "line", size = 1) +
    geom_smooth(aes(group = label3, color = label3), na.rm = TRUE, alpha = .2, method = "loess") +
    stat_summary(aes(color = label3), fun = mean, geom = "point", size = 5) +
    # stat_summary(aes(fill = label3), fun = mean, geom = "bar", width = .7) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(0, .5, .5, 1), "cm")) +
    scale_color_manual(values = c("#1B9E77", "#D95F02"),name = "Group") +
    scale_x_continuous(limits = c(.8, 3.2), breaks = 1:3,
                       labels = c("Precursor", "Peptide", "Protein")) +
    scale_y_continuous(limits = c(-.1, 1.15), breaks = seq(0, 1, .25),
                       labels = ~ . * 100, name = "PVCA (%)") +
    facet_wrap(~ label2, ncol = 2, scales = "free_y")
  
  p_i <- plot_grid(p_title, p_pvca_scatter, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_i)
})

supp10d <- plot_grid(plotlist = p_list, ncol = 2);supp10d


## combined supplementary figure10 -----------------------
supp10ab <- plot_grid(supp10a, supp10b,
                       ncol = 2,
                       labels = c("a", "b"), label_size = 24)
supp10 <- plot_grid(supp10ab, supp10c, supp10d,
                     nrow = 3, rel_heights = c(.5, 1, 1.05),
                     labels = c("", "c", "d"), label_size = 24)

ggsave("./results/figures/extended_figure10.pdf", supp10, width = 14, height = 16)


## combined figure5 -----------------------
figure5ab <- plot_grid(figure5a, figure5b,
                     ncol = 2,
                     labels = c("a", "b"), label_size = 24)
figure5 <- plot_grid(figure5ab, figure5c, figure5d,
                   nrow = 3, rel_heights = c(.5, 1, 1),
                   labels = c("", "c", "d"), label_size = 24)

ggsave("./results/figures/figure5.pdf", figure5, width = 14, height = 16)


## supplementary data table -----------------------
cv_objs <- readRDS("./results/tables/cv.rds")
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  cv_table <- tmp_table %>%
    select(dataset, scenario, data_level, correct_level, correct_method, quant_method, cv)
  return(cv_table)
})

sub_cv <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(dataset, scenario, correct_method, quant_method, data_level, correct_level) %>%
  dplyr::summarise(cv_mean = mean(cv), cv_sd = sd(cv)) %>%
  ungroup() 

sub_test <- sub_cv %>%
  filter(data_level %in% "protein") %>%
  filter(!grepl("log", correct_method)) %>%
  filter(quant_method %in% "maxlfq") %>%
  reshape2::dcast(., dataset + scenario + correct_method + quant_method ~ correct_level,
                  value.var = "cv_mean") %>%
  mutate(label = ifelse(protein < precursor & protein < peptide, "yes", "no")) %>%
  reshape2::dcast(., dataset + scenario ~ label, value.var = "label")
  
openxlsx::write.xlsx(sub_cv, "./results/tables/exdended_data_table1.xlsx")

pca_objs <- readRDS(paste("./results/tables/pca_snr.rds", sep = ""))
pcs_tables <- pblapply(pca_objs, function(pca_obj) {
  sub_tmp <- pca_obj$pcs_values
  sub_pcs <- data.frame(snr = pca_obj$snr_results$snr,
           dataset = pca_obj$dataset,
           scenario = pca_obj$scenario,
           data_level = pca_obj$data_level,
           correct_level = pca_obj$correct_level,
           correct_method = pca_obj$correct_method,
           quant_method = pca_obj$quant_method)
  
  return(sub_pcs)
})
sub_snr <- rbindlist(pcs_tables)

sub_test <- sub_snr %>%
  filter(data_level %in% "protein") %>%
  filter(!grepl("log", correct_method)) %>%
  filter(quant_method %in% "maxlfq") %>%
  reshape2::dcast(., dataset + scenario + correct_method + quant_method ~ correct_level,
                  value.var = "snr") %>%
  mutate(label = ifelse(protein > precursor & protein > peptide, "yes", "no")) %>%
  reshape2::dcast(., dataset + scenario ~ label, value.var = "label")

openxlsx::write.xlsx(sub_snr, "./results/tables/exdended_data_table2.xlsx")

mcc_objs <- readRDS("./results/tables/mcc.rds")
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) tmp_obj$mcc)
sub_mcc <- rbindlist(mcc_objs_tmp)

sub_test <- sub_mcc %>%
  filter(data_level %in% "protein") %>%
  filter(!grepl("log", correct_method)) %>%
  # filter(quant_method %in% "maxlfq") %>%
  reshape2::dcast(., dataset + scenario + correct_method + quant_method ~ correct_level,
                  value.var = "mcc") %>%
  mutate(label = ifelse(protein > precursor & protein > peptide, "yes", "no"))%>%
  reshape2::dcast(., dataset + scenario ~ label, value.var = "label")

openxlsx::write.xlsx(sub_mcc, "./results/tables/exdended_data_table3.xlsx")

sub_pvca <- data.table::fread("./results/tables/pvca.csv")

sub_test <- sub_pvca %>%
  filter(data_level %in% "protein") %>%
  filter(!grepl("log", correct_method)) %>%
  # filter(quant_method %in% "maxlfq") %>%
  filter(`Random Effect` %in% "sample") %>%
  reshape2::dcast(., dataset + scenario + correct_method + quant_method ~ correct_level,
                  value.var = "Proportion") %>%
  mutate(label = ifelse(protein > precursor & protein > peptide, "yes", "no"))%>%
  reshape2::dcast(., dataset + scenario ~ label, value.var = "label")

openxlsx::write.xlsx(sub_pvca, "./results/tables/exdended_data_table4.xlsx")




