# ' Title: Figure4-quality assessment-SNR.
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

p_format <- function(p_value) {
  # 设置一个非常小的阈值
  threshold <- .0001
  # 如果P值小于阈值，则使用 "<" 符号
  if (p_value < threshold) {
    return(paste0("**** (P < 0.0001)"))
  } else if (p_value < .001) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("*** (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else if (p_value < .01) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("** (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else if (p_value < .05) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("* (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else {
    return(paste0("NS (P = ", format(p_value, digits = 2), ")", sep = ""))
  }
}


## figure 4a ------------------------
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

sub_snr_tmp <- sub_pca_final %>%
  distinct(snr, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  mutate(label = ifelse(correct_method %in% "log",
                        "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_snr_stat <- sub_snr_tmp %>%
  group_by(label, scenario) %>%
  summarise(snr_mean = mean(snr),cv_sd = sd(snr), cv_median = median(snr),
            n = length(unique(correct_method)),
            n_total = length(snr))

p_snr_box <- ggplot(sub_snr_tmp, aes(x = label, y = snr)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  geom_signif(comparisons = list(c("Precursor-corrected", "Peptide-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = p_format,
              step_increase = .15,
              tip_length = .01,
              textsize = 2,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, 0, .5), "cm")) +
  scale_fill_manual(values = dictColorsLevel, name = "Correction level") +
  scale_y_continuous(name = "SNR", n.breaks = 5, 
                     expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap( ~ scenario, scales = "free_y", ncol = 2)

figure4a <- p_snr_box


## supplementary figure 6a ------------------
sub_snr_tmp <- sub_pca_final %>%
  distinct(snr, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  mutate(label = ifelse(correct_method %in% "log",
                        "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("correct_method", ~ ifelse(. %in% "log", "ratio_median_combat_ruv_harmony_waveica_normae", .)) %>%
  tidyr::separate_rows(correct_method, sep = "_") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_snr_stat <- sub_snr_tmp %>%
  group_by(correct_method, label) %>%
  summarise(snr_mean = mean(snr), snr_sd = sd(snr), snr_median = median(snr))

p_snr_box <- ggplot(sub_snr_tmp, aes(x = label, y = snr)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 5,  name = "SNR") +
  facet_wrap(~ correct_method, nrow = 1);p_snr_box

suppl6a <- p_snr_box


## supplementary figure 6b ------------------
sub_snr_tmp <- sub_pca_final %>%
  distinct(snr, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  group_by(dataset, scenario, quant_method) %>%
  mutate(cut_off = snr[correct_method %in% "log"]) %>%
  ungroup %>%
  filter(!correct_method %in% "log") %>%
  mutate(label = paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_snr_stat <- sub_snr_tmp %>%
  group_by(quant_method, scenario) %>%
  summarise(snr_mean = mean(snr), mcc_sd = sd(snr))

p_snr_box <- ggplot(sub_snr_tmp, aes(x = label, y = snr)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
  geom_signif(comparisons = list(c("Precursor-corrected", "Peptide-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = p_format,
              # y_position = 20,
              step_increase = .15,
              tip_length = .01,
              textsize = 2,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = dictColorsLevel, name = "Correction level") +
  scale_y_continuous(n.breaks = 5, name = "SNR",
                     expand = expansion(mult = c(0.05, 0.15))) +
  facet_grid(cols = vars(scenario), rows = vars(quant_method));p_snr_box

suppl6b <- p_snr_box


## combined supplementary figure6 -----------------------
suppl6 <- plot_grid(suppl6a + theme(plot.margin = unit(c(.5, 1.5, .5, .5), "cm")),
                    suppl6b + theme(plot.margin = unit(c(.5, .5, .5, .5), "cm")),
                   nrow = 2, rel_heights = c(.5, 1),
                   labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure6.pdf", suppl6,
       height = 260, width = 210, units = "mm")


## figure 4b ------------------
sub_snr_tmp <- sub_pca_final %>%
  distinct(snr, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  group_by(dataset, scenario, quant_method) %>%
  mutate(cut_off = snr[correct_method %in% "log"]) %>%
  ungroup %>%
  filter(!correct_method %in% "log") %>%
  mutate(label = paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

p_list <- pblapply(unique(sub_snr_tmp$scenario), function(scenario_id) {
  
  sub_snr_tmp_i <- sub_snr_tmp %>%
    filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_snr_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_rect(linewidth = .7),
          plot.margin = unit(c(.9, 1.4, .5, 1.5), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_snr_bar <- ggplot(sub_snr_tmp_i, aes(x = reorder(correct_method, snr), y = snr)) +
    geom_col(aes(fill = label), alpha = .7, width = .1, position = position_dodge(width = .9)) +
    geom_point(aes(color = label), size = 3, position = position_dodge(width = .9)) +
    geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(size = 10, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_blank(),
          axis.title.y = element_text(size =12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0, .5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsLevel) +
    scale_fill_manual(values = dictColorsLevel) +
    scale_y_continuous(n.breaks = 5, name = "SNR") +
    facet_grid(cols = vars(label), rows = vars(quant_method))
  
  if (scenario_id %in% "Quartet-B") {
    p_snr_bar <- p_snr_bar + theme(plot.margin = unit(c(0, .5, 1.3, .5), "cm"))
    p_snr_i <- plot_grid(p_title, p_snr_bar, nrow = 2, rel_heights = c(.11, 1))
  } else {
    p_snr_i <- plot_grid(p_title, p_snr_bar, nrow = 2, rel_heights = c(.15, 1))
  }
  
  return(p_snr_i)
})

figure4b <- p_list[[1]];figure4b


## supplementary figure 7a ------------------
supp7a <-  p_list[[2]]


## supplementary figure 7b ------------------
supp7b <-  p_list[[3]]


## supplementary figure 7b ------------------
supp7c <-  p_list[[4]]


## combined supplementary figure7 -----------------------
supp7 <- plot_grid(supp7a, supp7b, supp7c,
                   nrow = 3, labels = c("a", "b", "c"), label_size = 20)

ggsave("./results/figures/extended_figure7.pdf",
       supp7, width = 210, height = 297, units = "mm")


## supplementary figure 8a ------------------------
sub_pca_tmp <- sub_pca_final %>%
  filter(!correct_method %in% "log") %>%
  filter(!correct_method %in% "normae") %>%
  filter(data_level %in% "protein") %>%
  filter(correct_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate(title = sprintf("%s\n(%s-corrected)", quant_method, Hmisc::capitalize(correct_level))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate(label = sprintf("%s\n(SNR = %.2f)", correct_method, snr)) %>%
  mutate_at("batch", ~ dictLabelsBatch[.]) %>%
  mutate_at("sample", ~ factor(Hmisc::capitalize(.),
                               levels = c("D5", "D6", "F7", "M8",
                                          "Group1", "Group2", "Group3"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

lg_pca_shape1 <- ggplot(sub_pca_tmp %>% filter(grepl("Quartet", scenario))) +
  geom_point(aes(x = PC1, y = PC2, shape = batch), size = 4) +
  scale_shape_manual(values = dictShapesBatch) +
  labs(shape = "Batch") +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

lg_pca_color1 <- ggplot(sub_pca_tmp %>% filter(grepl("Quartet", scenario))) +
  geom_point(aes(x = PC1, y = PC2, color = sample), size = 4) +
  scale_color_manual(values = dictColorsSample) +
  labs(color = "Sample") +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

p_title <- ggplot(sub_pca_tmp) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(linewidth = .7),
        plot.margin = unit(c(.5, .5, .5, 1.6), "cm"))+
  facet_grid(cols = vars(title))

p_pca_list <- mclapply(unique(sub_pca_tmp$scenario), function(scenario_id) {
  
  sub_pca_tmp_i <- sub_pca_tmp %>%
    filter(scenario %in% scenario_id)
  
  sub_tmp <- sub_pca_tmp_i %>% distinct(correct_method, label)
  dictLabelsTmp <- setNames(sub_tmp$label, sub_tmp$correct_method)
  
  p_title_i <- ggplot(sub_pca_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 10),
          strip.background = element_rect(linewidth = .7),
          plot.margin = unit(c(.5, .5, .5, 1.6), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_tmp <- ggplot(sub_pca_tmp_i,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, shape = batch), size = 1) +
    theme_bw() +
    theme(legend.position = "none",
          legend.text = element_text(size = 10),
          title = element_text(size = 10),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 9),
          strip.background = element_blank(),
          plot.margin = unit(c(0, .5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsSample) +
    scale_shape_manual(values = dictShapesBatch) +
    labs(color = "Sample", shape = "Batch") +
    facet_wrap(vars(correct_method), labeller = as_labeller(dictLabelsTmp),
               scales = "free", ncol = 3)
  
  p_i <- plot_grid(p_title_i, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_i)
})

p_pca_all <- plot_grid(plotlist = p_pca_list, ncol = 2)
p_pca_final <- plot_grid(p_title, p_pca_all, nrow = 2, rel_heights = c(.11, 1))
lg_pca1 <- plot_grid(get_legend(lg_pca_shape1), get_legend(lg_pca_color1), ncol = 2)
supp8a <- plot_grid(p_pca_final, lg_pca1, nrow = 2, rel_heights = c(1, .07))


## supplementary figure 8b -----------------------------------
sub_pca_tmp <- sub_pca_final %>%
  filter(!correct_method %in% "log") %>%
  filter(!correct_method %in% "normae") %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  filter(dataset %in% "simulated") %>%
  filter(scenario %in% "confounded") %>%
  mutate(title = sprintf("%s\n(%s-corrected)", quant_method, Hmisc::capitalize(correct_level))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate(label = sprintf("%s\n(SNR = %.2f)", correct_method, snr)) %>%
  mutate_at("batch", ~ dictLabelsBatch[.]) %>%
  mutate_at("sample", ~ factor(Hmisc::capitalize(.),
                               levels = c("D5", "D6", "F7", "M8",
                                          "Group1", "Group2", "Group3"))) %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq"))) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

p_box3 <- ggplot(sub_pca_tmp, aes(x = correct_level, y = PC1)) +
  geom_boxplot(aes(fill = sample)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 12, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(.2, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsSample) +
  scale_y_continuous(n.breaks = 5, name = "PC1") +
  facet_wrap(~ scenario + correct_level, scales = "free_x", nrow = 1)

p_box4 <- ggplot(sub_pca_tmp, aes(x = correct_level, y = PC2)) +
  geom_boxplot(aes(fill = sample)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 12, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(.2, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsSample) +
  scale_y_continuous(n.breaks = 5, name = "PC2") +
  facet_wrap(~ scenario + correct_level, scales = "free_x", nrow = 1)

supp8b <- plot_grid(p_box3, p_box4,nrow = 1)


## combined supplementary figure8 -----------------------------------
supp8 <- plot_grid(supp8a, supp8b,
                   nrow = 2, rel_heights = c(1, .3),
                   labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure8.pdf",
       supp8, width = 210, height = 290, units = "mm")


## figure 4c -----------------------------------
sub_pvca_final <- data.table::fread("./results/tables/pvca.csv")

## bar
sub_pvca1 <- sub_pvca_final %>%
  filter(dataset %in% "simulated") %>%
  filter(scenario %in% "confounded") %>%
  filter(data_level %in% "protein") %>%
  # filter(quant_method %in% "maxlfq") %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate(label2 = apply(., 1, function(x) {
    if (as.character(x[1] %in% "sample")) {
      return("Biological effects")
    } else if (as.character(x[1]) %in% "Residual") {
      "Residual"
    } else if (grepl(":sample", x[1])|grepl("sample:", x[1])) {
      "Biological:Technical interative effects"
    } else {
      "Technical effects"
    }
  })) %>%
  group_by(quant_method) %>%
  mutate(cut_off = Proportion[correct_method %in% "log" & label2 %in% "Biological effects"]) %>%
  ungroup %>%
  filter(!correct_method %in% "log") %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate(label = paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")) %>%
  mutate(name = paste(scenario, label, correct_method, quant_method, sep = "_")) %>%
  mutate(name = as.numeric(factor(name))) %>%
  group_by(name, cut_off, correct_level, correct_method, quant_method, label, label2) %>%
  dplyr::summarise(proportion = sum(Proportion)) %>%
  ungroup() %>%
  dplyr::arrange(desc(proportion), .by_group = TRUE)

sub_pvca1$name <- factor(sub_pvca1$name, levels = unique(sub_pvca1$name[sub_pvca1$label2 %in% "Biological effects"]))
sub_pvca1$label <- factor(sub_pvca1$label, levels = names(dictColorsLevel))
sub_pvca1$label2 <- factor(sub_pvca1$label2,
                           levels = c("Residual",
                                      "Technical effects",
                                      "Biological:Technical interative effects",
                                      "Biological effects"))

p_pvca_bar <- ggplot(sub_pvca1, aes(x = name, y = proportion)) +
  geom_col(aes(fill = label2), width = .7) +
  # geom_hline(aes(yintercept = cut_off), lty = 2, size = 1) +
  theme_classic() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, margin = unit(c(0, 0, 0, 0), "cm")),
        axis.text.y = element_text(size = 10)) +
  labs(y = "Weighted proportion of variances (%)", fill = "Label") +
  scale_fill_manual(values = c("#F7F7F7", "#FDE0EF", "#B8E186", "#66BD63")) +
  scale_y_continuous(breaks = seq(0, 1, .1),
                     labels = ~ . * 100)

sub_pvca2 <- sub_pvca1 %>%
  distinct(name, correct_method, quant_method, label)

p_pvca_level <- ggplot(sub_pvca2, aes(x = name, y = 1)) +
  geom_col(aes(fill = label), width = .7) +
  scale_fill_manual(values = dictColorsLevel, name = "Correction Level") +
  theme_void() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = unit(c(0, .2, 0, .35), "cm"))) +
  labs(y = "Level")

p_pvca_method <- ggplot(sub_pvca2, aes(x = name, y = 1)) +
  geom_col(aes(fill = correct_method), width = .7) +
  scale_fill_manual(values = dictColorsMethod, name = "BECA") +
  theme_void() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = unit(c(0, .2, 0, .22), "cm"))) +
  labs(y = "BECA")

p_pvca_quant <- ggplot(sub_pvca2, aes(x = name, y = 1)) +
  geom_col(aes(fill = quant_method), width = .7) +
  scale_fill_manual(values = dictColorsQuantMethods, name = "QM") +
  theme_void() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = unit(c(0, .2, 0, .61), "cm"))) +
  labs(y = "QM")

p_pvca_annot <- plot_grid(p_pvca_level + theme(legend.position = "none"),
                          p_pvca_method + theme(legend.position = "none"),
                          p_pvca_quant + theme(legend.position = "none"),
                          nrow = 3) +
  theme(plot.margin = unit(c(0, .2, .7, 0), "cm"))

p_pvca_final <- plot_grid(p_pvca_bar + theme(legend.position = "none"),
                          p_pvca_annot, nrow = 2,
                          rel_heights = c(5, 1))

p_pvca_legend <- plot_grid(get_legend(p_pvca_bar),
                           get_legend(p_pvca_level),
                           get_legend(p_pvca_method),
                           get_legend(p_pvca_quant),
                           ncol = 4, rel_widths = c(1, .8, .7, .7)) +
  theme(plot.margin = unit(c(0, 0, .2, 1), "cm")) 

figure4c <- plot_grid(p_pvca_final, p_pvca_legend,
                      nrow = 2, rel_heights = c(3.7, 1)) +
  theme(plot.margin = unit(c(0, 1.2, .5, .2), "cm"))


## combined figure4 -----------------------
figure4ab <- plot_grid(figure4a, figure4b,
                       ncol = 2, rel_widths = c(.6, 1),
                       labels = c("a", "b"), label_size = 20)
figure4 <- plot_grid(figure4ab, figure4c,
                     nrow = 2, rel_heights = c(1, 1.2),
                     labels = c("", "c"), label_size = 20)

ggsave("./results/figures/figure4.pdf",
       figure4, width = 210, height = 297, units = "mm")


## supplementary figure 9 -----------------------
sub_pvca_tmp <- sub_pvca_final %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  filter(!grepl("_log", correct_method)) %>%
  mutate(label2 = `Random Effect`) %>%
  group_by(dataset, scenario, correct_level, correct_method, quant_method, label2) %>%
  dplyr::summarise(proportion = sum(Proportion)) %>%
  ungroup() %>%
  group_by(quant_method, label2, scenario, dataset) %>%
  mutate(cut_off = proportion[correct_method %in% "log"]) %>%
  ungroup %>%
  filter(!correct_method %in% "log") %>%
  mutate(label = paste(Hmisc::capitalize(correct_level), "-corrected", sep = "")) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("quant_method", ~ ifelse(. %in% "", "ibaq_maxlfq_toppep3", .)) %>%
  tidyr::separate_rows(quant_method, sep = "_") %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq")))

p_list <- pblapply(unique(sub_pvca_tmp$scenario), function(scenario_id) {
  
  sub_pvca_tmp_i <- sub_pvca_tmp %>%
    filter(label2 %in% c("batch", "sample", "sample:mode", "batch:sample",
                         "mode", "lab", "lab:mode", "sample:lab")) %>%
    mutate_at("label2", ~ factor(., levels = c("mode", "lab", "batch", "lab:mode",
                                               "sample","sample:mode", "batch:sample",
                                               "sample:lab"))) %>%
    filter(scenario %in% scenario_id)  
  
  p_title <- ggplot(sub_pvca_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_rect(linewidth = .7),
          plot.margin = unit(c(.5, .5, .5, 1.75), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_pvca_box <- ggplot(sub_pvca_tmp_i, aes(x = label, y = proportion)) +
    geom_boxplot(aes(fill = label), width = .7,
                 position = position_dodge(width = .8)) +
    geom_hline(aes(yintercept = cut_off), lty = 2, color = "red") +
    # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
    # geom_signif(comparisons = list(c("Precursor-corrected", "Protein-corrected"),
    #                                c("Precursor-corrected", "Protein-corrected"),
    #                                c("Peptide-corrected", "Protein-corrected")),
    #             map_signif_level = TRUE,
    #             y_position = 1,
    #             step_increase = .06,
    #             tip_length = .01,
    #             textsize = 4,
    #             test = "t.test") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 12, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_rect(fill = "white"),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
    scale_fill_manual(values = dictColorsLevel) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                       labels = ~ .*100, name = "PVCA (%)") +
    facet_wrap( ~ label2, ncol = 3, shrink = FALSE, scales = "fixed")
  
  p_final <- plot_grid(p_title, p_pvca_box, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
})

supp9 <- plot_grid(plotlist = p_list,
                   nrow = 2, rel_heights = c(1, .7),
                   labels = c("a", "b", "c", "d"), label_size = 20)

ggsave("./results/figures/extended_figure9.pdf", supp9,
       width = 210, height = 240, units = "mm")


## supplementary table ------------------------
sub_gpca_final <- fread("./results/tables/gpca.csv")

sub_gpca_tmp <- sub_gpca_final %>%
  distinct(gpca_delta, dataset, scenario, data_level, correct_level, correct_method, quant_method) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  mutate(label = ifelse(correct_method %in% "log",
                        "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_gpca_stat <- sub_gpca_tmp %>%
  group_by(label, scenario) %>%
  summarise(gpca_mean = mean(gpca_delta),gpca_sd = sd(gpca_delta),
            gpca_median = median(gpca_delta))

p_gpca_box <- ggplot(sub_gpca_tmp, aes(x = label, y = gpca_delta)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsLevel, name = "Correction level") +
  scale_y_continuous(name = "gPCA delta") +
  facet_wrap( ~ scenario);p_gpca_box

fwrite(sub_gpca_stat, "./results/tables/extended_table5.csv")




