# ' Title: Figure3-quality assessment-CV/MCC.
# ' Author: Qiaochu Chen
# ‘ Date: Jun 30th, 2025

library(parallel)
library(pbapply)
library(readxl)
library(openxlsx)
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


## source data: Test vs Reference --------------------
mcc_objs <- readRDS("./results/tables/mcc.rds")
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) {
  df_tmp <- tmp_obj$df %>%
    distinct(class, value, estimate, dataset, scenario, data_level,
             correct_level, correct_method, quant_method, label)
})
sub_ref_final <- mcc_objs_tmp %>%
  rbindlist %>%
  filter(scenario %in% "confounded") %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  filter(correct_method %in% c("ruv", "combat", "ratio")) %>%
  rename(Reference = value, Test = estimate) %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")])%>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(sub_ref_final, "./results/tables/7-TestReference-source_data.xlsx")


## figure 3a ------------------------
cv_objs <- readRDS("./results/tables/cv.rds")
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  tmp_table <- tmp_table %>%
    select(any_of(c("precursor", "peptide", "protein")),
           cv, dataset, scenario, data_level,
           correct_level, correct_method, quant_method) %>%
    unite("feature", any_of(c("precursor", "peptide", "protein")),
          sep = "_", remove = TRUE, na.rm = TRUE)
})

sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  filter(!grepl("_log", correct_method)) %>%
  filter(data_level %in% "protein") %>%
  filter(quant_method %in% "maxlfq") %>%
  mutate(label = ifelse(correct_method %in% "log",
                        "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel)))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(correct_level, scenario, label, quant_method) %>%
  summarise(cv_median = median(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv))

p_cv_box <- ggplot(sub_cv_tmp, aes(x = label, y = cv)) +
  geom_boxplot(aes(fill = label), width = .7, outlier.size = .5,
               position = position_dodge(width = .8)) +
  geom_signif(comparisons = list(c("Precursor-corrected", "Peptide-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = p_format,
              vjust = -.5,
              step_increase = .3,
              tip_length = .03,
              textsize = 2,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 5, name = "CV",
                     expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap( ~ scenario, ncol = 4)

figure3a <- p_cv_box


## supplementary figure 5a ------------------------
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  tmp_table <- tmp_table %>%
    select(any_of(c("precursor", "peptide", "protein")),
           cv, dataset, scenario, data_level,
           correct_level, correct_method, quant_method) %>%
    unite("feature", any_of(c("precursor", "peptide", "protein")),
          sep = "_", remove = TRUE, na.rm = TRUE)
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
  mutate(label = ifelse(correct_method %in% "log",
                        "Uncorrected",
                        paste(Hmisc::capitalize(correct_level), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("label", ~ factor(., levels = names(dictColorsLevel))) %>%
  mutate(x_axis = as.numeric(label))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(correct_level, scenario, label, correct_method) %>%
  summarise(cv_median = median(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv))

p_cv_bar <- ggplot(sub_cv_tmp, aes(x = label, y = cv)) +
  geom_boxplot(aes(x = label, y = cv, fill = label), width = .7,
               position = position_dodge(width = .8)) +
  # stat_summary(aes(fill = label), fun = mean, geom = "bar", width = .7) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Uncorrected", "Protein-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = p_format,
              # y_position = 1,
              step_increase = .2,
              tip_length = .05,
              textsize = 1.8,
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
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .9), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 10, name = "CV",
                     expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap( ~ correct_method,
              labeller = as_labeller(dictLabelsCorrectMethods), 
              nrow = 1)

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
  summarise(mcc_mean = mean(mcc), mcc_sd = sd(mcc), n = length(mcc))

p_mcc_bar <- ggplot(sub_mcc_tmp, aes(x = label, y = mcc)) +
  geom_boxplot(aes(fill = label), width = .7,
               position = position_dodge(width = .8)) +
  geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
  geom_signif(comparisons = list(c("Precursor-corrected", "Peptide-corrected"),
                                 c("Precursor-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Protein-corrected")),
              map_signif_level = p_format,
              # y_position = .75,
              step_increase = .2,
              tip_length = .05,
              textsize = 2.5,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .2, .5, .35), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 8, name = "MCC",
                     expand = expansion(mult = c(0.05, 0.2))) +
  facet_grid(cols = vars(quant_method), rows = vars(scenario))

suppl5b <- p_mcc_bar


## combind supplementary figure5 -----------------------
supp5 <- plot_grid(suppl5a, suppl5b,
                   nrow = 2, rel_heights = c(.6, 1),
                   labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure5.pdf",
       supp5, width = 210, height = 260, units = "mm")


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
  geom_point(aes(color = label), size = 3, position = position_dodge(width = .9)) +
  geom_hline(aes(yintercept = cut_off), lty = 2, col = "red") +
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
  scale_color_manual(values = dictColorsLevel) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 5, name = "MCC",
                     expand = expansion(mult = c(0.05, 0.1))) +
  facet_grid(cols = vars(correct_method), rows = vars(scenario))

figure3b <- p_mcc_bar


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
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_color_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 5, name = "Test") +
  scale_x_continuous(limits = c(-3.6, 3.6), name = "Reference") +
  facet_grid(cols = vars(label2), rows = vars(correct_method))

figure3c <- p_tmp


## combind figure3 -----------------------
figure3 <- plot_grid(figure3a + theme(panel.spacing = unit(.5, "cm"),
                                      plot.margin = unit(c(0, 1.5, .5, .7), "cm")) +
                       labs(fill = "Correction level"),
                     figure3b + theme(plot.margin = unit(c(.5, .4, .5, .35), "cm")),
                     figure3c + theme(panel.spacing = unit(.5, "cm"),
                                      plot.margin = unit(c(.5, .4, .5, .1), "cm")),
                     nrow = 3, rel_heights = c(.8, 1.1, 2),
                     labels = c("a", "b", "c"), label_size = 20)

ggsave(paste("./results/figures/figure3.pdf", sep = ""),
       figure3, width = 210, height = 297, units = "mm")


