# ' Title: Figure2-Batch effect diagnosis.
# ' Author: Qiaochu Chen
# ‘ Date: Jun 30th, 2025

library(parallel)
library(pbapply)
library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(tidyr)
library(reshape2)
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)
library(cowplot)


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

dictShapesBatch <- c(12:17, 1:3)
names(dictShapesBatch) <- c("L1-1", "L2-1", "L3", "L1-2", "L2-2", "L4",
                            "Batch1", "Batch2", "Batch3")

dictLabelsBatch <- c("L1-1", "L2-1", "L3", "L1-2", "L2-2", "L4", "Batch1", "Batch2", "Batch3")
names(dictLabelsBatch) <- c("DDA_APT", "DDA_FDU", "DDA_NVG", "DIA_APT", "DIA_BGI", "DIA_FDU",
                            "batch1", "batch2", "batch3")

dictLabelsCorrectMethods <- c("Ratio", "Med-C", "Combat", "RUV-III-C", "Harmony", "WaveICA2", "NormAE")
names(dictLabelsCorrectMethods) <- c("ratio", "median", "combat", "ruv", "harmony", "waveica", "normae")


## source data: CV --------------------
cv_objs <- readRDS("./results/tables/cv.rds")
cv_objs_tmp <- pblapply(cv_objs, function(tmp_table) {
  tmp_table <- tmp_table %>%
    select(any_of(c("precursor", "peptide", "protein")),
           cv, dataset, scenario, data_level,
           correct_level, correct_method, quant_method) %>%
    unite("feature", any_of(c("precursor", "peptide", "protein")),
          sep = "_", remove = TRUE, na.rm = TRUE)
})

sub_cv_stat <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(dataset, data_level, scenario,
           correct_level, correct_method, quant_method) %>%
  summarise(cv_median = median(cv), cv_mean = mean(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv), .groups = "drop") %>%
  ungroup %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(sub_cv_stat, "./results/tables/1-CV-source_data.xlsx")


## source data: MCC --------------------
mcc_objs <- readRDS("./results/tables/mcc.rds")
mcc_objs_tmp <- pblapply(mcc_objs, function(tmp_obj) {
  tmp_table <- tmp_obj$mcc
})
mcc_df <- rbindlist(mcc_objs_tmp)
mcc_df_final <- mcc_df[, 1:11] %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(mcc_df_final, "./results/tables/2-MCC-source_data.xlsx")


## source data: PCA --------------------
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
sub_pca_final_revised <- sub_pca_final %>%
  rename(run_id = library) %>%
  mutate_at("batch", ~ dictLabelsBatch[.]) %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(sub_pca_final_revised, "./results/tables/3-PCA-source_data.xlsx")


## source data: SNR --------------------
pca_objs <- readRDS(paste("./results/tables/pca_snr.rds", sep = ""))
snr_tables <- pblapply(pca_objs, function(pca_obj) {
  sub_snr <- data.frame(dataset = pca_obj$dataset,
                        scenario = pca_obj$scenario,
                        data_level = pca_obj$data_level,
                        correct_level = pca_obj$correct_level,
                        correct_method = pca_obj$correct_method,
                        quant_method = pca_obj$quant_method,
                        snr = pca_obj$snr_results$snr)
  
  return(sub_snr)
})
sub_snr_final <- rbindlist(snr_tables)
sub_snr_final_revised <- sub_snr_final %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(sub_snr_final_revised, "./results/tables/4-SNR-source_data.xlsx")


## source data: PVCA --------------------
sub_pvca_final <- fread(paste("./results/tables/pvca.csv", sep = ""))
sub_pvca_final_revised <- sub_pvca_final %>%
  mutate_at("data_level", ~ paste(Hmisc::capitalize(.), "-level", sep = "")) %>%
  mutate_at("correct_level", ~ ifelse(correct_method %in% "log",  "Uncorrected",
                                      paste(Hmisc::capitalize(.), "-corrected", sep = ""))) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")])%>%
  mutate_at("correct_method", ~ ifelse(grepl("_log", .), ., dictLabelsCorrectMethods[.]))

write.xlsx(sub_pvca_final_revised, "./results/tables/5-PVCA-source_data.xlsx")


## supplementary figure 1a --------------------
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  # filter(cv < 1) %>%
  filter(correct_method %in% "log") %>%
  filter(dataset %in% "quartet") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(data_level, scenario) %>%
  summarise(cv_mean = mean(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv))

p_cv_box <- ggplot(sub_cv_tmp, aes(x = data_level, y = cv)) +
  geom_boxplot(aes(fill = data_level), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(n.breaks = 5, name = "CV") +
  scale_x_discrete(labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  facet_wrap(~ scenario, ncol = 4)

supp1a <- p_cv_box


## supplementary figure 1b -------------------
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  # filter(cv < 1) %>%
  filter(dataset %in% c("simulated")) %>%
  filter(grepl("_log", correct_method)) %>%
  mutate_at("dataset", Hmisc::capitalize) %>%
  mutate_at("correct_method",
            ~ factor(., levels = c("gamma_log", "withbiol_log", "withscal_log", "withbatch_log")))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(data_level, scenario, correct_method) %>%
  summarise(cv_mean = mean(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv))

p_cv_box <- ggplot(sub_cv_tmp, aes(x = correct_method, y = cv)) +
  geom_boxplot(aes(fill = correct_method), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  scale_y_continuous(n.breaks = 10, name = "CV") +
  scale_x_discrete(labels = ~ gsub("_log", "", Hmisc::capitalize(.))) +
  facet_wrap(~ dataset, nrow = 1)

supp1b <- p_cv_box


## supplementary figure 1c -----------------------
sub_cv_tmp <- cv_objs_tmp %>%
  rbindlist %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  # filter(cv < 1) %>%
  filter(dataset %in% c("simulated")) %>%
  filter(correct_method %in% "log") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

sub_cv_stat <- sub_cv_tmp %>%
  filter(!(is.na(cv) | is.infinite(cv))) %>%
  group_by(data_level, scenario) %>%
  summarise(cv_mean = mean(cv), cv_sd = sd(cv),
            feature_n = length(unique(feature)),
            total_n = length(cv))

p_cv_box <- ggplot(sub_cv_tmp, aes(x = data_level, y = cv)) +
  geom_boxplot(aes(fill = data_level), width = .7,
               position = position_dodge(width = .8)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_brewer(palette = "Blues", direction = 1) +
  scale_y_continuous(n.breaks = 10, name = "CV") +
  scale_x_discrete(labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  facet_wrap(~ scenario)

supp1c <- p_cv_box


## supplementary figure 1d --------------------
sub_mcc_tmp <- mcc_df %>%
  filter(grepl("_log", correct_method)) %>%
  distinct(dataset, scenario, data_level, correct_level, correct_method, quant_method, mcc) %>%
  mutate_at("dataset", Hmisc::capitalize) %>%
  mutate_at("correct_method",
            ~ factor(., levels = c("gamma_log", "withbiol_log", "withscal_log", "withbatch_log")))

p_mcc_bar <- ggplot(sub_mcc_tmp, aes(x = correct_method, y = mcc)) +
  geom_col(aes(fill = correct_method), width = .06, position = position_dodge(width = .9)) +
  geom_point(aes(color = correct_method), size = 3, position = position_dodge(width = .9)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        panel.spacing = unit(.8, "cm"),
        plot.margin = unit(c(0, 0, .5, .5), "cm")) +
  scale_fill_brewer(palette = "Greys", name = "Data level",
                    labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_color_brewer(palette = "Greys", name = "Data level",
                     labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_y_continuous(name = "MCC", n.breaks = 5) +
  facet_wrap(~ dataset, scales = "fixed");p_mcc_bar

supp1d <- p_mcc_bar


## supplementary figure 1e --------------------
sub_mcc_tmp <- mcc_df %>%
  filter(correct_method %in% "log") %>%
  distinct(dataset, scenario, data_level, correct_level, correct_method, quant_method, mcc) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ ifelse(is.na(.) |. %in% "", "ibaq_maxlfq_toppep3", .)) %>%
  tidyr::separate_rows(quant_method, sep = "_") %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq")))

p_mcc_bar <- ggplot(sub_mcc_tmp, aes(x = quant_method, y = mcc)) +
  geom_col(aes(fill = data_level), width = .1, position = position_dodge(width = .9)) +
  geom_point(aes(color = data_level), size = 3, position = position_dodge(width = .9)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        panel.spacing = unit(.8, "cm"),
        plot.margin = unit(c(0, 0, .5, .5), "cm")) +
  scale_fill_brewer(palette = "Blues", name = "Data level",
                    labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_color_brewer(palette = "Blues", name = "Data level",
                     labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_y_continuous(limits = c(-.1, 1), name = "MCC", n.breaks = 5) +
  facet_wrap(~ scenario, scales = "fixed")

supp1e <- p_mcc_bar


## combined supplementary figure1 -----------------------
supp1bc <- plot_grid(supp1b + theme(plot.margin = unit(c(.5, .5, 1.3, .5), "cm")),
                     supp1c, ncol = 2, rel_widths = c(.55, 1),
                     labels = c("b", "c"), label_size = 20)

supp1de <- plot_grid(supp1d + theme(plot.margin = unit(c(0, .5, 0, .5), "cm")),
                     supp1e + theme(plot.margin = unit(c(0, 0, 2, .5), "cm")),
                     ncol = 2, rel_widths = c(.5, 1),
                     labels = c("d", "e"), label_size = 20)

supp1 <- plot_grid(supp1a + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                   supp1bc, supp1de + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                   nrow = 3, labels = c("a", ""), label_size = 20)

ggsave("./results/figures/extended_figure1.pdf", supp1,
       width = 210, height = 270, units = "mm")


## supplementary figure 2a ----------------
sub_pca_tmp <- sub_pca_final %>%
  filter(grepl("_log", correct_method)) %>%
  filter(is.na(quant_method)|quant_method %in% "maxlfq") %>%
  mutate(label = sprintf("%s\n(SNR = %.2f)",
                         Hmisc::capitalize(gsub("_log", "", correct_method)), snr)) %>%
  mutate_at("dataset", Hmisc::capitalize) %>%
  mutate_at("scenario", Hmisc::capitalize) %>%
  mutate_at("batch", ~ dictLabelsBatch[.]) %>%
  mutate_at("sample", ~ factor(Hmisc::capitalize(.),
                               levels = c("D5", "D6", "F7", "M8",
                                          "Group1", "Group2", "Group3"))) %>%
  mutate_at("correct_method", ~ factor(., levels = c("withbiol_log", "withscal_log",
                                                     "withbatch_log")))

sub_tmp <- sub_pca_tmp %>% distinct(correct_method, label)
dictLabelsTmp <- setNames(sub_tmp$label, sub_tmp$correct_method)

p_title <- ggplot(sub_pca_tmp) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 16),
        plot.margin = unit(c(.5, .5, 0, 1.8), "cm"))+
  facet_grid(cols = vars(dataset))

p_tmp <- ggplot(sub_pca_tmp, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample, shape = batch), size = 3) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        plot.margin = unit(c(0, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsSample) +
  scale_shape_manual(values = dictShapesBatch) +
  labs(color = "Sample", shape = "Batch") +
  facet_wrap(~ correct_method, labeller = as_labeller(dictLabelsTmp),
             scales = "free") +
  guides(shape = guide_legend(nrow = 1, title.position = "top"),
         color = guide_legend(nrow = 1, title.position = "top"))

supp2a <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.14, 1))


## supplementary figure 2b -----------------------
sub_pvca_tmp <- sub_pvca_final %>%
  filter(grepl("_log", correct_method)) %>%
  # filter(`Random Effect` %in% "batch_id:group_id") %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

sub_pvca_tmp_i <- sub_pvca_tmp %>%
  filter(`Random Effect` %in% c("batch", "sample", "batch:sample")) %>%
  mutate_at("Random Effect", ~ factor(., levels = c("batch", "sample", "batch:sample")))

p_pvca_box <- ggplot(sub_pvca_tmp_i, aes(x = correct_method, y = Proportion)) +
  geom_col(aes(fill = Label), width = .1, position = position_dodge(width = .9)) +
  geom_point(aes(color = Label), size = 5, position = position_dodge(width = .9)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 16, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                     labels = ~ . * 100,
                     name = "Weighted average\nproportion of variance (%)") +
  scale_x_discrete(labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  facet_wrap(~ `Random Effect`, ncol = 3, shrink = FALSE, scales = "fixed")

supp2b <- p_pvca_box


## combined supplementary figure2 -----------------------
supp2 <- plot_grid(supp2a, supp2b,
                   nrow = 2, rel_heights = c(1, .9),
                   labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure2.pdf",
       supp2, width = 210, height = 230, units = "mm")


## figure 2a ------------------------
sub_pca_tmp <- sub_pca_final %>%
  filter(correct_method %in% "log") %>%
  filter(is.na(quant_method)|quant_method %in% "maxlfq") %>%
  mutate(label = sprintf("%s-level\n(SNR = %.2f)", Hmisc::capitalize(data_level), snr)) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("batch", ~ dictLabelsBatch[.]) %>%
  mutate_at("sample", ~ factor(Hmisc::capitalize(.),
                               levels = c("D5", "D6", "F7", "M8",
                                          "Group1", "Group2", "Group3"))) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

p_list <- pblapply(unique(sub_pca_tmp$scenario), function(scenario_id) {
  
  sub_pca_tmp_j <- sub_pca_tmp %>%
    filter(scenario %in% scenario_id)
  
  sub_tmp <- sub_pca_tmp_j %>% distinct(data_level, label)
  dictLabelsTmp <- setNames(sub_tmp$label, sub_tmp$data_level)
  
  p_title_j <- ggplot(sub_pca_tmp_j) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 14),
          plot.margin = unit(c(.5, .5, 0, 1.66), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_tmp <- ggplot(sub_pca_tmp_j, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, shape = batch), size = 3) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          title = element_text(size = 14),
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          plot.margin = unit(c(0, .5, 0, .5), "cm")) +
    scale_color_manual(values = dictColorsSample) +
    scale_shape_manual(values = dictShapesBatch) +
    labs(color = "Sample", shape = "Batch") +
    facet_wrap(~ data_level, labeller = as_labeller(dictLabelsTmp),
               scales = "free") +
    guides(shape = guide_legend(nrow = 1, title.position = "top"),
           color = guide_legend(nrow = 1, title.position = "top"))
  
  if (scenario_id %in% "Simulated-B") p_tmp <- p_tmp + theme(legend.position = "none")
  
  p_j <- plot_grid(p_title_j, p_tmp, nrow = 2, rel_heights = c(.2, 1))
  
  return(p_j)
})

figure2a <- p_list[[1]]


## supplementary figure 3 -----------
supp3 <- plot_grid(plotlist = p_list[2:4],
                   nrow = 3, rel_heights = c(1, .8, 1),
                   labels = c("a", "b", "c"), label_size = 20)

ggsave("./results/figures/extended_figure3.pdf",
       supp3, width = 210, height = 290, units = "mm")


## figure 2b ------------------------
sub_snr_tmp <- sub_pca_final %>%
  filter(correct_method %in% "log") %>%
  distinct(dataset, scenario, data_level, correct_level, correct_method, quant_method, snr) %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein"))) %>%
  mutate_at("quant_method", ~ ifelse(is.na(.), "ibaq_maxlfq_toppep3", .)) %>%
  tidyr::separate_rows(quant_method, sep = "_") %>%
  mutate_at("quant_method", ~ factor(., levels = c("ibaq", "toppep3", "maxlfq")))

p_snr_bar <- ggplot(sub_snr_tmp, aes(x = quant_method, y = snr)) +
  geom_col(aes(fill = data_level), width = .1, position = position_dodge(width = .9)) +
  geom_point(aes(color = data_level), size = 3, position = position_dodge(width = .9)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.spacing = unit(.8, "cm"),
        plot.margin = unit(c(0, 0, 0, .5), "cm")) +
  scale_fill_brewer(palette = "Blues", name = "Data level",
                    labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_color_brewer(palette = "Blues", name = "Data level",
                     labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
  scale_y_continuous(limits = c(-1, 10), name = "SNR", n.breaks = 5) +
  facet_wrap(~ scenario, ncol = 4, scales = "fixed") +
  guides(fill = guide_legend(nrow = 1, title.position = "top"))

figure2b <- p_snr_bar


## figure 2c -----------------------
sub_pvca_tmp <- sub_pvca_final %>%
  filter(correct_method %in% "log") %>%
  # filter(`Random Effect` %in% "batch_id:group_id") %>%
  mutate_at("scenario", ~ dictLabelsScenario[paste(dataset, ., sep = "_")]) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

p_list <- pblapply(unique(sub_pvca_tmp$scenario), function(scenario_id) {
  
  sub_pvca_tmp_i <- sub_pvca_tmp %>%
    filter(scenario %in% scenario_id) %>%
    filter(`Random Effect` %in% c("batch", "sample", "sample:mode", "batch:sample",
                                  "mode", "lab", "lab:mode", "sample:lab")) %>%
    mutate_at("Random Effect", ~ factor(., levels = c("mode", "lab", "batch", "sample",
                                                      "lab:mode", "sample:mode",
                                                      "batch:sample", "sample:lab")))
  
  p_title <- ggplot(sub_pvca_tmp_i) +
    theme_classic() +
    theme(line = element_blank(),
          strip.text = element_text(size = 14),
          plot.margin = unit(c(.5, .5, .5, 1.95), "cm"))+
    facet_grid(cols = vars(scenario))
  
  p_pvca_box <- ggplot(sub_pvca_tmp_i, aes(x = data_level, y = Proportion)) +
    stat_summary(aes(fill = data_level), fun = mean, geom = "bar", width = .7) +
    geom_point(aes(fill = data_level), shape = 21, color = "black", size = 1.5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 11, margin = unit(rep(.3, 4), "cm")),
          strip.background = element_rect(fill = "white"),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(0, .5, .5, .0), "cm")) +
    scale_fill_brewer(palette = "Blues", name = "Data level",
                      labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                       labels = ~ . * 100,
                       name = "Weighted average\nproportion of variance (%)") +
    scale_x_discrete(labels = ~ paste(Hmisc::capitalize(.), "-level", sep = "")) +
    facet_wrap(~ `Random Effect`, ncol = 3, shrink = FALSE, scales = "fixed") +
    guides(fill = guide_legend(nrow = 1, title.position = "top"))
  
  p_final <- plot_grid(p_title, p_pvca_box, nrow = 2, rel_heights = c(.1, 1));p_final
  
  return(p_final)
})

figure2c <- plot_grid(p_list[[1]], p_list[[3]], ncol = 2)


## combined figure2 -----------------------
figure2 <- plot_grid(figure2a + theme(plot.margin = unit(c(0, 0, 0, .5), "cm")),
                     figure2b + theme(plot.margin = unit(c(.5, .5, .5, .5), "cm")),
                     figure2c + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                     nrow = 3, rel_heights = c(1.1, .9, 1),
                     label_y = c(1, 1, 1.1),
                     labels = c("a", "b", "c"), label_size = 20)

ggsave(paste("./results/figures/figure2.pdf", sep = ""),
       figure2, width = 210, height = 297, units = "mm")


## supplementary figure 4 ---------------------
supp4 <- plot_grid(p_list[[2]], p_list[[4]],
                   nrow = 2, rel_heights = c(1, .6),
                   labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure4.pdf",
       supp4, width = 210, height = 240, units = "mm")

