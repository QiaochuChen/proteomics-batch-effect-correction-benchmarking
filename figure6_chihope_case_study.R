# ' Title: Figure5-correction vs quantification
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
dictColorsSample <- c("#4D9221","#C51B7D", "black", "#999")
names(dictColorsSample) <- c("QC sample (P10)", "QC sample (P11)", "QC sample (PM)", "Study sample")

dictShapesSample <- c(15, 17, 16, 16)
names(dictShapesSample) <- c("QC sample (P10)", "QC sample (P11)", "QC sample (PM)", "Study sample")

dictColorsBatch <- c("#FB9A99", "#FDBF6F", "#1F78B4")
names(dictColorsBatch) <- c("Batch1", "Batch2", "Batch3")

dictColorsSex <- c("F" = "#e98989", "M" = "#78b7ee")
dictColorsWeek <- c("Baseline" = "#1B9E77", "Week24" = "#D95F02")

dictColorsTV <- c("#7db7e6", "#e6e384")
names(dictColorsTV) <- 1:2

dictClassNames <- c("Training", "Validation")
names(dictClassNames) <- 1:2

dictLabelsCorrectMethods <- c("LOESS", "Ratio", "Med-C", "Combat", "RUV-III-C", "Harmony")
names(dictLabelsCorrectMethods) <- c("loess", "ratio", "median", "combat", "ruv", "harmony")

dictLabelsCorrectMethods2 <- paste("LOESS +", dictLabelsCorrectMethods[2:6])
names(dictLabelsCorrectMethods2) <- paste(names(dictLabelsCorrectMethods)[2:6], "_loess", sep = "")

dictLabelsCorrectMethods <- c(dictLabelsCorrectMethods, dictLabelsCorrectMethods2)
dictLabelsCorrectMethods <- dictLabelsCorrectMethods[c(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11)]

dictLabelsCorrectMethods <- c(dictLabelsCorrectMethods, "negctrl" = "Negative Control")

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


## source data: MCC & R square --------------------
df_mcc_rsquare <- fread("./results/tables/mcc_rsquare_chihope.csv")
df_mcc_rsquare_revised <- df_mcc_rsquare %>%
  mutate_at("correct_method", ~ ifelse(grepl("log", .), "Uncorrected", dictLabelsCorrectMethods[.]))

write.xlsx(df_mcc_rsquare_revised, "./results/tables/7-ChiHOPE-source_data.xlsx")


## supplementary figure 11-------------------------
all_files <- list.files("./ChiHOPE", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("protein/maxlfq/expdata", all_files)]

expr_list <- pblapply(queried_files, function(file_id) {
  
  dataset <- "ChiHOPE"
  scenario <- "ChiHOPE"
  data_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", scenario))
  quant_method <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=\\/.)", data_level))
  correct_level <- str_extract(file_id, sprintf("(?<=.%s/).+?(?=_corrected.)", quant_method))
  correct_level <- ifelse(is.na(correct_level), data_level, correct_level)
  correct_method <- str_extract(file_id,"(?<=expdata_).+?(?=.csv)")
  correct_method <- gsub("_cutoff_NAs", "", correct_method)
  
  df_expr <- fread(file_id)
  df_meta <- fread("./ChiHOPE/meta.csv")
  
  sub_ms <- df_expr %>%
    reshape2::melt(., na.rm = TRUE, id = 1, variable.name = "run_id") %>%
    group_by(run_id) %>%
    summarise(intensity = median(value, na.rm = TRUE)) %>%
    left_join(., df_meta, by = "run_id") %>%
    mutate(dataset = dataset,
           scenario = scenario,
           data_level = data_level,
           correct_level = correct_level,
           correct_method = correct_method,
           quant_method = quant_method)
  
  return(sub_ms)
})
sub_ms_final <- rbindlist(expr_list)

sub_ms_tmp <- sub_ms_final %>%
  filter(!run_id %in% "ExpB19") %>%
  # filter(correct_method %in% c("log", "loess")) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.), "Uncorrected", paste(., "corrected"))) %>%
  mutate_at("correct_method", ~ factor(., levels = c("Uncorrected", paste(dictLabelsCorrectMethods, "corrected"))))

sub_ms_final <- sub_ms_tmp[, c(1:5, 27:32)]
write.xlsx(sub_ms_final, "./results/tables/8-ChiHOPE-intensity-source_data.xlsx")

x_breaks <- c(503, 969)

p_drift <- ggplot(sub_ms_tmp, aes(x = order, y = intensity)) +
  geom_point(aes(color = order), size = 1) +
  geom_smooth(aes(group = batch), color = "darkgrey") +
  geom_vline(xintercept = x_breaks, lty = 2, color = "grey") +
  scale_color_gradientn(colors = c("#DF65B0", "#FD8D3C", "#74C476", "#6A51A3"),
                        values = scales::rescale(c(0, x_breaks, 1495)), name = "Batch",
                        breaks = c(1, x_breaks, 1495),
                        labels = c(1, x_breaks, 1495),
                        guide = "colorbar") +
  scale_x_continuous(breaks = c(0, x_breaks, 1495), name = "Injection order") +
  scale_y_continuous(n.breaks = 5, name = "Intensity (log transformed)") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(size = 14, hjust = .5, vjust = 1),
        legend.text = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1, title.position = "top")) +
  facet_wrap(~ correct_method)

p_drift_qc <- ggplot(sub_ms_tmp, aes(x = order, y = intensity)) +
  geom_point(aes(color = sample, shape = sample), size = 1, alpha = .7) +
  geom_smooth(aes(group = batch), color = "#2171B5", fill = "#6BAED6") +
  geom_vline(xintercept = x_breaks, lty = 2, color = "red") +
  scale_size_manual(values = c(rep(3, 3), 2)) +
  scale_shape_manual(values = dictShapesSample) +
  scale_color_manual(values = dictColorsSample) +
  scale_x_continuous(breaks = c(0, x_breaks, 1495), name = "Injection order") +
  scale_y_continuous(n.breaks = 5, name = "Protein quantities (transformed)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.spacing = unit(.5, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1, title.position = "top")) +
  facet_wrap(~ correct_method, ncol = 2, scale = "free")

supp11 <- p_drift_qc

ggsave("./results/figures/extended_figure11.pdf", supp11,
       width = 210, height = 280, units = "mm")


## supplementary figure 12a ------------------------
pca_objs <- readRDS(paste("./results/tables/pca_snr_chihope.rds", sep = ""))
pcs_tables <- pblapply(pca_objs, function(pca_obj) {
  sub_tmp <- pca_obj$pcs_values
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

sub_pca_tmp <- sub_pca_final %>%
  filter(!library %in% "ExpJ1") %>%
  filter(data_level %in% "protein") %>%
  filter(is.na(quant_method)|quant_method %in% "maxlfq") %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.), "Uncorrected", paste(., "corrected"))) %>%
  mutate_at("correct_method", ~ factor(., levels = c("Uncorrected", paste(dictLabelsCorrectMethods, "corrected")))) %>%
  mutate(label = sprintf("%s\n(SNR = %.2f)", correct_method, snr)) %>%
  mutate_at("data_level", ~ factor(., levels = c("precursor", "peptide", "protein")))

write.xlsx(sub_pca_tmp, "./results/tables/9-ChiHOPE-pca-source_data.xlsx")

sub_tmp <- sub_pca_tmp %>% distinct(correct_method, label)
dictLabelsTmp <- setNames(sub_tmp$label, sub_tmp$correct_method)

p_tmp <- ggplot(sub_pca_tmp, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample, shape = batch), size = 1, alpha = .7) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 7),
        strip.background = element_blank()) +
  scale_color_manual(values = dictColorsSample) +
  labs(color = "Sample", shape = "Batch") +
  facet_wrap(~ correct_method, labeller = as_labeller(dictLabelsTmp),
             scales = "free", ncol = 2) +
  guides(shape = guide_legend(nrow = 1, title.position = "top"),
         color = "none")

supp12a <- p_tmp


## supplementary figure 12b ------------------------
df_meta <- fread("./ChiHOPE/meta.csv")
umap_objs <- readRDS("./results/tables/umap_chihope.rds")
sub_umap_final <- rbindlist(umap_objs)

sub_umap_tmp <- sub_umap_final %>%
  # filter(!library %in% c("ExpJ1", "ExpB19")) %>%
  # filter(correct_method %in% "ruv_loess") %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.), "Uncorrected", paste(., "corrected"))) %>%
  mutate_at("correct_method", ~ factor(., levels = c("Uncorrected", paste(dictLabelsCorrectMethods, "corrected"))))

sub_umap_final <- sub_umap_tmp[, c(1:6, 28:33)]
write.xlsx(sub_umap_final, "./results/tables/10-ChiHOPE-umap-source_data.xlsx")

p_tmp <- ggplot(sub_umap_tmp, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = sample, shape = batch), size = 1, alpha = .7) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 7),
        strip.background = element_blank()) +
  scale_color_manual(values = dictColorsSample, name = "Sample") +
  facet_wrap(~ correct_method, scales = "free", ncol = 2) +
  guides(shape = "none",
         color = guide_legend(nrow = 2, title.position = "top"))

supp12b <- p_tmp


## combined supplementary figure 12-------------------------
supp12 <- plot_grid(supp12a + theme(plot.margin = unit(c(.2, .5, .8, .5), "cm")),
                    supp12b + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                    ncol = 2,
                    labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure12.pdf", supp12,
       width = 210, height = 297, units = "mm")


## figure6a -------------------------
df_mcc_rsquare <- fread("./results/tables/mcc_rsquare_chihope.csv")

df_test <- df_mcc_rsquare %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.),
                                       "Uncorrected", .)) %>%
  mutate_at("correct_method", ~ ifelse(. %in% c("Negative Control", "Uncorrected"),
                                       ., paste(., "corrected"))) %>%
  arrange(mcc)

df_test1 <- df_test %>%
  filter(correct_method %in% c("Uncorrected", "Negative Control") | grepl("LOESS", correct_method))
  
dictColorsTmp <- brewer.pal(8, "BrBG")
names(dictColorsTmp) <- df_test1$correct_method

p_mcc_bar <- ggplot(df_test1, aes(x = reorder(correct_method, mcc), y = mcc)) +
  geom_col(aes(fill = correct_method), alpha = .7, width = .1) +
  geom_point(aes(color = correct_method), size = 3) +
  geom_text(aes(label = sprintf("%.2f", mcc)), size = 3, vjust = -1) +
  geom_hline(yintercept = .3, lty = 2, col = "red") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsTmp) +
  scale_fill_manual(values = dictColorsTmp) +
  scale_y_continuous(limits = c(0, .5), n.breaks =5, name = "MCC")

figure6a <- p_mcc_bar


## figure6b -------------------------
df_test2 <- df_mcc_rsquare %>%
  select(correct_method, FN, FP, TN, TP) %>%
  reshape2::melt(., id = 1, variable.name = "class") %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.),
                                       "Uncorrected", .)) %>%
  mutate_at("correct_method", ~ ifelse(. %in% c("Negative Control", "Uncorrected"),
                                       ., paste(., "corrected"))) %>%
  mutate_at("correct_method", ~ factor(., levels = df_test$correct_method))

df_test2 <- df_test2 %>%
  filter(correct_method %in% c("Uncorrected", "Negative Control") | grepl("LOESS", correct_method))

p_mcc_bar <- ggplot(df_test2, aes(x = correct_method, y = value)) +
  geom_col(aes(fill = correct_method), alpha = .7, width = .1) +
  geom_point(aes(color = correct_method), size = 3) +
  # geom_text(aes(label = value), size = 5, vjust = -1) +
  # geom_hline(yintercept = .3, lty = 2, col = "red") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_rect(fill = "white"),
        panel.spacing.x = unit(3.7, "lines"),
        panel.spacing.y = unit(1, "lines"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .35), "cm")) +
  scale_color_manual(values = dictColorsTmp) +
  scale_fill_manual(values = dictColorsTmp) +
  scale_y_continuous(n.breaks = 5, name = "Number of Samples",
                     expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~ class, scales = "free_y")

figure6b <- p_mcc_bar


## figure6c -------------------------
df_test3 <- df_mcc_rsquare %>%
  select(correct_method, r_square) %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.),
                                       "Uncorrected", .)) %>%
  mutate_at("correct_method", ~ ifelse(. %in% c("Negative Control", "Uncorrected"),
                                       ., paste(., "corrected"))) %>%
  mutate_at("correct_method", ~ factor(., levels = df_test$correct_method))

df_test3 <- df_test3 %>%
  filter(correct_method %in% c("Uncorrected", "Negative Control") | grepl("LOESS", correct_method))

p_r_bar <- ggplot(df_test3, aes(x = correct_method, y = r_square)) +
  geom_col(aes(fill = correct_method), alpha = .7, width = .1) +
  geom_point(aes(color = correct_method), size = 3) +
  geom_text(aes(label = sprintf("%.2f", r_square),
                vjust = ifelse(r_square > 0, -1, 2)), size = 3) +
  # geom_hline(yintercept = .3, lty = 2, col = "red") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsTmp) +
  scale_fill_manual(values = dictColorsTmp) +
  scale_y_continuous(limits = c(-.1, .3), n.breaks = 5, name = "R square")

figure6c <- p_r_bar


## figure6d -------------------------
all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
queried_files <- all_files[grepl("predictions", all_files)]
queried_files <- queried_files[grepl("1431", queried_files)]

predicted_tables <- pblapply(queried_files, function(file_id) {
  df_predicted <- read_excel(file_id)
  df_predicted_test <- df_predicted %>%
    filter(Training.Validation == 2) %>%
    select(correct_method, Age, prediction_age)
  
  if(grepl("negctrl", file_id)) df_predicted_test$correct_method <- "negctrl"
  
  return(df_predicted_test)
})

df_test4 <- predicted_tables %>%
  rbindlist %>%
  mutate_at("correct_method", ~ dictLabelsCorrectMethods[.]) %>%
  mutate_at("correct_method", ~ ifelse(is.na(.),
                                       "Uncorrected", .)) %>%
  mutate_at("correct_method", ~ ifelse(. %in% c("Negative Control", "Uncorrected"),
                                       ., paste(., "corrected"))) %>%
  filter(correct_method %in% df_test1$correct_method) %>%
  mutate_at("correct_method", ~ factor(., levels = df_test1$correct_method))

df_test4 <- df_test4 %>%
  filter(correct_method %in% c("Uncorrected", "Negative Control") | grepl("LOESS", correct_method))

write.xlsx(df_test4, "./results/tables/7-ChiHOPE-TestReference-source_data.xlsx")

p_tmp <- ggplot(df_test4, aes(x = Age, y = prediction_age)) +
  geom_point(aes(color = correct_method), size = .8, alpha = .5, shape = 16) +
  geom_smooth(aes(color = correct_method), method = "lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           p.digits = 2, size = 3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsTmp) +
  scale_y_continuous(limits = c(20, 75), n.breaks = 6, name = "Predicted Age") +
  scale_x_continuous(limits = c(20, 75), n.breaks = 6, name = "Age") +
  facet_wrap(~ correct_method, ncol = 4)

figure6d <- p_tmp


## combined figure 6-------------------------
figure6ac <- plot_grid(figure6a, figure6c,
                       ncol = 2, rel_widths = c(1, 1),
                       labels = c("a", "c"), label_size = 20)

figure6bd <- plot_grid(figure6b, figure6d,
                       nrow = 2, rel_heights = c(1.2, 1),
                       labels = c("b", "d"), label_size = 20)

figure6 <- plot_grid(figure6ac, figure6bd,
                     nrow = 2, rel_heights = c(.21, 1))

ggsave("./results/figures/figure6.pdf", figure6,
       width = 210, height = 297, units = "mm")


## supplementary figure 13a -------------------------
df_test5 <- df_test %>%
  filter(! correct_method %in% c("Uncorrected", "Negative Control", "LOESS corrected")) %>%
  mutate(group = ifelse(grepl("LOESS", correct_method), "LOESS-dependent", "LOESS-independent"))

p_mcc_box <- ggplot(df_test5, aes(x = group, y = mcc)) +
  geom_boxplot(aes(fill = group), width = .7, alpha = .8,
               position = position_dodge(width = .8)) +
  geom_signif(comparisons = list(c("LOESS-dependent", "LOESS-independent")),
              map_signif_level = p_format,
              vjust = -.5,
              y_position = .45,
              tip_length = .1,
              textsize = 6,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c("#1B7837", "#762A83")) +
  scale_y_continuous(limits = c(.3, .5), n.breaks = 5, name = "MCC")

sp_figure13a <- p_mcc_box


## supplementary figure 13b -------------------------
df_test5 <- df_test %>%
  filter(! correct_method %in% c("Uncorrected", "Negative Control", "LOESS corrected")) %>%
  mutate(group = ifelse(grepl("LOESS", correct_method), "LOESS-dependent", "LOESS-independent"))

sub_stat <- df_test5 %>%
  group_by(label, correct_method) %>%
  summarise(snr_mean = mean(snr),snr_sd = sd(snr), snr_median = median(snr),
            n_total = length(snr))

p_rsquare_box <- ggplot(df_test5, aes(x = group, y = r_square)) +
  geom_boxplot(aes(fill = group), width = .7, alpha = .8,
               position = position_dodge(width = .8)) +
  geom_signif(comparisons = list(c("LOESS-dependent", "LOESS-independent")),
              map_signif_level = p_format,
              vjust = -.5,
              y_position = .25,
              tip_length = .25,
              textsize = 6,
              test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c("#1B7837", "#762A83")) +
  scale_y_continuous(limits = c(.1, .3), n.breaks = 5, name = "R square")

sp_figure13b <- p_rsquare_box


## combined supplementary figure 13-------------------------
sp_figure13 <- plot_grid(sp_figure13a, sp_figure13b,
                       ncol = 2, rel_widths = c(1, 1),
                       labels = c("a", "b"), label_size = 20)

ggsave("./results/figures/extended_figure13.pdf", sp_figure13,
       width = 210, height = 120, units = "mm")









































