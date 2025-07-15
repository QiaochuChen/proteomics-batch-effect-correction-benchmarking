# Author: Qiaochu Chen
# Date: Feb 4, 2023

library(dplyr)
library(umap)
library(ggplot2)
# library(ggbiplot)
library(ggExtra)
library(lme4)
library(parallel)

calculate_pca <- function(exprdata_t, metadata, center, scale, group = NULL) {

  meta <- metadata %>% tibble::rownames_to_column("library")

  pcs_promp_object <- prcomp(exprdata_t, retx = TRUE, center = center, scale. = scale)
  pca_results <- pcs_promp_object %>% summary
  pcs_props <- pca_results$importance %>% as.data.frame
  pcs_values <-  pca_results$x

  pcs_values <- pcs_values %>%
    as.data.frame %>%
    tibble::rownames_to_column("library") %>%
    right_join(meta, ., by = "library")

  if (!is.null(group)) {
    pcs_values <- pcs_values %>% dplyr::rename(sample = sym(group))
  }

  results <- list(
    pcs_promp_object = pcs_promp_object,
    pcs_values = pcs_values,
    pcs_props = pcs_props
  )

  return(results)
}

calculate_snr <- function(pcs_values, pcs_props, pct_threshold) {

  if(!is.null(pct_threshold)) {
    pcs_cut <- colnames(pcs_props)[which(pcs_props[3, ] > pct_threshold)][1]
    pcs_num <- match(pcs_cut, colnames(pcs_props))
  } else {
    # # 目前不考虑方差贡献率是否大于50%，仅选择前两个主成分
    pcs_num <- 2
  }

  pcs_props <- pcs_props[, 1:pcs_num]
  pcs_values <- pcs_values %>% 
    dplyr::select(any_of(c("library", "sample", colnames(pcs_props))))

  libs <- pcs_values$library
  sams <- pcs_values$sample

  pair_dt <- data.frame(lib1 = rep(libs, each = length(libs))) %>%
    mutate(lib2 = rep(libs, time = length(libs))) %>%
    mutate(grp1 = sams[match(lib1, libs)]) %>%
    mutate(grp2 = sams[match(lib2, libs)]) %>%
    mutate(type = apply(., 1, function(x) ifelse(x[3] == x[4], "Intra", "Inter"))) %>%
    filter(lib1 != lib2)

  dist_dt <- pair_dt %>%
    mutate(distance = apply(., 1, function(a) {

      pcs <- pcs_values %>%
        filter(library %in% a[1:2]) %>%
        dplyr::select(any_of(colnames(pcs_props)))

      dis <- apply(pcs, 2, function(b) (b[1] - b[2]) ^ 2)

      dis_weighted <-  rbind(dis, pcs_props[2, ]) %>%
        apply(., 2, function(c) c[1] * c[2]) %>%
        sum

      return(dis_weighted)

    })) %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(distance_mean = mean(distance, na.rm = TRUE))

  snr_results <- dist_dt %>%
    tibble::column_to_rownames("type") %>%
    t %>%
    as.data.frame %>%
    mutate(snr = 10 * log10(Inter / Intra))

  return(snr_results)
}

plot_pca <- function(pcs_values, pcs_props, dictColors, group, plot_title, snr_value = NULL) {

  axis_x <- sprintf("PC1 (%.2f%%)", pcs_props$PC1[2] * 100)
  axis_y <- sprintf("PC2 (%.2f%%)", pcs_props$PC2[2] * 100)

  if (!is.null(snr_value)) {
    subtitle <- sprintf("SNR = %.3f", snr_value)
  } else {
    subtitle <- ""
  }

  p <- ggplot(pcs_values, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 2) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = dictColors) +
    labs(subtitle = subtitle) +
    labs(x = axis_x, y = axis_y, colour = group, title = plot_title)

  # # 目前不考虑增加边际密度分布
  # p <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)

  return(p)
}

plot_pca_biplot <- function(pcs_promp_object, pcs_values, pcs_props, dictColors, dictGroups, group, plot_title = 'Principal Components Analysis (PCA) - Biplot', snr_value = NULL) {

  axis_x <- sprintf("PC1 (%.2f%%)", pcs_props$PC1[2] * 100)
  axis_y <- sprintf("PC2 (%.2f%%)", pcs_props$PC2[2] * 100)

  range_x <- range(pcs_values$PC1)
  range_y <- range(pcs_values$PC2)

  ratio <- (range_x[2] - range_x[1]) / (range_y[2] - range_y[1])

  if (!is.null(snr_value)) {
    subtitle <- sprintf("SNR = %.3f", snr_value)
  } else {
    subtitle <- ""
  }

  if(ncol(pcs_promp_object$x)>=5){
    with_vars <- F
    print('Sorry, there are TOO MANY variables which will not be labeled or drawn with arrows.')
  }else {
    with_vars <- T
  }

  p <- ggbiplot(pcs_promp_object, obs.scale = 1, groups = dictGroups, var.axes = with_vars, var.scale = 1, ellipse = T, circle = T, circle.prob = 0.95) +
    theme_bw() +
    coord_fixed(ratio = ratio) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") +
    scale_color_manual(values = dictColors) +
    labs(subtitle = subtitle) +
    labs(x = axis_x, y = axis_y, color = group, title = plot_title)

  return(p)
}

main_pca <- function(exprdata_t, metadata, dictColors, group, dictGroups = NULL, center = FALSE, scale = FALSE, pct_threshold = NULL, snr = FALSE, plot = TRUE, plot_title = "PCA", biplot = FALSE) {

  pcs_list <- calculate_pca(exprdata_t, metadata, center, scale, group)
  pcs_values <- pcs_list$pcs_values
  pcs_props <- pcs_list$pcs_props
  pcs_promp_object <- pcs_list$pcs_promp_object

  if (snr) {
    snr_results <- calculate_snr(pcs_values, pcs_props, pct_threshold)
    snr_value <- snr_results$snr
  } else {
    snr_results <- NULL
    snr_value <- NULL
  }

  if (biplot) {
    pca_plot <- plot_pca_biplot(pcs_promp_object, pcs_values, pcs_props, dictColors, dictGroups, group, plot_title, snr_value)
  } else if (plot) {
    pca_plot <- plot_pca(pcs_values, pcs_props, dictColors, group, plot_title, snr_value)
  } else {
    pca_plot <- NULL
  }

  result_list <- list(
    p = pca_plot,
    pcs_values = pcs_values,
    pcs_props = pcs_props,
    snr_results = snr_results
  )

  return(result_list)
}

calculate_umap <- function(exprdata_t, metadata, group, center, scale, pct_threshold = .6){

  meta <- metadata %>% tibble::rownames_to_column("library")

  ## PCA
  pcs_list <- calculate_pca(exprdata_t, metadata, center = center, scale = scale)
  pcs_values <- pcs_list$pcs_values
  pcs_props <- pcs_list$pcs_props

  ## 此时考虑方差贡献率
  pcs_cut <- colnames(pcs_props)[which(pcs_props[3, ] > pct_threshold)][1]
  pcs_num <- match(pcs_cut, colnames(pcs_props))

  pcs_props <- pcs_props[, 1:pcs_num]
  pcs_values <- pcs_values %>%
    tibble::column_to_rownames("library") %>%
    dplyr::select(any_of(colnames(pcs_props)))

  ## UMAP
  umap_settings <- umap.defaults
  umap_results <- umap(pcs_values, config = umap_settings)
  umap_layouts <- umap_results$layout %>% as.data.frame

  umap_final <- umap_layouts %>%
    tibble::rownames_to_column("library") %>%
    left_join(., meta, by = "library") %>%
    dplyr::rename(sample = sym(group))

  return(umap_final)
}

plot_umap <- function(umap_final, group, plot_title) {

  p <- ggplot(data = umap_final, aes(x = V1, y = V2))+
    geom_point(aes(color = sample), alpha=0.8,size=4) +
    labs(x = 'UMAP 1',y = 'UMAP 2') +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    theme_light() +
    theme(legend.position = "bottom") +
    labs(colour = group, title = plot_title)

  p_final <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)

  return(p_final)
}

main_umap <- function(exprdata_t, metadata, group, center = TRUE, scale = TRUE, pct_threshold = .6, plot = TRUE, plot_title = "UMAP") {

  umap_final <- calculate_umap(exprdata_t, metadata, group, center = center, scale = scale, pct_threshold = pct_threshold)

  if (plot) {
    p_final <- plot_umap(umap_final, group, plot_title)
  } else {
    p_final <- NULL
  }

  result_list <- list(
    p = p_final,
    table = umap_final
  )

  return(result_list)
}

calculate_pvca <- function(exprdata_t, metadata, pct_threshold, center, scale) {

  ## PCA
  pcs_list <- calculate_pca(exprdata_t, metadata, center = center, scale = scale)
  pcs_values <- pcs_list$pcs_values
  pcs_props <- pcs_list$pcs_props

  ## 此时考虑方差贡献率或只选择前三个主成分
  if(!is.null(pct_threshold)) {
    pcs_cut <- colnames(pcs_props)[which(pcs_props[3, ] > pct_threshold)][1]
    pcs_num <- match(pcs_cut, colnames(pcs_props))
  } else {
    pcs_num <- 3
  }

  pcs_props <- pcs_props[, 1:pcs_num]
  pcs_values <- pcs_values %>%
    dplyr::select(any_of(c(colnames(metadata), colnames(pcs_props))))

  ## 独立 & 交互影响因素
  fact_indep <- colnames(metadata)
  if (length(fact_indep) >= 2) {
    fact_inter <- combn(fact_indep, 2) %>% apply(., 2, function(x) paste0(x, collapse = ":"))
    fact_final <- c(fact_indep, fact_inter)
    x_formula <- paste("(1|", paste(fact_final,collapse = ") + (1|"), ")", sep = "")
  } else {
    x_formula <- paste("(1|", fact_indep, ")", sep = "")
  }

  ## 混合线性模型估计随机效应
  pcs_values[is.na(pcs_values)] <- ""
  pcs_data <- pcs_values %>% mutate_if(is.character, function(x) as.numeric(as.factor(x)))

  pcs_lmer <- mclapply(colnames(pcs_props), function(pc) {
    lmer_formula <- formula(paste(pc, x_formula, sep = " ~ "))
    lmer_results <- lmer(lmer_formula, pcs_data, REML = TRUE, verbose = FALSE)
    pcs_lmer <- lmer_results %>%
      VarCorr %>%
      as.data.frame %>%
      mutate(PC = pc)

    return(pcs_lmer)
  },
  mc.cores = 3)

  ## 计算随机效应对主成分贡献率的加权平均值
  pcs_lmer_data <- pcs_lmer %>%
    data.table::rbindlist(.) %>%
    reshape2::dcast(., grp ~ PC, value.var = "vcov") %>%
    tibble::column_to_rownames("grp") %>%
    apply(., 2, function(x) x / sum(x))

  lmer_weight <- pcs_lmer_data %>%
    apply(., 1, function(a) {
      lmer_weighted <- rbind(a, pcs_props) %>%
        apply(., 2, function(b) b[1] * b[3])
    })

  lmer_final <- lmer_weight %>%
    apply(., 2, function(x) sum(x) / sum(lmer_weight)) %>%
    data.frame(Proportion = .) %>%
    tibble::rownames_to_column("Random Effect") %>%
    arrange(desc(Proportion))

  orderFact <- c(setdiff(lmer_final$`Random Effect`, "Residual"), "Residual")

  sub_pvca <- lmer_final %>%
    mutate(Label = sapply(`Random Effect`, function(x) {
      if (grepl(":", x)) {
        return("Interactive")
      } else if (grepl("Residual", x)) {
        return("Residual")
      } else return("Independent")
    })) %>%
    mutate(`Random Effect` = factor(`Random Effect`, levels = orderFact))

  return(sub_pvca)

}

plot_pvca <- function(sub_pvca, plot_title) {

  p <- ggplot(sub_pvca, aes(x = `Random Effect`, y = Proportion)) +
    geom_bar(aes(fill = Label), stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(vjust = 1, hjust = 1, size = 16, angle = 45),
          legend.title = element_text(size = 12),
          legend.position = "right") +
    labs(y = "weighted average\nproportion of variance ",
         title = plot_title, fill = "Factor Type")

  return(p)

}

main_pvca <- function(exprdata_t, metadata, pct_threshold = NULL, center = TRUE, scale = TRUE, plot = TRUE, plot_title = "Principal Variance Compoment Analysis (PVCA)", show_legend = TRUE) {

  pvca_result <- calculate_pvca(exprdata_t, metadata, pct_threshold = pct_threshold, center = center, scale = scale)

  pvca_plot <- plot_pvca(pvca_result, plot_title = plot_title)

  if (plot) {
    pvca_plot <- plot_pvca(pvca_result, plot_title = plot_title)

    if (show_legend == FALSE) {
    pvca_plot <- pvca_plot + theme(legend.position = "none")
    }
  } else {
    pvca_plot <- NULL
  }

  result_list <- list(
    p = pvca_plot,
    table = pvca_result
  )

  return(result_list)

}
