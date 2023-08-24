library(tidyverse)
library(rstatix)
library(ggstatsplot)
library(ComplexHeatmap)
library(circlize)
library(envalysis)

param_dist <- read.csv("parameter_distribution.csv")
param <- read.csv("parameters.csv")

dunn <- function(data, params, parameter_name) {
  filtered_data <- data %>%
    filter(parameter == parameter_name) %>% 
    filter(value < 2) %>% 
    select(-c(parameter))
  
  pwc <- filtered_data %>% 
    dunn_test(value ~ exp, p.adjust.method = "BH")
  pwc_sig4 <- pwc %>% filter(p.adj.signif == "****")
  pwc_sig4 <- pwc_sig4 %>% add_xy_position(x = "exp")
  pwc <- pwc %>% add_xy_position(x = "exp")
  
  range <- max(filtered_data$value) - min(filtered_data$value)
  
  min_ypos <- max(filtered_data$value) + (range*0.05)
  by_val <- range*0.1
  pwc_sig4 <- pwc_sig4 %>% arrange(y.position)
  nrows_sig4 <- nrow(pwc_sig4)
  new_ypos <- seq(from = min_ypos, by = by_val, length.out = nrows_sig4)
  pwc_sig4$y.position <- new_ypos
  
  if (nrow(pwc_sig4) > 0) {
    plot <- ggbetweenstats(filtered_data, x = exp, y = value,
                   type = "nonparametric",
                   p.adjust.method = "BH",
                   conf.level = 0.99,
                   pairwise.comparisons = FALSE) +
      ggsignif::geom_signif(
        comparisons = pwc_sig4$groups,
        map_signif_level = TRUE,
        annotations = pwc_sig4$p.adj.signif,
        y_position = pwc_sig4$y.position,
        test = NULL,
        na.rm = TRUE
      ) +
      labs(y = "Parameter value", x = "Experiment") +
      ggtitle(paste("Parameter:", parameter_name)) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    plot <- ggbetweenstats(filtered_data, x = exp, y = value,
                   type = "nonparametric",
                   p.adjust.method = "BH",
                   conf.level = 0.99,
                   pairwise.comparisons = FALSE) +
      labs(y = "Parameter value", x = "Experiment") +
      ggtitle(paste("Parameter:", parameter_name)) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  xlabs <- ggplot_build(plot)$layout$panel_params[[1]]$x$scale$labels
  split_vector <- strsplit(xlabs, "\n")
  expslabs <- sapply(split_vector, "[", 1)
  nlabs <- sapply(split_vector, "[", 2)
  expslabs <- toupper(expslabs)
  expslabs <- ifelse(expslabs != "ALL", paste0("\u0394", expslabs), expslabs)
  xlabs <- paste(expslabs, nlabs, sep = "\n")
  
  param_vals <- param %>% filter(parameter == parameter_name) %>% 
    pull(value) %>% round(digits = 5)
  xlabs_new <- paste(xlabs, "\nValue = ", param_vals, sep = "")
  
  if (!dir.exists("results/")) {
    dir.create("results/", recursive = TRUE)
  }
  
  plot <- plot + scale_x_discrete(labels = xlabs_new) +
    theme_publish() +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    ggtitle(expression("Parameter"~k[4]~"Lactose degradation rate"))
  
  
  plot_name <- paste("results/", parameter_name, "_distribution.png", sep = "")
  
  ggsave(plot_name, plot, 
         width = 30, height = 30*9/16, units = "cm", dpi = 500, bg = "white")
}

ratio <- function(data, params, parameter_name) {
  filtered_data <- data %>%
    filter(parameter == parameter_name) %>% 
    filter(value < 2) %>% 
    select(-c(parameter))
  
  pwc <- filtered_data %>% 
    dunn_test(value ~ exp, p.adjust.method = "BH")
  pwc_sig4 <- pwc %>% filter(p.adj.signif == "****")
  pwc_sig4 <- pwc_sig4 %>% add_xy_position(x = "exp")
  pwc <- pwc %>% add_xy_position(x = "exp")
  
  filtered_params <- params %>% 
    filter(parameter == parameter_name)
  
  experiments <- unique(filtered_params$exp)
  
  ratio_matrix <- matrix(NA, nrow = length(experiments), ncol = length(experiments), 
                         dimnames = list(experiments, experiments))
  
  for (i in 1:(length(experiments) - 1)) {
    for (j in (i + 1):length(experiments)) {
      exp1 <- experiments[i]
      exp2 <- experiments[j]
      value1 <- filtered_params$value[filtered_params$exp == exp1]
      value2 <- filtered_params$value[filtered_params$exp == exp2]
      ratio <- value1 / value2
      ratio_matrix[i, j] <- ratio
      ratio_matrix[j, i] <- 1 / ratio
    }
  }
  
  sig_comp1 <- pwc %>% 
    select(groups, p.adj.signif)
  
  sig_comp1 <- sig_comp1 %>%
    mutate(groups = map(groups, ~ setNames(.x, c("groups_1", "groups_2")))) %>%
    unnest_wider(groups)
  
  sig_comp2 <- sig_comp1 %>%
    mutate(groups = groups_2) %>%
    select(groups, groups_1, p.adj.signif) %>%
    rename(groups_2 = groups_1,
           groups_1 = groups)
  
  sig_comp <- rbind(sig_comp1, sig_comp2)
  
  annotation_matrix <- sig_comp %>% 
    pivot_wider(names_from = groups_1, values_from = p.adj.signif) %>% 
    arrange(groups_2) %>% 
    mutate_all(~ifelse(is.na(.), "ns", .)) %>% 
    column_to_rownames("groups_2") %>% 
    as.matrix()
  
  annotation_matrix <- gsub("ns", " ", annotation_matrix)
  
  range <- max(ratio_matrix, na.rm = TRUE) - min(ratio_matrix, na.rm = TRUE)
  
  if (range > 15) {
    ratio_matrix <- log(ratio_matrix, base = 10)
    legend_title <- expression(bold(Ratio~Log[10](k[j]/k[i])))
    color_fun = colorRamp2(c(min(ratio_matrix, na.rm = TRUE), 0, max(ratio_matrix, na.rm = TRUE)), 
                           c("#6e94e6", "white", "#a4302a"))
  } else {
    legend_title <- expression(bold(Ratio~k[j]/k[i]))
    color_fun = colorRamp2(c(min(ratio_matrix, na.rm = TRUE), 1, max(ratio_matrix, na.rm = TRUE)), 
                           c("#6e94e6", "white", "#a4302a"))
  }
  
  colnames(ratio_matrix) <- toupper(colnames(ratio_matrix))
  colnames(ratio_matrix) <- ifelse(colnames(ratio_matrix) != "ALL", paste0("\u0394",colnames(ratio_matrix)), colnames(ratio_matrix))
  rownames(ratio_matrix) <- toupper(rownames(ratio_matrix))
  rownames(ratio_matrix) <- ifelse(rownames(ratio_matrix) != "ALL", paste0("\u0394",rownames(ratio_matrix)), rownames(ratio_matrix))
  
  hm <- Heatmap(t(ratio_matrix),
          na_col = "white",
          column_names_side = "top",
          row_names_side = "left",
          column_names_rot = 0,
          column_names_centered = TRUE,
          row_names_centered = TRUE,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          col = color_fun,
          cell_fun = function(j, i, x, y, width, height, fill) {
            avg_color <- col2rgb(fill)
            if (mean(avg_color) < 130) {
              text_color <- "white"
            } else {
              text_color <- "black"
            }
            grid.rect(x, y, width = width, height = height, gp = gpar(fill = fill, col = "grey"))
            grid.text(label = annotation_matrix[i, j], x = x, y = y,
                      gp = gpar(fontface = "bold", fontsize = 13, col = text_color,
                                fontfamily = "mono"))
          },
          # row_names_gp = gpar(fontface = "bold"),
          # column_names_gp = gpar(fontface = "bold"),
          heatmap_legend_param = list(
            title = legend_title
          ),
          column_title = expression("Parameter"~k[4]~"Lactose degradation rate"),
          column_title_gp = gpar(fontface = "bold"))
  
  hm_name <- paste("results/", parameter_name, "_heatmap.png", sep = "")
  png(filename = hm_name, width = 1920, height = 1080, units = "px", bg = "white", res = 150)
  draw(hm)
  dev.off()
  
}

for (parameter in unique(param$parameter)) {
  dunn(param_dist, param, parameter)
  ratio(param_dist, param, parameter)
}
