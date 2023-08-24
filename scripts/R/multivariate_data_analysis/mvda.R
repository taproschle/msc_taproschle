# MVDA
#### Libraries & data ####

library(tidyverse)
library(factoextra)
library(ggrepel)
library(ggforce)
library(ggalt)
library(grid)
library(gridExtra)
library(Hmisc)

data <- read.csv("reactor_data.csv", header = TRUE, sep = ";")

#### PCA ####

neuac_data_exp <- data %>% 
  filter(measurement == "neuac" & (time == 8 | time == 12)) %>%
  group_by(experiment, measurement) %>% 
  mutate(mean_value = sum(value)/n(), .after = value,
         time = 10) %>% 
  select(-c(value, replicate)) %>% 
  unique()

neuac_data_sta <- data %>% 
  filter(measurement == "neuac" & (time == 16 | time == 20)) %>%
  group_by(experiment, measurement) %>% 
  mutate(mean_value = sum(value)/n(), .after = value,
         time = 18) %>% 
  select(-c(value, replicate)) %>% 
  unique()

# PCA-ready data frame
df.exp <- data %>% 
  group_by(experiment, measurement, time) %>% 
  mutate(mean_value = sum(value)/n(), .after = value) %>% 
  select(-c(value, replicate)) %>% 
  unique() %>% 
  filter(time == 10) %>%
  bind_rows(neuac_data_exp) %>% 
  arrange(experiment, measurement) %>% 
  select(-c(unit)) %>% 
  pivot_wider(names_from = experiment, values_from = mean_value) %>% 
  ungroup(c(time)) %>% 
  select(-c(time)) %>% 
  column_to_rownames(var = "measurement") %>% 
  as.data.frame.matrix() %>% 
  t()

df.sta <- data %>% 
  group_by(experiment, measurement, time) %>% 
  mutate(mean_value = sum(value)/n(), .after = value) %>% 
  select(-c(value, replicate)) %>% 
  unique() %>%
  filter(time == 18) %>%
  # Produces null column in PCA
  filter(measurement != "2fl" & measurement != "3sl" & measurement != "lnt") %>% 
  bind_rows(neuac_data_sta) %>% 
  arrange(experiment, measurement) %>% 
  select(-c(unit)) %>% 
  pivot_wider(names_from = experiment, values_from = mean_value) %>% 
  ungroup(c(time)) %>% 
  select(-c(time)) %>% 
  column_to_rownames(var = "measurement") %>% 
  as.data.frame.matrix() %>% 
  t()

# Exponential time
pca.df.exp <- prcomp(df.exp, scale. = TRUE, center = TRUE)
pca.exp.vars <- fviz_pca_var(pca.df.exp, col.var = "steelblue") +
  theme_classic() +
  ggtitle("")

if (!dir.exists("figures")) {
  dir.create("figures")
}

ggsave("figures/pca_exp_vars.png", pca.exp.vars, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

library(paran)
paran(df.exp, iterations = 10000, quietly = FALSE,
      status = FALSE, all = TRUE, cfa = FALSE, graph = TRUE,
      color = TRUE, col = c("black", "red", "blue"),
      lty = c(1,2,3), lwd = 1, legend = TRUE, file = "",
      width = 640, height = 640, grdevice = "png", seed = 0,
      mat = NA, n = NA)

fviz_contrib(pca.df.exp, "var", axes = 1, xtickslab.rt =  90)
fviz_contrib(pca.df.exp, "var", axes = 2, xtickslab.rt =  90)
fviz_contrib(pca.df.exp, "var", axes = 3, xtickslab.rt =  90)

# Quantity of ranks defined by the paran function is 3
pca.df.exp <- prcomp(df.exp, scale. = TRUE, center = TRUE, rank. = 3)
n.clust.exp <- fviz_nbclust(pca.df.exp$x, FUNcluster = cluster::pam, k.max = 7,
                            method = "silhouette")

plot.pca.vars.nclust.exp <- grid.arrange(pca.exp.vars, n.clust.exp, ncol = 2, 
             top = textGrob("Exponential", gp = gpar(fontsize = 15, font = 2)))

ggsave("figures/vars_nclust_exp.png", plot.pca.vars.nclust.exp, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

pam <- eclust(pca.df.exp$x, "pam", k = 5)
fviz_silhouette(pam)
cluster.exp <- pam$cluster

pca.df.exp.data <- data.frame(pca.df.exp$x)
pca.df.exp.data$x <- pca.df.exp.data[,1]
pca.df.exp.data$y <- pca.df.exp.data[,2]
pca.df.exp.data$z <- pca.df.exp.data[,3]
pca.df.exp.data$group <- as.factor(cluster.exp)
pca.df.exp.data$label <- rownames(pca.df.exp.data)
pca.exp.ev1 <- round(get_eigenvalue(pca.df.exp)[1,2], 1)
pca.exp.ev2 <- round(get_eigenvalue(pca.df.exp)[2,2], 1)
pca.exp.ev3 <- round(get_eigenvalue(pca.df.exp)[3,2], 1)
pca.exp.totalvar <- pca.exp.ev1 + pca.exp.ev2 + pca.exp.ev3

pca.df.exp.data <- as_tibble(pca.df.exp.data)
pca.df.exp.data$label <- toupper(pca.df.exp.data$label)
pca.df.exp.data$label <- ifelse(pca.df.exp.data$label != "ALL", paste("\u0394", pca.df.exp.data$label, sep = ""), pca.df.exp.data$label)

pca.exp.xy <- pca.df.exp.data %>% 
  ggplot(aes(x = x, y = y, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC1 - ", pca.exp.ev1, "%", sep="")) +
  ylab(paste("PC2 - ", pca.exp.ev2, "%", sep=""))

pca.exp.xz <- pca.df.exp.data %>% 
  ggplot(aes(x = x, y = z, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC1 - ", pca.exp.ev1, "%", sep="")) +
  ylab(paste("PC3 - ", pca.exp.ev3, "%", sep=""))

pca.exp.yz <- pca.df.exp.data %>% 
  ggplot(aes(x = y, y = z, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC2 - ", pca.exp.ev2, "%", sep="")) +
  ylab(paste("PC3 - ", pca.exp.ev3, "%", sep=""))

pca.plot.exp <- grid.arrange(pca.exp.xy, pca.exp.xz, pca.exp.yz, ncol = 3)

ggsave("figures/pca_plot_exp.png", pca.plot.exp, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

# Stationary time
pca.df.sta <- prcomp(df.sta, scale. = TRUE, center = TRUE)
pca.sta.vars <- fviz_pca_var(pca.df.sta, col.var = "steelblue")

ggsave("figures/pca_sta_vars.png", pca.sta.vars, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

paran(df.sta, iterations = 10000, quietly = FALSE,
      status = FALSE, all = TRUE, cfa = FALSE, graph = TRUE,
      color = TRUE, col = c("black", "red", "blue"),
      lty = c(1,2,3), lwd = 1, legend = TRUE, file = "",
      width = 640, height = 640, grdevice = "png", seed = 0,
      mat = NA, n = NA)

fviz_contrib(pca.df.sta, "var", axes = 1, xtickslab.rt =  90)
fviz_contrib(pca.df.sta, "var", axes = 2, xtickslab.rt =  90)
fviz_contrib(pca.df.sta, "var", axes = 3, xtickslab.rt =  90)

# Quantity of ranks defined by the paran function is 2 but 3 will be used
pca.df.sta <- prcomp(df.sta, scale. = TRUE, center = TRUE, rank. = 3)
n.clust.sta <- fviz_nbclust(pca.df.sta$x, FUNcluster = cluster::pam, k.max = 7)

ggsave("figures/n_clust_sta.png", n.clust.sta, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

plot.pca.vars.nclust.sta <- grid.arrange(pca.sta.vars, n.clust.sta, ncol = 2, 
                                     top = textGrob("Stationary", gp = gpar(fontsize = 15, font = 2)))

ggsave("figures/vars_nclust_sta.png", plot.pca.vars.nclust.sta, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

pam <- eclust(pca.df.sta$x, "pam", k = 2)
fviz_silhouette(pam)
cluster.sta <- pam$cluster

pca.df.sta.data <- data.frame(pca.df.sta$x)
pca.df.sta.data$x <- pca.df.sta.data[,1]
pca.df.sta.data$y <- pca.df.sta.data[,2]
pca.df.sta.data$z <- pca.df.sta.data[,3]
pca.df.sta.data$group <- as.factor(cluster.sta)
pca.df.sta.data$label <- rownames(pca.df.sta.data)
pca.sta.ev1 <- round(get_eigenvalue(pca.df.sta)[1,2], 1)
pca.sta.ev2 <- round(get_eigenvalue(pca.df.sta)[2,2], 1)
pca.sta.ev3 <- round(get_eigenvalue(pca.df.sta)[3,2], 1)
pca.sta.totalvar <- pca.sta.ev1 + pca.sta.ev2 + pca.sta.ev3

pca.df.sta.data <- as_tibble(pca.df.sta.data)
pca.df.sta.data$label <- toupper(pca.df.sta.data$label)
pca.df.sta.data$label <- ifelse(pca.df.sta.data$label != "ALL", paste("\u0394", pca.df.sta.data$label, sep = ""), pca.df.sta.data$label)


pca.sta.xy <- pca.df.sta.data %>% 
  ggplot(aes(x = x, y = y, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC1 - ", pca.sta.ev1, "%", sep="")) +
  ylab(paste("PC2 - ", pca.sta.ev2, "%", sep=""))

pca.sta.xz <- pca.df.sta.data %>% 
  ggplot(aes(x = x, y = z, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC1 - ", pca.sta.ev1, "%", sep="")) +
  ylab(paste("PC3 - ", pca.sta.ev3, "%", sep=""))

pca.sta.yz <- pca.df.sta.data %>% 
  ggplot(aes(x = y, y = z, color = group)) +
  geom_mark_ellipse(aes(fill = group)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = label), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none") +
  xlab(paste("PC2 - ", pca.sta.ev2, "%", sep="")) +
  ylab(paste("PC3 - ", pca.sta.ev3, "%", sep=""))

pca.plot.sta <- grid.arrange(pca.sta.xy, pca.sta.xz, pca.sta.yz, ncol = 3)

ggsave("figures/pca_plot_sta.png", pca.plot.sta, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

#### Correlation network ####
# Spearman correlation at 5% significance level
cor_exp <- rcorr(df.exp, type = "spearman")
cor_exp.0.05 <- cor_exp$r * (cor_exp$P < 0.05)
cor_exp.0.05[is.na(cor_exp.0.05)] <- 1

cor_exp.0.05 <- reshape2::melt(cor_exp.0.05)
names(cor_exp.0.05) <- c("Source", "Target", "Correlation")

cor_exp.0.05 <- cor_exp.0.05 %>% 
  filter(abs(Correlation) >= 0.6)
write.table(cor_exp.0.05, file = "cor_exp_cys_0.05.csv", sep = ";", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

cor_sta <- rcorr(df.sta, type = "spearman")
cor_sta.0.05 <- cor_sta$r * (cor_sta$P < 0.05)
cor_sta.0.05[is.na(cor_sta.0.05)] <- 1

cor_sta.0.05 <- reshape2::melt(cor_sta.0.05)
names(cor_sta.0.05) <- c("Source", "Target", "Correlation")

cor_sta.0.05 <- cor_sta.0.05 %>% 
  filter(abs(Correlation) >= 0.6)
write.table(cor_sta.0.05, file = "cor_sta_cys_0.05.csv", sep = ";", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
