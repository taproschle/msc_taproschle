#### Libraries & data ####

library(tidyverse)
library(gridExtra)
library(grid)
library(envalysis)
library(lemon)

data <- read.csv("reactor_data.csv", header = TRUE, sep = ";")

#### Growth data ####

experiments <- unique(data$experiment)

# Zwietering modification of Gompertz equation
gompertz_zwietering <- function(t, A, mu_max, lambda, B){
  B + (A - B) * exp(-exp(((mu_max * exp(1))/(A - B)) * (lambda - t) + 1))
}

fit.table <- data.frame()
pred.df <- data.frame()

for (exp_name in experiments) {
  # Select the data
  growth <- data %>% 
    filter(measurement == "biomass") %>% 
    filter(experiment == exp_name) %>% 
    select(time, value, replicate)
  
  for (i in seq(1:2)) {
    growth_i <- growth %>% 
      filter(replicate == i) %>% 
      select(time, value)
    
    # Fit the model using nonlinear regression
    fit <- nls(value ~ gompertz_zwietering(time, A, mu_max, lambda, B = growth_i$value[1]),
               data = growth_i,
               start = list(A = 2, mu_max = 0.3, lambda = 7))
    
    # Extract the estimated parameter values
    A <- coef(fit)["A"]
    mu_max <- coef(fit)["mu_max"]
    lambda <- coef(fit)["lambda"]
    B <- growth$value[1]
    t_i <- lambda + (A - B)/(exp(1)*mu_max)
    mid <-(A - B)/2
    
    # Calculate the root mean squared error
    rmse <- sqrt(mean(residuals(fit)^2))
    
    # Get the final and maximum OD value
    final_od <- tail(growth_i$value, n = 1)
    max_od <- max(growth_i$value)
    
    # Create a data frame
    data_add <- data.frame(experiment = exp_name,
                           mu_max = mu_max,
                           lambda = lambda,
                           t_i = t_i,
                           mid = mid,
                           max_od = max_od,
                           rmse = rmse,
                           replicate = i,
                           row.names = NULL)
    
    fit.table <- rbind(fit.table, data_add)
    
    # # Predicted values
    # time_seq <- seq(0, 24, by = 0.1)
    # fitted_values <- predict(fit, newdata = data.frame(time = time_seq))
    # 
    # pred_add <- data.frame(experiment = exp_name,
    #                        time = time_seq,
    #                        value = fitted_values,
    #                        replicate = i)
    # 
    # pred.df <- rbind(pred.df, pred_add)
    
    
  }
}

write.table(fit.table, "growth_table.csv", sep = ";", row.names = FALSE,
            quote = FALSE)

growthSE <- data %>% 
  filter(measurement == "biomass") %>% 
  group_by(experiment, time) %>% 
  mutate(
    mean_value = mean(value),
    se = sd(value)/n()
  ) %>% 
  select(experiment, time, mean_value, se) %>% 
  unique()

fitSE <- fit.table %>% 
  group_by(experiment) %>% 
  mutate(
    mu_max_mv = mean(mu_max),
    mu_max_se = sd(mu_max)/n(),
    lambda_mv = mean(lambda),
    lambda_se = sd(lambda)/n(),
    max_od_mv = mean(max_od),
    max_od_se = sd(max_od)/n(),
  ) %>% 
  select(experiment,
         mu_max_mv, mu_max_se,
         lambda_mv, lambda_se,
         max_od_mv, max_od_se) %>% 
  unique()

custom_colors = c("black", "#0b84a5", "#f6c85f", "#9dd866", "#ffa056",
                  "#8dddd0", "#ca472f", "#6f4e7c")

growthSE$experiment <- toupper(growthSE$experiment)
growthSE$experiment <- ifelse(growthSE$experiment != "ALL", paste0("\u0394", growthSE$experiment), growthSE$experiment)
growthSE$experiment <- relevel(as.factor(growthSE$experiment), ref = "ALL")

growth_plot <- ggplot(growthSE, aes(x = time, y = mean_value, color = experiment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.5) +
  labs(x = "Time", y = "Value", color = "Experiment") +
  theme_publish() +
  facet_rep_wrap(~experiment, ncol = 4, repeat.tick.labels = TRUE) +
  scale_color_manual(values = custom_colors) +
  guides(color = "none") +
  labs(x = "Time (h)", y = "OD600")

fitSE_mu_max <- fitSE %>%
  select(experiment, mu_max_mv, mu_max_se) %>% 
  arrange(desc(mu_max_mv)) %>% 
  arrange(desc(experiment == "all")) %>% 
  mutate(experiment = as.factor(experiment),
         label = experiment)

fitSE_mu_max$label <- toupper(fitSE_mu_max$label)
fitSE_mu_max$label <- ifelse(fitSE_mu_max$label != "ALL", paste0("\u0394", fitSE_mu_max$label), fitSE_mu_max$label)

fitSE_mu_max$experiment <- factor(fitSE_mu_max$experiment,
                                  levels = fitSE_mu_max$experiment)

fitSE_max_od <- fitSE %>%
  select(experiment, max_od_mv, max_od_se) %>% 
  arrange(desc(max_od_mv)) %>% 
  arrange(desc(experiment == "all")) %>% 
  mutate(experiment = as.factor(experiment),
         label = experiment)

fitSE_max_od$label <- toupper(fitSE_max_od$label)
fitSE_max_od$label <- ifelse(fitSE_max_od$label != "ALL", paste0("\u0394", fitSE_max_od$label), fitSE_max_od$label)

fitSE_max_od$experiment <- factor(fitSE_max_od$experiment,
                                  levels = fitSE_max_od$experiment)

fitSE_lambda <- fitSE %>%
  select(experiment, lambda_mv, lambda_se) %>% 
  arrange(lambda_mv) %>% 
  arrange(desc(experiment == "all")) %>% 
  mutate(experiment = as.factor(experiment),
         label = experiment)

fitSE_lambda$label <- toupper(fitSE_lambda$label)
fitSE_lambda$label <- ifelse(fitSE_lambda$label != "ALL", paste0("\u0394", fitSE_lambda$label), fitSE_lambda$label)

fitSE_lambda$experiment <- factor(fitSE_lambda$experiment,
                                  levels = fitSE_lambda$experiment)

mu_max_table <- as.matrix(TukeyHSD(aov(mu_max ~ experiment, data = fit.table))$experiment)
vector0.05 <- mu_max_table[,4] < 0.05
vector0.01 <- mu_max_table[,4] < 0.01
pairs0.05 <- rownames(table)[vector0.05]
pairs0.01 <- rownames(table)[vector0.01]
pairs0.05 <- pairs0.05[!pairs0.05 %in% pairs0.01]
df0.01 <- data.frame(pairs = pairs0.01, sig = "**")
df0.05 <- data.frame(pairs = pairs0.05, sig = "*")
sig_df <- rbind(df0.01, df0.05)

df <- data.frame()

for (i in seq(nrow(sig_df))) {
  experiment <- strsplit(sig_df[i,]$pairs, "-")[[1]][1]
  exp <- experiment
  color <- strsplit(sig_df[i,]$pairs, "-")[[1]][2]
  sig <- sig_df[i,]$sig
  pos <- sum(df$experiment == experiment) + 1
  
  mu_max_mv <- filter(fitSE_mu_max, experiment == exp)$mu_max_mv
  mu_max_se <- filter(fitSE_mu_max, experiment == exp)$mu_max_se
  max_height <- max(fitSE_mu_max$mu_max_mv + fitSE_mu_max$mu_max_se)
  
  y <- mu_max_mv + mu_max_se + .05*(max_height + pos)
  
  df_add <- data.frame(experiment = experiment, color = color, sig = sig, pos = pos, y = y)
  df <- rbind(df, df_add)
}

df$experiment <- as.factor(df$experiment)
df$color <- as.factor(df$color)
df <- as_tibble(df)

mu_max <- ggplot(fitSE_mu_max, aes(x = experiment, y = mu_max_mv, color = experiment)) +
  geom_errorbar(aes(ymin = mu_max_mv,
                    ymax = mu_max_mv + mu_max_se),
                width = 0.4,
                color = "black") +
  geom_col(fill = "white", color = "black") +
  geom_text(aes(label = label, y = mu_max_mv/2,
                fontface = "bold")) +
  geom_text(aes(label = round(mu_max_mv, 2),
                y = mu_max_mv + mu_max_se + .05*max(mu_max_mv + mu_max_se),
                fontface = "bold")) +
  theme_void() +
  labs(caption = expression(μ[max])) +
  theme(strip.text.x = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("all" = "black",
                                "bbif" = "#0b84a5",
                                "bbre" = "#f6c85f",
                                "binf" = "#9dd866",
                                "bthe" = "#ffa056",
                                "ecol" = "#8dddd0",
                                "lsym" = "#ca472f",
                                "paci" = "#6f4e7c")) +
  geom_text(data = df, 
            aes(x = experiment, y = y, label = sig, color = color, 
                fontface = "bold"), size = 10)

max_od <- ggplot(fitSE_max_od, aes(x = experiment, y = max_od_mv, color = experiment)) +
  geom_errorbar(aes(ymin = max_od_mv,
                    ymax = max_od_mv + max_od_se),
                width = 0.4,
                color = "black") +
  geom_col(fill = "white", color = "black") +
  geom_text(aes(label = label, y = max_od_mv/2,
                fontface = "bold")) +
  geom_text(aes(label = round(max_od_mv, 2),
                y = max_od_mv + max_od_se + .05*max(max_od_mv + max_od_se),
                fontface = "bold")) +
  theme_void() +
  labs(caption = expression(OD600[max])) +
  theme(strip.text.x = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("all" = "black",
                                "bbif" = "#0b84a5",
                                "bbre" = "#f6c85f",
                                "binf" = "#9dd866",
                                "bthe" = "#ffa056",
                                "ecol" = "#8dddd0",
                                "lsym" = "#ca472f",
                                "paci" = "#6f4e7c"))

lambda_table <- as.matrix(TukeyHSD(aov(lambda ~ experiment, data = fit.table))$experiment)
vector0.05 <- lambda_table[,4] < 0.05
vector0.01 <- lambda_table[,4] < 0.01
pairs0.05 <- rownames(table)[vector0.05]
pairs0.01 <- rownames(table)[vector0.01]
pairs0.05 <- pairs0.05[!pairs0.05 %in% pairs0.01]
df0.01 <- data.frame(pairs = pairs0.01, sig = "**")
df0.05 <- data.frame(pairs = pairs0.05, sig = "*")
sig_df <- rbind(df0.01, df0.05)

df <- data.frame()

for (i in seq(nrow(sig_df))) {
  experiment <- strsplit(sig_df[i,]$pairs, "-")[[1]][1]
  exp <- experiment
  color <- strsplit(sig_df[i,]$pairs, "-")[[1]][2]
  sig <- sig_df[i,]$sig
  pos <- sum(df$experiment == experiment) + 1
  
  lambda_mv <- filter(fitSE_lambda, experiment == exp)$lambda_mv
  lambda_se <- filter(fitSE_lambda, experiment == exp)$lambda_se
  max_height <- max(fitSE_lambda$lambda_mv + fitSE_lambda$lambda_se)
  
  y <- lambda_mv + lambda_se + .05*(max_height + 15*pos)
  
  df_add <- data.frame(experiment = experiment, color = color, sig = sig, pos = pos, y = y)
  df <- rbind(df, df_add)
}

df$experiment <- as.factor(df$experiment)
df$color <- as.factor(df$color)
df <- as_tibble(df)

lambda <- ggplot(fitSE_lambda, aes(x = experiment, y = lambda_mv, color = experiment)) +
  geom_errorbar(aes(ymin = lambda_mv,
                    ymax = lambda_mv + lambda_se),
                width = 0.4,
                color = "black") +
  geom_col(fill = "white", color = "black") +
  geom_text(aes(label = label, y = lambda_mv/2,
                fontface = "bold")) +
  geom_text(aes(label = round(lambda_mv, 2),
                y = lambda_mv + lambda_se + .05*max(lambda_mv + lambda_se),
                fontface = "bold")) +
  theme_void() +
  labs(caption = expression(T[lag])) +
  theme(strip.text.x = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("all" = "black",
                                "bbif" = "#0b84a5",
                                "bbre" = "#f6c85f",
                                "binf" = "#9dd866",
                                "bthe" = "#ffa056",
                                "ecol" = "#8dddd0",
                                "lsym" = "#ca472f",
                                "paci" = "#6f4e7c")) +
  geom_text(data = df, 
            aes(x = experiment, y = y, label = sig, color = color, 
                fontface = "bold"), size = 10)

param_plot <- grid.arrange(mu_max, max_od, lambda, ncol = 3)
growth_param_plot <- grid.arrange(growth_plot, param_plot, ncol = 1, heights = c(0.6, 0.4))

ggsave("growth_param_plot.png", growth_param_plot, bg = "white", dpi = 500, 
       width = 45, height = 25.3125, units = 'cm')

#### Mets ####

mets_data <- data %>% 
  filter(unit == "g/L")

mets_data$experiment <- toupper(mets_data$experiment)
mets_data$experiment <- ifelse(mets_data$experiment != "ALL", paste("\u0394", mets_data$experiment, sep = ""), mets_data$experiment)
mets_data$experiment <- relevel(as.factor(mets_data$experiment), ref = "ALL")
mets_data$measurement <- as.factor(mets_data$measurement)
levels(mets_data$measurement) <- c("2FL", "3SL", "Acetate", "Butyrate", "Fucose",
                                   "Galactose", "Glucose", "Lactate", "Lactose",
                                   "LNT", "Neu5Ac", "Succinate")
mets_data$measurement <- factor(mets_data$measurement,
                                levels = c("LNT", "2FL", "3SL", "Lactose", 
                                           "Galactose", "Glucose", "Fucose", "Neu5Ac",
                                           "Acetate", "Butyrate", "Lactate", "Succinate"))

mets_data <- mets_data %>%
  mutate(value = if_else(measurement == "Acetate", 1000 * value / 59.044, value),
         unit = if_else(measurement == "Lactate", "mmol/L", unit)) %>% 
  mutate(value = if_else(measurement == "Lactate", 1000 * value / 90.08, value),
         unit = if_else(measurement == "Butyrate", "mmol/L", unit)) %>% 
  mutate(value = if_else(measurement == "Butyrate", 1000 * value / 88.11, value),
         unit = if_else(measurement == "Butyrate", "mmol/L", unit)) %>% 
  mutate(value = if_else(measurement == "Succinate", 1000 * value / 118.09, value),
         unit = if_else(measurement == "Succinate", "mmol/L", unit))

metsSE <- mets_data %>% 
  group_by(experiment, time, measurement) %>% 
  mutate(
    mean_value = mean(value),
    se = sd(value)/n()
  ) %>% 
  select(experiment, time, measurement, mean_value, se) %>% 
  unique()

mets_plot <- metsSE %>% 
  filter(!measurement %in% c("Acetate", "Lactate", "Butyrate", "Succinate", "Neu5Ac")) %>% 
  ggplot(aes(x = time, y = mean_value, color = measurement)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.5) +
  facet_rep_wrap(~ experiment, ncol = 4, repeat.tick.labels = TRUE) +
  labs(x = "Time (h)", y = "g/L", color = "Metabolite") +
  theme_publish() +
  theme(legend.position = "top")

ggsave("mets_plot.png", mets_plot, bg = "white", dpi = 500, 
       width = 45, height = 25.3125, units = 'cm')

neuac_plot <- metsSE %>% 
  filter(measurement == "Neu5Ac") %>% 
  ggplot(aes(x = time, y = mean_value, color = measurement)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), width = 0.5) +
  facet_rep_wrap(~ experiment, ncol = 4, repeat.tick.labels = TRUE) +
  labs(x = "Time (h)", y = "g/L", color = "Metabolite") +
  theme_publish() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

ggsave("neuac_plot.png", neuac_plot, bg = "white", dpi = 500, 
       width = 45, height = 25.3125, units = 'cm')

acetate <- metsSE %>% 
  filter(measurement == "Acetate") %>%
  mutate(growth = ifelse(time == 10, "Exponential", "Stationary")) %>% 
  ggplot(aes(x = experiment, y = mean_value, fill = growth)) +
  geom_errorbar(aes(ymin = mean_value, ymax = mean_value + se), width = 0.5) +
  geom_col(color = "black") +
  facet_rep_wrap(~ growth, ncol = 2, repeat.tick.labels = TRUE) +
  ggtitle("Acetate") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_publish() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "mmol/L")

lactate <- metsSE %>% 
  filter(measurement == "Lactate") %>%
  mutate(growth = ifelse(time == 10, "Exponential", "Stationary")) %>% 
  ggplot(aes(x = experiment, y = mean_value, fill = growth)) +
  geom_errorbar(aes(ymin = mean_value, ymax = mean_value + se), width = 0.5) +
  geom_col(color = "black") +
  facet_rep_wrap(~ growth, ncol = 2, repeat.tick.labels = TRUE) +
  ggtitle("Lactate") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_publish() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "mmol/L")

butyrate <- metsSE %>% 
  filter(measurement == "Butyrate") %>%
  mutate(growth = ifelse(time == 10, "Exponential", "Stationary")) %>% 
  ggplot(aes(x = experiment, y = mean_value, fill = growth)) +
  geom_errorbar(aes(ymin = mean_value, ymax = mean_value + se), width = 0.5) +
  geom_col(color = "black") +
  facet_rep_wrap(~ growth, ncol = 2, repeat.tick.labels = TRUE) +
  ggtitle("Butyrate") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_publish() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "mmol/L")

succinate <- metsSE %>% 
  filter(measurement == "Succinate") %>%
  mutate(growth = ifelse(time == 10, "Exponential", "Stationary")) %>% 
  ggplot(aes(x = experiment, y = mean_value, fill = growth)) +
  geom_errorbar(aes(ymin = mean_value, ymax = mean_value + se), width = 0.5) +
  geom_col(color = "black") +
  facet_rep_wrap(~ growth, ncol = 2, repeat.tick.labels = TRUE) +
  ggtitle("Succinate") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_publish() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = "mmol/L")

hplc <- grid.arrange(acetate, lactate, butyrate, succinate, ncol = 2)

ggsave("hplc.png", hplc, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')
