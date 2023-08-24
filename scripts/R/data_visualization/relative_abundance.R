library(tidyverse)
library(envalysis)
library(lemon)
library(scales)

data <- read.csv("reactor_data.csv", header = TRUE, sep = ";")

species <- c("bbif", "bbre", "binf", "bthe", "ecol", "lsym", "paci")
colors <-  c("#0b84a5", "#f6c85f", "#9dd866", "#ffa056",
                 "#8dddd0", "#ca472f", "#6f4e7c")

data <- data %>% 
  filter(measurement %in% species)

data$experiment <- toupper(data$experiment)
data$experiment <- ifelse(data$experiment != "ALL", paste("\u0394", data$experiment, sep = ""), data$experiment)
data$experiment <- as.factor(data$experiment)
data$experiment <- relevel(data$experiment, "ALL")

data <- data %>% 
  mutate(measurement = case_when(
    measurement == "bbif" ~ "B. bifidum",
    measurement == "bbre" ~ "B. breve",
    measurement == "binf" ~ "B. infantis",
    measurement == "bthe" ~ "B. thetaiotaomicron",
    measurement == "ecol" ~ "E. coli",
    measurement == "lsym" ~ "L. symbiosum",
    measurement == "paci" ~ "P. acidilactici",
    TRUE ~ measurement
  ))

relab <- data %>% 
  ggplot(aes(x = factor(time), y = value, fill = factor(measurement))) +
  geom_bar(stat = "identity") +
  facet_rep_wrap(~ experiment, ncol = 4, repeat.tick.labels = TRUE) +
  theme_minimal() +
  scale_fill_manual(values = colors,
                    guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  labs(fill = "Species", x = "Time (h)", y = "") +
  theme_publish() +
  theme(legend.position = "top",
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(labels = percent_format(),
                     expand = c(0, 0), limits = c(0, NA))

ggsave("relative_abundance.png", relab, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')
