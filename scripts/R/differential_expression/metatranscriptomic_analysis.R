# Metatranscriptomic analysis
#### Libraries ####

# Load needed libraries
library(tidyverse)
library(stringr)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(PoiClaClu)
library(apeglm)
library(gridExtra)
library(factoextra)
library(ggrepel)
library(ggalt)
library(ggforce)
library(envalysis)

#### Composition data ####

comp_data <- read.csv("comp_data.tsv", header=TRUE, sep="\t")

genomic_features <- tibble(
  species = c("bbif", "bbre", "binf", "bthe", "lsym", "ecol", "paci"),
  num_genomic_features = c(2018, 2088, 2603, 4965, 4838, 4666, 2130),
  genome_length = c(2216167, 2269415,  2832748, 6293399, 4916964, 4639675, 2044083)
)

comp_data <- comp_data %>%
  left_join(genomic_features, by = "species") %>%
  mutate(normalized_num = num / (num_genomic_features*genome_length))

comp_data <- comp_data %>% 
  group_by(exp, replicate) %>%
  mutate(rel_num = as.numeric(normalized_num) / sum(as.numeric(normalized_num)))

comp_data$exp <- str_to_title(comp_data$exp)
comp_data$replicate <- str_to_title(comp_data$replicate)

species_mapping <- data.frame(
  species = unique(comp_data$species),
  custom_name = c("B. bifidum", "B. breve", "B. infantis", 
                  "B. thetaiotaomicron", "L. symbiosum", "E. coli", 
                  "P. acidilactici"),
  custom_color = c("#0b84a5", "#f6c85f", "#9dd866", "#ffa056",
                   "#8dddd0", "#ca472f", "#6f4e7c")
)
comp_data <- left_join(comp_data, species_mapping, by = "species")

comp_data$exp <- toupper(comp_data$exp)
comp_data$exp <- ifelse(comp_data$exp != "ALL", 
                        paste0("\u0394", comp_data$exp), comp_data$exp)
comp_data$exp <- as.factor(comp_data$exp)
comp_data$exp <- relevel(comp_data$exp, ref = "ALL")

comp <- ggplot(comp_data, aes(x = replicate, y = rel_num, fill = custom_name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ exp, ncol = 4) +
  scale_fill_manual(values = species_mapping$custom_color,
                    guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  theme_publish() +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "Species")

if (!dir.exists("results")) {
  dir.create("results")
}

ggsave("results/library_composition.png", comp, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

#### Counts data load ####

data <- read.csv("counts.tsv", header=TRUE, sep="\t", row.names = 1)

# Separate columns of mapping information
data_info <- data %>% 
  select(c("locus_tag", "product", "gene"))

data_info_lt <- data_info %>% 
  filter(locus_tag != "") %>% 
  rownames_to_column(var = "patric_id")

data_info_lt$patric_id <- gsub("fig\\|", "", data_info_lt$patric_id)

# Save table of information
write.table(data_info_lt, file = "data_info.tsv", sep = "\t", row.names = FALSE)

# Create data frame with counts
data_counts <- data %>% 
  select(-c("locus_tag", "product", "gene")) %>% 
  replace(is.na(.), 0) %>% 
  filter_all(any_vars(. != 0)) %>%
  filter(rowSums(.) >= 10)

# Data formatting
colnames(data_counts) <- toupper(colnames(data_counts))

#### PCA ####

# Set up needed data for DESeq object
colData <- data.frame(exp = as.factor(sub("*_R.", "", colnames(data_counts))),
                      row.names = colnames(data_counts))
dds <- DESeqDataSetFromMatrix(countData = data_counts,
                              colData = colData,
                              design = ~ exp)

# testing
vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = 'local')

# Identify constant and zero columns
pca.matrix <- t(assay(vsd))
constant_cols <- apply(pca.matrix, 2, function(x) all(x == x[1]))
zero_cols <- apply(pca.matrix, 2, function(x) all(x == 0))
pca.matrix <- pca.matrix[, !(constant_cols | zero_cols)]

pca.df <- prcomp(pca.matrix, scale. = TRUE, center = TRUE, rank. = 3)
n.clust <- fviz_nbclust(pca.df$x, FUNcluster = cluster::pam, k.max = 23)
pam <- eclust(pca.df$x, "pam", k = 8)
cluster <- pam$cluster

pca.df.data <- data.frame(pca.df$x)
pca.df.data$x <- pca.df.data[,1]
pca.df.data$y <- pca.df.data[,2]
pca.df.data$z <- pca.df.data[,3]
pca.df.data$group <- as.factor(cluster)
pca.df.data$exp <- gsub("_R.", "", rownames(pca.df.data))
pca.df.data$replicate <- gsub("..._", "", gsub("...._", "", rownames(pca.df.data)))
pca.ev1 <- round(get_eigenvalue(pca.df)[1,2], 1)
pca.ev2 <- round(get_eigenvalue(pca.df)[2,2], 1)
pca.ev3 <- round(get_eigenvalue(pca.df)[3,2], 1)
pca.totalvar <- pca.ev1 + pca.ev2 + pca.ev3
pca.df.data$exp <- ifelse(pca.df.data$exp != "ALL", paste("\u0394", pca.df.data$exp, sep = ""), pca.df.data$exp)

pca.xy <- pca.df.data %>% 
  ggplot(aes(x = x, y = y, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC1 - ", pca.ev1, "%", sep="")) +
  ylab(paste("PC2 - ", pca.ev2, "%", sep=""))

pca.xz <- pca.df.data %>% 
  ggplot(aes(x = x, y = z, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC1 - ", pca.ev1, "%", sep="")) +
  ylab(paste("PC3 - ", pca.ev3, "%", sep=""))

pca.yz <- pca.df.data %>% 
  ggplot(aes(x = y, y = z, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC2 - ", pca.ev2, "%", sep="")) +
  ylab(paste("PC3 - ", pca.ev3, "%", sep=""))

pca.plot.bf <- grid.arrange(pca.xy, pca.xz, pca.yz, ncol = 3)

if (!dir.exists("results")) {
  dir.create("results")
}

ggsave("results/PCA_bf.png", pca.plot.bf, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

# Remove the most distant replica
remove_outliers <- function(df) {
  df_grouped <- df %>% 
    group_by(exp)
  
  df_distance <- df_grouped %>% 
    mutate(distance = sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2 + (PC3 - mean(PC3))^2))
  
  df_max_distance <- df_distance %>% 
    filter(distance == max(distance))
  
  df_filtered <- anti_join(df, df_max_distance, by = c("exp", "replicate"))
  
  return(df_filtered)
}

# Select the columns 
data_counts <- data_counts %>% 
  select(rownames(remove_outliers(pca.df.data)))

# PCA analysis with the filtered data
colData <- data.frame(exp = as.factor(sub("*_R.", "", colnames(data_counts))),
                      row.names = colnames(data_counts))
dds <- DESeqDataSetFromMatrix(countData = data_counts,
                              colData = colData,
                              design = ~ exp)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = 'local')

# Identify constant and zero columns
pca.matrix <- t(assay(vsd))
constant_cols <- apply(pca.matrix, 2, function(x) all(x == x[1]))
zero_cols <- apply(pca.matrix, 2, function(x) all(x == 0))
pca.matrix <- pca.matrix[, !(constant_cols | zero_cols)]

pca.df <- prcomp(pca.matrix, scale. = TRUE, center = TRUE, rank. = 3)
n.clust <- fviz_nbclust(pca.df$x, FUNcluster = cluster::pam, k.max = 15)
pam <- eclust(pca.df$x, "pam", k = 8)
cluster <- pam$cluster

pca.df.data <- data.frame(pca.df$x)
pca.df.data$x <- pca.df.data[,1]
pca.df.data$y <- pca.df.data[,2]
pca.df.data$z <- pca.df.data[,3]
pca.df.data$group <- as.factor(cluster)
pca.df.data$exp <- gsub("_R.", "", rownames(pca.df.data))
pca.df.data$replicate <- gsub("..._", "", gsub("...._", "", rownames(pca.df.data)))
pca.ev1 <- round(get_eigenvalue(pca.df)[1,2], 1)
pca.ev2 <- round(get_eigenvalue(pca.df)[2,2], 1)
pca.ev3 <- round(get_eigenvalue(pca.df)[3,2], 1)
pca.totalvar <- pca.ev1 + pca.ev2 + pca.ev3
pca.df.data$exp <- ifelse(pca.df.data$exp != "ALL", paste("\u0394", pca.df.data$exp, sep = ""), pca.df.data$exp)

pca.xy <- pca.df.data %>% 
  ggplot(aes(x = x, y = y, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC1 - ", pca.ev1, "%", sep="")) +
  ylab(paste("PC2 - ", pca.ev2, "%", sep=""))

pca.xz <- pca.df.data %>% 
  ggplot(aes(x = x, y = z, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC1 - ", pca.ev1, "%", sep="")) +
  ylab(paste("PC3 - ", pca.ev3, "%", sep=""))

pca.yz <- pca.df.data %>% 
  ggplot(aes(x = y, y = z, color = group, shape = replicate)) +
  geom_mark_ellipse(aes(group = group, fill = group)) +
  geom_point(color = "black", size = 2.5) +
  geom_text_repel(aes(label = exp), fontface = "bold", show.legend = FALSE, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none", fill = "none", shape = "none") +
  xlab(paste("PC2 - ", pca.ev2, "%", sep="")) +
  ylab(paste("PC3 - ", pca.ev3, "%", sep=""))

pca.plot.af <- grid.arrange(pca.xy, pca.xz, pca.yz, ncol = 3)

ggsave("results/PCA_af.png", pca.plot.af, bg = "white", dpi = 500, 
       width = 30, height = 20, units = 'cm')

#### Differential expression analysis ####

# Table normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
locus_tags <- data_info[rownames(normalized_counts), 1]
has_lt <- str_length(locus_tags) != 0
rownames(normalized_counts) <- locus_tags
normalized_counts <- normalized_counts[has_lt, ]

# Save table
write.table(normalized_counts, file = "normalized_counts.tsv", sep = "\t", row.names = FALSE)

# Function to iterate the workflow
process_results <- function(dds, exp_coef, exp, ref, lfc_thr, padj_thr){

# Obtain results
res <- lfcShrink(dds, coef=exp_coef[exp,], type="apeglm")

# Convert results table to data frame
resOrdered <- as_tibble(res[order(res$pvalue),],
                        rownames = "ID",
                        rownames_to_column = TRUE) %>% 
  filter(!str_detect(ID, "^__")) %>% 
  mutate(microorganism = case_when(
    str_detect(ID, "^fig\\|398514\\.7") ~ "B. bifidum",
    str_detect(ID, "^fig\\|518634\\.20") ~ "B. breve",
    str_detect(ID, "^fig\\|391904\\.5") ~ "B. infantis",
    str_detect(ID, "^fig\\|226186\\.12") ~ "B. thetaiotaomicron",
    str_detect(ID, "^fig\\|742741\\.3") ~ "L. symbiosum",
    str_detect(ID, "^fig\\|511145\\.12") ~ "E. coli",
    str_detect(ID, "^fig\\|1254\\.353") ~ "P. acidilactici",
    TRUE ~ NA_character_),
    mo = case_when(
      str_detect(microorganism, "B. bifidum") ~ "BBIF",
      str_detect(microorganism, "B. breve") ~ "BBRE",
      str_detect(microorganism, "B. infantis") ~ "BINF",
      str_detect(microorganism, "B. thetaiotaomicron") ~ "BTHE",
      str_detect(microorganism, "L. symbiosum") ~ "LSYM",
      str_detect(microorganism, "E. coli") ~ "ECOL",
      str_detect(microorganism, "P. acidilactici") ~ "PACI",
      TRUE ~ NA_character_)
    ) %>% 
  na.omit() %>%  
  mutate(diffexpr = case_when(
    log2FoldChange >= lfc_thr & padj <= padj_thr ~ "Upregulated",
    log2FoldChange <= -lfc_thr & padj <= padj_thr ~ "Downregulated",
    .default = "Not Significant"
  )) %>% filter(mo != exp & mo != ref)

# Search for locus tags
resOrdered$locus_tag <- data_info[resOrdered$ID, 1]
resOrdered <- resOrdered %>% 
  filter(locus_tag != "")

up_color <- "#a4302a"
down_color <- "#6e94e6"
mycolors <- c(up_color, "gray", down_color)
names(mycolors) <- c("Upregulated", "Not Significant", "Downregulated")

if (!dir.exists("results")) {
  dir.create("results")
}

if (!dir.exists(paste("results/", ref, "_deseq", sep = ""))) {
  dir.create(paste("results/", ref, "_deseq", sep = ""))
}

# Make plots
resOrdered$diffexpr <- factor(resOrdered$diffexpr, 
                 levels = c("Upregulated", "Downregulated", "Not Significant"))

if (length(unique(resOrdered$mo)) == 6){
  shape_values = c(15, 16, 17,0, 1, 2)
} else {
  shape_values = c(15, 16, 17,0, 1)
}

exp_title <- ifelse(startsWith(exp, "ALL"), exp, 
                    paste0("\u0394", exp))
ref_title <- ifelse(startsWith(ref, "ALL"), ref, 
                    paste0("\u0394", ref))

vp <- resOrdered %>% ggplot(
  aes(x = log2FoldChange, y = -log10(padj), color = diffexpr,
      shape = microorganism)
  ) + geom_point(size = 2) +
  scale_shape_manual("Microorganism", values = shape_values) +
  scale_color_manual("Differential Expression", values = mycolors) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) + ggtitle(sprintf("%s vs %s", exp_title, ref_title)) +
  xlab(expression(Log[2](Fold~change))) +
  ylab(expression(-Log[10](P~adjusted))) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2,
                         label.theme = element_text(face = "italic", size = 10))
  )

g_counts <- resOrdered %>% 
  group_by(microorganism, diffexpr) %>% 
  summarize(count = n(), .groups = "drop") %>% 
  filter(diffexpr != "Not Significant")

g_counts_total <- g_counts %>%
  group_by(microorganism) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count)) 

g_counts$microorganism <- factor(g_counts$microorganism, levels = g_counts_total$microorganism)
g_counts$diffexpr <- factor(g_counts$diffexpr, levels = c("Upregulated", "Downregulated"))

cplot <- ggplot(g_counts, aes(x = microorganism, y = count, fill = diffexpr)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(up_color, down_color)) +
  labs(y = "Counts", x = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(face = "italic")
  )

plot <- grid.arrange(vp, cplot, ncol = 1, heights = c(0.7, 0.3))

image = paste("results/", ref, "_deseq/", exp, "_vs_", ref, "_vplot_", lfc_thr, "_", padj_thr, ".png", sep="")
ggsave(image, plot, bg = "white", dpi = 500, 
       width = 25, height = 20, units = 'cm')

# Export data
DESeq_export_filtered <- resOrdered %>% 
  filter(diffexpr != "Not Significant") %>% 
  select(microorganism, locus_tag, log2FoldChange, padj, diffexpr)

tsv_file = paste("results/", ref, "_deseq/", exp, "_vs_", ref, "_deseq_res_", lfc_thr, "_", padj_thr, ".tsv", sep="")
write.table(DESeq_export_filtered, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

DESeq_export <- resOrdered %>% 
  select(microorganism, locus_tag, log2FoldChange, padj, diffexpr)

tsv_file = paste("results/", ref, "_deseq/", exp, "_vs_", ref, "_deseq_res.tsv", sep="")
write.table(DESeq_export, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

}

# Generate all possible combinations
exp_list <- unique(gsub("_R.", "", colnames(dds)))
combinations <- expand.grid(ref = exp_list, exp = exp_list)
combinations <- combinations[combinations$ref != combinations$exp, ]

# Iterate through each combination
for (i in 1:nrow(combinations)) {
  exp <- as.character(combinations$exp[i])
  ref <- as.character(combinations$ref[i])
  
  # Set the reference
  dds_work <- dds
  dds_work$exp <- relevel(dds_work$exp, ref = ref)
  
  # Run DE analysis
  dds_work <- DESeq(dds_work)
  exp_coef = data.frame(coef = resultsNames(dds_work)[-1],
                        row.names =  gsub(paste("_vs_", ref, sep = ""), "", gsub("exp_", "", resultsNames(dds_work)[-1])))
  
  # Process the results
  process_results(dds = dds_work,
                  exp_coef = exp_coef,
                  exp = exp,
                  ref = ref,
                  lfc_thr = 2,
                  padj_thr = 0.01
                  )
}

#### HCA ####

up_color <- "#a4302a"
down_color <- "#6e94e6"

for (ref_name in exp_list){
  list <- combinations %>% filter(ref == ref_name)
  list <- as.character(list$exp)
  
  for (exp_name in list){
    tsv_file <- paste("results/", ref_name, "_deseq/", exp_name, "_vs_",
                     ref_name, "_deseq_res_2_0.01.tsv", sep = "")
    data <- read.csv(tsv_file, header=TRUE, sep="\t")
    if (exp_name == list[1]){
      deseq_data <- data.frame(lt = data$locus_tag,
                               mo = data$microorganism,
                               name = data$log2FoldChange)
      deseq_data <- deseq_data %>% rename("name" = exp_name)
    } else {
      data_join <- data.frame(lt = data$locus_tag,
                              mo = data$microorganism,
                              name = data$log2FoldChange)
      data_join <- data_join %>% rename("name" = exp_name)
      deseq_data <- merge(deseq_data, data_join, by = c("lt", "mo"), all = TRUE)
    }
  }
  
  rownames(deseq_data) <- deseq_data$lt
  row_groups <- deseq_data$mo
  deseq_data <- deseq_data[, -1:-2]
  deseq_data <- as.matrix(deseq_data)
  deseq_data[is.na(deseq_data)] <- 0
  df <- deseq_data
  colnames(df) <- ifelse(!startsWith(colnames(df), "ALL"),
                         paste0("\u0394", colnames(df)), colnames(df))
  
  ref_name_title <- ifelse(!startsWith(ref_name, "ALL"), paste0("\u0394", ref_name), ref_name)
  ref_name_title <- paste0(ref_name_title, " as reference")
  
  hm <- Heatmap(df,
                show_row_names = FALSE,
                column_names_side = "top",
                column_title = ref_name_title,
                column_title_gp = gpar(fontface = "bold", fontsize = 18),
                column_names_rot = 0,
                column_names_centered = TRUE,
                cluster_rows = cluster_within_group(t(df), row_groups),
                row_title_gp = gpar(fontface = "italic"),
                column_names_gp = gpar(fontface = "bold"),
                row_title_rot = 0,
                row_split = 7,
                border = TRUE,
                row_dend_side = "right",
                column_dend_side = "bottom",
                col = colorRamp2(c(min(df), 0, max(df)), c(down_color, "white", up_color)),
                heatmap_legend_param = list(title = "Log2FC"))
  
  hm_filename <- paste("results/", ref_name, "_HCA.png", sep = "")
  png(filename = hm_filename, width = 25, height = 20, units = "cm", bg = "white", res = 500)
  draw(hm)
  dev.off()
}

ggsave(image, plot, bg = "white", dpi = 500, 
       width = 25, height = 20, units = 'cm')