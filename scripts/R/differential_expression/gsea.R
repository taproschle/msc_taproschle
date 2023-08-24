library(tidyverse)
library(clusterProfiler)
library(httr)
library(enrichplot)

experiments <- c("ALL", "BBIF", "BBRE", "BINF", 
                 "BTHE", "ECOL", "LSYM", "PACI")
mos <- c("B. bifidum", "B. breve", "B. infantis", "B. thetaiotaomicron",
         "E. coli", "L. symbiosum", "P. acidilactici")

data <- data.frame()

for (experiment in experiments){
  
  ref <- experiment
  exps <- experiments[experiments != ref]
  
  for (exp in exps){
    
    file <- paste0("results/", ref, "_deseq/", exp, "_vs_", ref, "_deseq_res.tsv")
    data_add <- read.table(file, header = TRUE, sep = "\t")
    
    data_add <- data_add %>% 
      mutate(experiment = exp,
             reference = ref) %>% 
      select(microorganism, experiment, reference, locus_tag,
             log2FoldChange, padj)
    
    data <- rbind(data, data_add)
    
  }
}

eggnogg <- read.csv2("eggnogg.csv", header = TRUE)
eggnogg$query <- sub("^[^|]*\\|[^|]*\\|", "", eggnogg$query)
eggnogg <- subset(eggnogg, !grepl("^fig\\|", query))

eggnogg_add <- eggnogg %>% 
  mutate(locus_tag = query) %>% 
  select(locus_tag, Description, Preferred_name, EC, KEGG_Pathway)

data <- merge(data, eggnogg_add, by = "locus_tag", all.x = TRUE)
data <- replace(data, is.na(data), "-")
data$KEGG_Pathway <- gsub("ko.....,", "", data$KEGG_Pathway)

kegg_data <- read.table(text = content(GET("https://rest.kegg.jp/list/pathway"), 
                                       "text"), sep = "\t", header = FALSE)
colnames(kegg_data) <- c("term", "name")

data_ranked <- data %>% 
  mutate(ranking = sign(log2FoldChange)*-log10(padj)) %>% 
  arrange(desc(ranking)) 

gsea_calc <- function(subdata, kegg_data){
  
  geneList <- subdata$ranking
  names(geneList) <- subdata$locus_tag
  
  subpws <- paste(unique(subdata$KEGG_Pathway[subdata$KEGG_Pathway != "-"]), collapse = ",")
  subpws <- unique(strsplit(subpws, ",")[[1]])
  
  TERM2NAME <- data.frame(term = subpws)
  TERM2NAME <- merge(TERM2NAME, kegg_data, by = "term", all.x = TRUE)
  TERM2NAME <- TERM2NAME[complete.cases(TERM2NAME), ]
  
  TERM2NAME <- TERM2NAME %>% 
    filter(!grepl("map03010", term) & !grepl("map01120", term) &
           !grepl("map03010", term) & !grepl("map00710", term) &
           !grepl("map01100", term) & !grepl("map03018", term))
  
  TERM2GENE <- subdata %>%
    select(locus_tag, KEGG_Pathway) %>% 
    separate_rows(KEGG_Pathway, sep = ",") %>%
    mutate(KEGG_Pathway = ifelse(KEGG_Pathway == "-", NA, KEGG_Pathway)) %>%
    na.omit() %>% 
    rename(term = KEGG_Pathway, gene = locus_tag) %>% 
    select(term, gene)
  
  TERM2GENE <- merge(TERM2GENE, TERM2NAME, by = "term", all.x = TRUE)
  TERM2GENE <- TERM2GENE %>% 
    filter(!is.na(name)) %>% 
    select(term, gene)
  
  gseaResult <- GSEA(
    geneList = geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 300,
    eps = 1e-10,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea"
  )
  
  return(gseaResult)
}

permutations <- expand.grid(experiments, experiments)
permutations <- permutations[permutations$Var1 != permutations$Var2, ]
colnames(permutations) <- c("experiment", "reference")
rownames(permutations) <- seq(nrow(permutations))
permutations <- as.data.frame(permutations)

if (!dir.exists("results/GSEA")) {
  dir.create("results/GSEA")
}

gsea_df <- data.frame()

for (i in seq(nrow(permutations))) {
  exp <- as.character(permutations[i,]$experiment)
  ref <- as.character(permutations[i,]$reference)
  
  if (!dir.exists(paste0("results/GSEA/", ref))) {
    dir.create(paste0("results/GSEA/", ref))
  }
  
  for (j in seq(length(mos))) {
    mo <- mos[j]
    mo_name <- gsub(". ", "_", mo)
    exp_title <- ifelse(startsWith(exp, "ALL"), exp, 
                        paste0("\u0394", exp))
    ref_title <- ifelse(startsWith(ref, "ALL"), ref, 
                        paste0("\u0394", ref))
    plot_title <- paste(mo, "in", exp_title, "with", ref_title, "as reference.")
    
    
    subdata <- data_ranked %>% 
      filter(experiment == exp & reference == ref & microorganism == mo)
    
    if (nrow(subdata) != 0){
      
      gseaResult <- gsea_calc(subdata, kegg_data)
      gseaResult_df <- gseaResult@result
      
      if (nrow(gseaResult_df) > 5){
        n = 5
        gsea_plot <- gseaplot2(gseaResult, geneSetID = 1:n, title = plot_title)
      } else if (nrow(gseaResult_df) > 1){
        n = nrow(gseaResult_df)
        gsea_plot <- gseaplot2(gseaResult, geneSetID = 1:n, title = plot_title)
      }
      
      if (exists("gsea_plot")){
        plot_name <- paste0("results/GSEA/", ref, "/", mo_name, "_", exp, ".png")
        ggsave(plot_name, gsea_plot, bg = "white", dpi = 500, 
               width = 30, height = 20, units = 'cm')
      }
      gseaResult_df <- gseaResult_df %>% 
        mutate(microorganism = mo, experiment = exp, reference = ref, .before = ID)
      
      gsea_df <- rbind(gsea_df, gseaResult_df)
    }
    
    if (exists("gsea_plot")){
      rm(gsea_plot)
    }
  }
}

write.table(gsea_df, file = "results/GSEA/gsea.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
