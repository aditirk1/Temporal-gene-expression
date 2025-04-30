library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(edgeR)
library(DESeq2)
library(PCAtools)
library(Glimma)
library(DGEobj.utils)
library(ComplexHeatmap)
library(Hmisc)
library(corrplot)
library(gplots)
library(RColorBrewer)
library(NMF)
library(ashr)
library(vsn)
library(tibble)
library(DEGreport)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library(ReactomePA)
library(sva)
# List of IDs for your datasets
ids <- c("GSE132369", "GSE132370", "GSE130727", "GSE158743", "GSE258291", "A549") 

# Get the list of directories matching the pattern 
dirs <- list.dirs(path = "C:/Users/Aditi//Desktop/Work/Hypoxia_EMT_mit/R/DGE Analysis Main/data", recursive = FALSE, full.names = TRUE)

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each directory
for (dir in dirs) {
  # Identify the text file in the current directory
  file <- list.files(path = dir, pattern = "*.tsv", full.names = TRUE)
  
  # Ensure there's exactly one file in the directory
  if (length(file) == 1) {
    # Read the dataset
    dataset <- read.table(
      file = file,
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
    # Append to the list
    data_frames[[basename(dir)]] <- dataset
  } else {
    warning(sprintf("Skipped directory: %s (no or multiple .tsv files found)", dir))
  }
}

# Assuming your list of datasets is called 'dataset_list'
# Merge all datasets preserving row names
merge_datasets <- function(data_frames) {
  # Get all unique row names across all datasets
  all_genes <- unique(unlist(lapply(data_frames, rownames)))
  
  # Initialize an empty matrix with the right dimensions
  result <- matrix(NA, 
                   nrow = length(all_genes), 
                   ncol = sum(sapply(data_frames, ncol)))
  
  # Set row names
  rownames(result) <- all_genes
  colnames(result) <- unlist(lapply(data_frames, colnames))
  result <- as.data.frame(result)
  current_col <- 1
  for (dataset in data_frames) {
    for (col in 1:ncol(dataset)) {
      result[rownames(dataset), current_col] <- dataset[, col]
      current_col <- current_col + 1
    }
  }
  
  return(result)
}
merged_data <- merge_datasets(data_frames)
dim(merged_data)

library(tidyverse)

# Create the metadata dataframe
metadata <- data.frame(
  data_frames = character(),
  cell_type = character(),
  time = character(),
  stringsAsFactors = FALSE
)

# Load necessary library
library(dplyr)

# Read the metadata for GSE81513
gse81513_meta <- data.frame(
  sample_id = colnames(data_frames[["GSE81513"]]),
  dataset = "GSE81513",
  cell_type = "HCT116",
  time = as.numeric(c(rep(c(0, 1, 2, 24), each = 3))),
  stringsAsFactors = FALSE
)
# Create metadata for GSE53510
gse53510_meta <- data.frame(
  sample_id = colnames(data_frames[["GSE53510"]]),
  dataset = "GSE53510",
  cell_type = "HUVECs",
  time = as.numeric(c(rep(0, 3), rep(8, 3), rep(24, 3), rep(48, 2))),
  stringsAsFactors = FALSE
)
# Create metadata for GSE100099
gse100099_meta <- data.frame(
  sample_id = c(paste0("GSM26712", 16:28), "GSM2671285", "GSM2671288", paste0("GSM26712", 38:50), "GSM2671284", "GSM2671286", "GSM2671287", "GSM2671215", "GSM2671259", "GSM2671237", "GSM2671271", "GSM2671283"),
  dataset = "GSE100099",
  cell_type = "MCF7",
  time = as.numeric(c(1:12, 24, 2.5, 7.5, 1:12, 24, 1, 4, 5, rep(0, 5))),
  stringsAsFactors = FALSE
)

# Create metadata for GSE132624
gse132624_meta <- data.frame(
  sample_id = c(paste0("GSM38823", 37:48), paste0("GSM38823", 49:60), paste0("GSM38823", 61:66), paste0("GSM38823", c(67, 69)), paste0("GSM38823", 70:72)),
  dataset = "GSE132624",
  cell_type = c(rep("501", 12), rep("IGR37", 12), rep("IGR39", 11)),
  time = as.numeric(c(rep(c(0, 12, 24, 48), each = 3), rep(c(0, 12, 24, 48), each = 3), rep(c(0, 12), each = 3), rep(24, 2), rep(48, 3))),
  stringsAsFactors = FALSE
)

# Create metadata for GSE95280
gse95280_meta <- data.frame(
  sample_id = paste0("GSM2501", 230:256),
  dataset = "GSE95280",
  cell_type = c(rep("501", 12), rep("IGR37", 6), rep("IGR39", 9)),
  time = as.numeric(c(rep(c(0, 12, 24, 48), each = 3), rep(c(0, 12), each = 3), rep(c(0, 12, 24), each = 3))),
  stringsAsFactors = FALSE
)

# Combine all metadata data frames
combined_meta <- bind_rows(gse81513_meta, gse53510_meta, gse132624_meta, gse95280_meta, gse100099_meta)


# Convert sample_id column into row names for combined_metadata
combined_meta <- combined_meta %>%
  column_to_rownames(var = "sample_id")

# Write the combined metadata to a CSV file
write.csv(combined_meta, file = "combined_metadata.csv", row.names = TRUE)

# Keep only the columns in merged_data that are row names of combined_metadata
merged_data <- merged_data[, rownames(combined_meta)]
all(colnames(merged_data) %in% rownames(combined_meta))
combined_meta$time <- factor(combined_meta$time)
combined_meta$dataset <- factor(combined_meta$dataset)
combined_meta$cell_type <- factor(combined_meta$cell_type)
write.csv(merged_data, "merged_data.csv", row.names = TRUE)


RawCounts_htseq <- merged_data
time_points <- c("0", "1", "2", "8", "12", "24", "48")
selected_rows <- which(combined_meta$time %in% time_points)
RawCounts_htseq <- RawCounts_htseq[, selected_rows, drop = FALSE]
combined_meta_subset <- combined_meta[selected_rows, , drop = FALSE]

# Print the dimensions of the subsets to verify
dim(RawCounts_htseq)
dim(combined_meta_subset)
sample.info<-combined_meta_subset
sample.info$time<-as.factor(sample.info$time)
sample.info$dataset<-as.factor(sample.info$dataset)
Raw_Count_Plot<-ggplot(RawCounts_htseq) +
  geom_histogram(aes(x = GSM3882337), stat = "bin", bins = 200) +
  xlim(-5, 500)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")+
  ggtitle("Raw Count Distribution")
Raw_Count_Plot
mean_counts <- apply(RawCounts_htseq, 1, mean)
variance_counts <- apply(RawCounts_htseq, 1, var)
df <- data.frame(mean_counts, variance_counts)

MVR<-ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()+
  ggtitle("Mean-Variance Relation of the Raw Counts")
MVR

#Density Plot of the log-transformed Raw Count data 
df_log2<-log2(RawCounts_htseq)
df_log2_melted<-melt(df_log2)
Raw_Count_Density_Plot<-ggplot(df_log2_melted, aes(x=value, fill=variable))+geom_density(alpha=.25)
Raw_Count_Density_Plot+labs(x='Log2 Counts',y='Density',title='Raw Count Distribution',fill='Sample')

myCPM <- cpm(RawCounts_htseq)

#plot(myCPM[,1],RawCounts_htseq[,1])
plot(myCPM[,1],RawCounts_htseq[,1],ylim=c(0,50),xlim=c(0,3),xlab="CPM", ylab="Raw Counts")
abline(v=0.37)
abline(h=10)
#filtering low count genes
thresh <- myCPM > 0.37
Keep <- as.vector(rowSums(thresh) == ncol(RawCounts_htseq) )
filtered_counts<-RawCounts_htseq[Keep,]
dim(filtered_counts)

#Plotting counts after filtering low counts
Filtered_Count_Plot<-ggplot(filtered_counts) +
  geom_histogram(aes(x =GSM3882337), stat = "bin", bins = 200) +
  xlim(-5, 1000)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")
Filtered_Count_Plot

#Plotting Sample Density Distributions after low count filtering
df_log2<-log2(filtered_counts)
df_log2_melted<-melt(df_log2)
Filtered_Count_Density_Plot<-ggplot(df_log2_melted, aes(x=value, fill=variable))+geom_density(alpha=.25)
Filtered_Count_Density_Plot+labs(x='Log2 Intensity',y='Density',title='After Filtering',fill='Sample')

dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = sample.info,
                              design= ~ time+dataset)
#dds <- DESeq(dds)
vst <- assay(vst(dds))
p <- pca(vst, metadata = sample.info)
biplot(p,showLoadings = F, colby = 'dataset',lab = NULL, legendPosition = "right", legendLabSize = 10)
biplot(p,showLoadings = F, colby = 'time',lab = NULL, legendPosition = "right", legendLabSize = 10)
eigencorplot(p,
             components = getComponents(p, 1:8),  # Adjust components as needed
             metavars = c('time', 'dataset', 'cell_type'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2, 
             fontCorval = 2, 
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE, 
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))

pairsplot(p,
          components = getComponents(p, c(1:6)),
          hline = 0, vline = 0,
          pointSize = 0.8,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'dataset',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

adjusted_counts <- ComBat_seq(as.matrix(filtered_counts), batch=sample.info$dataset)
adjusted_counts<-as.data.frame(adjusted_counts)
#Plotting Sample Density Distributions after batch correction
df_log2<-log2(adjusted_counts)
df_log2_melted<-melt(df_log2)
Adjusted_Count_Density_Plot<-ggplot(df_log2_melted, aes(x=value, fill=variable))+geom_density(alpha=.25)
Adjusted_Count_Density_Plot+labs(x='Log2 Intensity',y='Density',title='After Batch Correction',fill='Sample')

dds_adj <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                  colData = sample.info,
                                  design= ~ time+dataset)
#dds <- DESeq(dds)
vst <- assay(vst(dds_adj))
p_adj <- pca(vst, metadata = sample.info)
biplot(p_adj,showLoadings = F, colby = 'dataset',lab = NULL, legendPosition = "right", legendLabSize = 10)
biplot(p_adj,showLoadings = F, colby = 'cell_type', lab = NULL, legendPosition = "right", legendLabSize = 10)
# Compute pairwise correlation values
dds_lrt <- DESeq(dds_adj, test="LRT", reduced = ~ dataset)
# Extract results
res_LRT <- results(dds_lrt)
# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.01
sigLRT_genes <- res_LRT_tb %>% dplyr::filter(padj < 0.01)

cluster_vst <- vst[sigLRT_genes$gene, ]
sample.info$time <-factor(sample.info$time)
clusters <- degPatterns(cluster_vst, metadata = sample.info, time = "time", col = "dataset")
ggsave("degPatterns_Plot.png", plot = clusters$plot, width = 20, height = 18, dpi = 300)
save(clusters, cluster_vst, file = "clusters_col.RData")


# Define the fix function
fix <- function(x) {
  x$gene <- rownames(x)
  x$genes <- NULL
  rownames(x) <- NULL
  return(x)
}
rownames(clusters$df) <- sub("^X", "", rownames(clusters$df))
# Initialize an empty list to store the groups
group_list <- list()

for (i in 1:50) {
  group_list[[paste0("group", i)]] <- clusters$df %>%
    dplyr::filter(cluster == i)
  group_list[[paste0("group", i)]] <- fix(group_list[[paste0("group", i)]])
}
group_list <- group_list[sapply(group_list, function(x) nrow(x) > 0)]
for (name in names(group_list)) {
  cluster_num <- unique(group_list[[name]]$cluster)  # Get the cluster number
  file_name <- paste0("cluster", cluster_num, ".csv")  # Create the correct file name
  write.csv(group_list[[name]], file_name, row.names = FALSE)  # Save as CSV
}



# Initialize an empty list to store enrichment results
enrichment_results <- list()

# Loop through the group_list and perform pathway enrichment
for (name in names(group_list)) {
  genes <- group_list[[name]]$gene
  cluster_num <- unique(group_list[[name]]$cluster)  # Get the actual cluster number
  
  enrichment <- enrichPathway(gene = genes,
                              organism = "human",
                              pvalueCutoff = 0.01,
                              pAdjustMethod = "BH",
                              readable = TRUE)
  
  if (nrow(enrichment@result) > 0) {
    enrichment_results[[paste0("cluster", cluster_num)]] <- enrichment
    # Write the enrichment results to a CSV file using the correct cluster number
    write.csv(enrichment@result, file = paste0("enrichment_results_cluster", cluster_num, ".csv"), row.names = FALSE)
  }
}

# Create dotplots for each cluster with significant results
for (name in names(enrichment_results)) {
  if (length(enrichment_results[[name]]) > 0 && nrow(enrichment_results[[name]]@result) > 0) {
    # Filter for significant results
    significant_results <- enrichment_results[[name]]@result %>% 
      filter(p.adjust < 0.01)
    
    if (nrow(significant_results) > 0) {
      p <- dotplot(enrichment_results[[name]], showCategory = 12, label_format = 70) +
        ggtitle(paste("Enrichment Dotplot for", name))
      
      # Save dotplot for the corresponding cluster
      ggsave(paste0("enrichment_dotplot_", name, ".png"), plot = p, width = 12, height = 7, dpi = 300)
    } else {
      message(paste("No significant pathways (p.adjust < 0.01) found for", name))
    }
  } else {
    message(paste("No enrichment results found for", name))
  }
}

# Create a combined dotplot for all clusters
all_results <- do.call(rbind, lapply(names(enrichment_results), function(name) {
  result <- enrichment_results[[name]]@result
  if (nrow(result) > 0) {
    significant_results <- result %>% filter(p.adjust < 0.01)
    if (nrow(significant_results) > 0) {
      significant_results$Cluster <- name  # Use the actual cluster name
      return(significant_results)
    }
  }
  return(NULL)
}))

# If significant pathways exist, create a combined dotplot
if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results <- all_results[order(all_results$p.adjust), ]
  top_pathways <- all_results %>%
    group_by(Cluster) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  p_combined <- ggplot(top_pathways, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cluster", y = "Pathway", size = "Gene Count", color = "Adjusted p-value") +
    ggtitle("Top 5 Enriched Pathways Across All Clusters")
  
  ggsave("combined_enrichment_dotplot.png", plot = p_combined, width = 20, height = 12, dpi = 300)
} else {
  message("No significant enrichment results (p.adjust < 0.01) found across all clusters")
}


#subset mitochondrial genes
mito_genes <- read.table("C:/Users/Aditi/Documents/Desktop/R/DGE Analysis Main/data/HumanMitoCarta.txt", header = TRUE, sep = "\t")
entrez_ids <- mito_genes$HumanGeneID
sigLRT_genes_mit <- sigLRT_genes %>% filter(sigLRT_genes$gene %in% entrez_ids)
cluster_vst <- vst[sigLRT_genes_mit$gene, ]
clusters <- degPatterns(cluster_vst, metadata = sample.info, time = "time", col = "cell_type") 
# Save the plot generated by degPatterns
save(clusters, cluster_vst, file = "mito_clusters.RData")

ggsave("degPatterns_Plot_mit.png", plot = clusters$plot, width = 20, height = 18, dpi = 300)

# Define the fix function 
fix <- function(x) {
  x$gene <- rownames(x)
  x$genes <- NULL
  rownames(x) <- NULL
  return(x)
}
rownames(clusters$df) <- sub("^X", "", rownames(clusters$df))
# Initialize an empty list to store the groups
group_list <- list()

for (i in 1:30) {
  group_list[[paste0("group", i)]] <- clusters$df %>%
    dplyr::filter(cluster == i)
  group_list[[paste0("group", i)]] <- fix(group_list[[paste0("group", i)]])
}
group_list <- group_list[sapply(group_list, function(x) nrow(x) > 0)]
for (name in names(group_list)) {
  cluster_num <- unique(group_list[[name]]$cluster)  # Get the cluster number
  file_name <- paste0("cluster_mit", cluster_num, ".csv")  # Create the correct file name
  write.csv(group_list[[name]], file_name, row.names = FALSE)  # Save as CSV
}
# Initialize an empty list to store enrichment results
enrichment_results <- list()

# Perform pathway enrichment and save results
for (name in names(group_list)) {
  genes <- group_list[[name]]$gene
  cluster_num <- unique(group_list[[name]]$cluster)  # Get actual cluster number
  
  enrichment <- enrichPathway(gene = genes,
                              organism = "human",
                              pvalueCutoff = 0.01,
                              pAdjustMethod = "BH",
                              readable = TRUE)
  
  if (nrow(enrichment@result) > 0) {
    enrichment_results[[paste0("cluster", cluster_num)]] <- enrichment
    write.csv(enrichment@result, file = paste0("enrichment_results_cluster", cluster_num, ".csv"), row.names = FALSE)
  }
}

# Create dotplots for each cluster
for (name in names(enrichment_results)) {
  significant_results <- enrichment_results[[name]]@result %>% filter(p.adjust < 0.01)
  
  if (nrow(significant_results) > 0) {
    p <- dotplot(enrichment_results[[name]], showCategory = 12, label_format = 70) +
      ggtitle(paste("Enrichment Dotplot for", name))
    ggsave(paste0("enrichment_dotplot_", name, ".png"), plot = p, width = 12, height = 7, dpi = 300)
  } else {
    message(paste("No significant pathways (p.adjust < 0.01) found for", name))
  }
}

# Create a combined dotplot for all clusters
all_results <- do.call(rbind, lapply(names(enrichment_results), function(name) {
  enrichment_results[[name]]@result %>% filter(p.adjust < 0.01) %>% 
    mutate(Cluster = name)
}))

if (!is.null(all_results) && nrow(all_results) > 0) {
  top_pathways <- all_results %>% group_by(Cluster) %>% slice_head(n = 5)
  
  p_combined <- ggplot(top_pathways, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cluster", y = "Pathway", size = "Gene Count", color = "Adjusted p-value") +
    ggtitle("Top 5 Enriched Pathways Across All Clusters")
  
  ggsave("combined_enrichment_dotplot.png", plot = p_combined, width = 20, height = 12, dpi = 300)
} else {
  message("No significant enrichment results (p.adjust < 0.01) found across all clusters")
}












