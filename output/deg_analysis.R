# Load required libraries
library(tidyverse)
library(DESeq2)

# Load count data
counts <- read_csv("gene_count_matrix.csv", col_names = TRUE)

# Prepare metadata
metadata <- tibble(
  SampleID = c("1BS", "2BS", "3BS", "1DS", "2DS", "3DS"),
  Treatment = c("B", "B", "B", "D", "D", "D"),
  Replicate = c(1, 2, 3, 1, 2, 3),
  Season = rep("S", 6)
)
metadata <- column_to_rownames(metadata, var = "SampleID")
rownames(counts) <- counts$gene_id  # Set gene_id as row names
counts <- counts[,-1]  # Remove gene_id column from the count matrix
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]

stopifnot(all(rownames(metadata) == colnames(counts)))  # Stop if there's a mismatch

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Treatment
)
dds <- DESeq(dds)

# Perform differential expression analysis (B vs D)
res <- results(dds, contrast = c("Treatment", "B", "D"))

# Filter significant results (padj < 0.05)
res_sig <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  filter(padj < 0.05 & !is.na(padj))

# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(normalized_counts) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "Sample", values_to = "Count") %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample") %>%
  filter(gene_id %in% res_sig$gene_id)


# Volcano plot using ggplot2
volcano_data <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(Significant = ifelse(padj < 0.05 & !is.na(padj), "Significant", "Not Significant"))

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_classic()

ggsave("volcano_plot_B_vs_D.png", width = 8, height = 6, dpi = 300)
