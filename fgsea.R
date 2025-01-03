library('tidyverse')
library('fgsea')

# File paths
id2_gene_file <- "data/human_id2gene.txt"
gene_set_file <- "data/c2.all.v2023.1.Hs.symbols.gmt"
deseq_file <- "data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt"

# Load data with proper column names and ensuring all columns are read as strings
id2gene_df <- read.table(id2_gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("gene_symbol", "ensembl_id"))
deseq_df <- read.csv(deseq_file, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

# Function to create ranked log2FC
make_ranked_log2fc <- function(labeled_results, id2gene_df) {
  merged_data <- merge(labeled_results, id2gene_df, by.x = "symbol", by.y = "ensembl_id", all.x = TRUE)
  merged_data$symbol <- as.character(merged_data$symbol)  # Ensure 'symbol' is a string
  ranked_log2fc_df <- merged_data[order(merged_data$log2FoldChange, decreasing = TRUE), ]
  # Remove duplicates by keeping the gene with the highest log2FoldChange
  ranked_log2fc_df <- ranked_log2fc_df[!duplicated(ranked_log2fc_df$symbol), ]
  
  ranked_log2fc_named <- setNames(ranked_log2fc_df$log2FoldChange, ranked_log2fc_df$symbol)
  return(ranked_log2fc_named)
}

rnk_list <- make_ranked_log2fc(deseq_df, id2gene_df)

# Clean up NAs in ranked_genes (remove entries with NA in either names or values)
rnk_list <- rnk_list[!is.na(names(rnk_list)) & !is.na(rnk_list)]

# Function to run fgsea
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # Make sure there are no NA values in the stats before running fgsea
  rnk_list <- rnk_list[is.finite(rnk_list)]
  
  gene_sets <- gmtPathways(gmt_file_path)
  fgsea_results <- fgsea(pathways = gene_sets, 
                         stats = rnk_list, 
                         minSize = min_size, 
                         maxSize = max_size)
  
  return(as_tibble(fgsea_results))
}

fgsea_result <- run_fgsea(gene_set_file, rnk_list, 15, 1115)

# Convert list column to a comma-separated string
fgsea_result$leadingEdge <- sapply(fgsea_result$leadingEdge, paste, collapse = ", ")

# Ensure all columns in fgsea_result are character strings
fgsea_result <- fgsea_result %>% mutate(across(everything(), as.character))

# Print the head of the result for inspection
print(head(fgsea_result))

# Write to CSV with strings as columns
write.csv(fgsea_result, "data/fgsea_results.csv", row.names = FALSE, quote = TRUE)
