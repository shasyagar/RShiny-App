## Huntington's Disease RNA-Seq Data: Insights and Analysis
### Final Project for BF591 (R for Biological Sciences)

This project involved creating an RShiny application to facilitate the exploration of a differential expression dataset derived from post-mortem Huntington’s Disease brain samples. The data is publicly available via GEO, and the project references findings published in this paper. The app's design is customized with CSS, bypassing the default RShiny themes for a tailored aesthetic.



<img width="1799" alt="Screenshot 2024-12-16 at 01 53 17" src="https://github.com/user-attachments/assets/635addca-2db0-4b23-856b-091c3df968f4" />



Project Objectives and Dataset

The app provides interactive tools for analyzing a dataset related to Huntington’s Disease. Various input files support the app’s functionality, including:

Metadata: Information about samples/experiments (GSE64810_series_matrix.txt).

Normalized Counts: Processed expression data (GSE64810_mlhd_DESeq2_norm_counts_adjust.csv).

Differential Expression: Outlier-trimmed analysis results (GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv).

Gene Conversion: File mapping Ensembl IDs to gene symbols (human_id2gene.txt).

Pathway Analysis Results: FGSEA outcomes (fgsea_results.csv), generated using a hallmark gene set (c2.all.v2023.1.Hs.symbols.gmt).


### App Features and Tabs

#### 1. Sample Metadata
This section allows users to explore metadata with filterable columns, a search bar, and dynamic visualizations:

Summary Tab: Displays column properties, including data types and calculated statistics for numerical fields.

Metadata Tab: Provides an interactive view of the complete metadata file with horizontal scrolling.

Violin Plot Tab: Generates violin plots based on user-selected continuous variables via dropdowns.

<img width="1799" alt="Screenshot 2024-12-16 at 01 54 58" src="https://github.com/user-attachments/assets/76885c26-506c-4a59-8463-8b0a24b9ec8f" />


<img width="1799" alt="Screenshot 2024-12-16 at 01 55 09" src="https://github.com/user-attachments/assets/0353b99d-0884-44b8-bc02-31db43dd61be" />



#### 2. Counts Matrix
This tab facilitates exploration of normalized counts data through:

Filtered Data: A summary of genes passing/failing user-defined variance thresholds.

Scatter Plot: Visualizes Ranked Median vs Log Variance and Ranked Median vs Number of Zeros, with color-coded data points.

Heatmap: Displays expression patterns across samples.

PCA: Visualizes principal components with adjustable settings for the number of PCs.


<img width="1799" alt="Screenshot 2024-12-16 at 01 57 17" src="https://github.com/user-attachments/assets/cc013d9f-a442-4d75-a933-a24de92c9ba6" />

<img width="1799" alt="Screenshot 2024-12-16 at 01 57 26" src="https://github.com/user-attachments/assets/8b8ba208-39d4-4bfb-9bd1-41f7368b39b0" />

<img width="1799" alt="Screenshot 2024-12-16 at 01 57 35" src="https://github.com/user-attachments/assets/83002e3f-89f1-4977-bfde-0e887a9f7b4b" />


#### 3. Differential Expression
Users can interactively explore differentially expressed genes with:

Metadata Tab: Full data display with filtering and search capabilities.

Plot Tab Includes:

  Volcano Plot: Customizable padj thresholds, axis selection, and dynamic color settings.

Filtered Table: Displays genes passing the user-defined filters.

<img width="1799" alt="Screenshot 2024-12-16 at 02 14 27" src="https://github.com/user-attachments/assets/badc56a1-5031-43bb-8366-80352abb6880" />

<img width="1799" alt="Screenshot 2024-12-16 at 02 14 47" src="https://github.com/user-attachments/assets/685196e0-c2d5-46d7-aabb-554a3cf3f36b" />


#### 4. GSEA (Gene Set Enrichment Analysis)
This section highlights FGSEA results:

Barplot Tab: Displays pathway enrichment scores filtered by user-defined thresholds.

Table Tab: Lists pathways passing the threshold, with options for pathway type and downloadable results.

Scatter Plot Tab: Plots normalized enrichment scores (NES) against -log10(padj) values.

<img width="1799" alt="Screenshot 2024-12-16 at 02 17 53" src="https://github.com/user-attachments/assets/71f28e65-87d5-4a2d-ab9c-eeb72300ec25" />

<img width="1799" alt="Screenshot 2024-12-16 at 02 18 18" src="https://github.com/user-attachments/assets/e4f9d5e7-6c5e-4d33-bb75-b1e04ae34034" />

<img width="1799" alt="Screenshot 2024-12-16 at 02 19 33" src="https://github.com/user-attachments/assets/24bad97b-5836-4c03-89ea-5c5e07260a2e" />

