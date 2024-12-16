#### R code for GWENA analysis and data visualization
#### Written by Pei-Wei Sun, 2024 

#### Part1: Read expression data 
#### Part2: Read phenotype data
#### Part3: Data processing and module construction 
#### Part4: Module information visualization
#### Part5: Phenotype association and visualization
#### Part6: Venn diagram, Upset plot, heatmap, and PCA analysis


##### Read DEG data for each tissue
leaf=read.csv("project/raw_data/dge_lists/leaf.csv",header = T)
leaf_sig=subset(leaf, (log2FoldChange > 2 | log2FoldChange < -2) & padj < 0.05)
leaf_sig_up=subset(leaf, (log2FoldChange > 2 ) & padj < 0.05)
leaf_sig_down=subset(leaf, (log2FoldChange < -2) & padj < 0.05)

stem=read.csv("project/raw_data/dge_lists/stem.csv",header = T)
stem_sig=subset(stem, (log2FoldChange > 2 | log2FoldChange < -2) & padj < 0.05)
stem_sig_up=subset(stem, (log2FoldChange > 2 ) & padj < 0.05)
stem_sig_down=subset(stem, (log2FoldChange < -2) & padj < 0.05)

root=read.csv("project/raw_data/dge_lists/root.csv",header = T)
root_sig=subset(root, (log2FoldChange > 2 | log2FoldChange < -2) & padj < 0.05)
root_sig_up=subset(root, (log2FoldChange > 2 ) & padj < 0.05)
root_sig_down=subset(root, (log2FoldChange < -2) & padj < 0.05)

# Generate the list of common genes based on the "gene" column
common_genes <- Reduce(intersect, list(leaf_sig$gene, root_sig$gene, stem_sig$gene))
###try work on specific DEGs for root, stem, or leaf
common_genes <- root_sig$gene
common_genes <- stem_sig$gene
common_genes <- leaf_sig$gene


##### Do the GWENA analysis
#### Part1: Read expression data 
BiocManager::install("Kumquatum/GWENA")
library(GWENA)
library(magrittr)
#threads_to_use <- 3
# Import expression table
# read gene expression data:
kuehne_expr = as.data.frame(t(read.csv("project/raw_data/counts/tpm_matrix.csv", header=TRUE, row.names=1)))
# xtract the first part before the first underscore
rownames(kuehne_expr) <- sub("^X", "", sub("_.*", "", rownames(kuehne_expr)))

#Add an underscore between the number and the letter
rownames(kuehne_expr) <- gsub("([0-9]+)([A-Za-z])", "\\1_\\2", rownames(kuehne_expr))

# Filter rows in data to keep only those where the gene column matches common_genes
kuehne_expr <- kuehne_expr[, colnames(kuehne_expr) %in% common_genes, drop = FALSE]
###only work with root, leaf, or stem once a time
# Filter rows (samples) where row names contain "_R", "_L", or "_S"
kuehne_expr <- kuehne_expr[grepl("_R", rownames(kuehne_expr)), ]
kuehne_expr <- kuehne_expr[grepl("_L", rownames(kuehne_expr)), ]
kuehne_expr <- kuehne_expr[grepl("_S", rownames(kuehne_expr)), ]

# Number of genes
ncol(kuehne_expr)
# Number of samples
nrow(kuehne_expr)

# Overview of expression table
kuehne_expr[1:5,1:5]


####part2: Read phenotype data 
# Import phenotype table: 
kuehne_traits = read.csv("project/raw_data/soil_properties_GWENA.csv", header=TRUE, row.names=1)
# Order the rows by the alphabetical order of row names
kuehne_traits <- kuehne_traits[order(rownames(kuehne_traits)), ]

###only work with root, leaf, or stem once a time
kuehne_traits <- kuehne_traits[kuehne_traits$tissue == "root", ]
kuehne_traits <- kuehne_traits[kuehne_traits$tissue == "leaf", ]
kuehne_traits <- kuehne_traits[kuehne_traits$tissue == "stem", ]

# Phenotype
unique(kuehne_traits$tissue)
# Overview of traits table
kuehne_traits[1:5,]

#### Part3: Data processing and module construction 
### Discard genes below the threshold
kuehne_expr_filtered <- filter_low_var(kuehne_expr, pct = 0.7, type = "median")

# Remaining number of genes
ncol(kuehne_expr_filtered)

# Order the rows by the alphabetical order of row names
kuehne_expr_filtered <- kuehne_expr_filtered[order(rownames(kuehne_expr_filtered)), ]


#### Replace original gene name to A. thaliana gene name
# Load the CSV file with old and new gene names
gene_mapping <- read.table("new_gene.csv", stringsAsFactors = FALSE, header = TRUE)

# Add ".v7.1" to the end of the first column
gene_mapping[, 1] <- paste0(gene_mapping[, 1], ".v7.1")

# Ensure the columns are correctly named or refer to them directly
# Assuming the first column is "old_name" and the second is "new_name"
old_names <- gene_mapping[, 1]
new_names <- gene_mapping[, 2]

# Match column names in kuehne_expr_filtered with old_names and replace with new_names
colnames(kuehne_expr_filtered) <- new_names[match(colnames(kuehne_expr_filtered), old_names)]
write.table(colnames(kuehne_expr_filtered),file="root_gene_name",row.names = FALSE)
# Remove columns with NA or "None" in the column names
kuehne_expr_filtered <- kuehne_expr_filtered[, !is.na(colnames(kuehne_expr_filtered)) & colnames(kuehne_expr_filtered) != "None"]

# Remove columns with duplicated base names (ignore .1, .2, etc.)
unique_cols <- !duplicated(gsub("\\.\\d+$", "", colnames(kuehne_expr_filtered)))
kuehne_expr_filtered <- kuehne_expr_filtered[, unique_cols]


# Construct network using annotated genes
net <- build_net(kuehne_expr_filtered, cor_func = "spearman", 
                 n_threads = 3)

# Power selected :
net$metadata$power

# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

# Module construction
modules <- detect_modules(kuehne_expr_filtered, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)
# Number of modules before merging :
length(unique(modules$modules_premerge))

# Number of modules after merging: 
length(unique(modules$modules))


#### Part4: Module information visualization 
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)
# Save the plot to a PDF
pdf("root_module.pdf", width = 7, height = 5)

# Create the bar plot with numbers on top of each bar
ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) +
  ggplot2::stat_count() +
  ggplot2::geom_text(
    stat = "count", 
    ggplot2::aes(label = ..count..), 
    vjust = -0.5,  # Position the labels just above the bars
    size = 4       # Adjust the size of the labels
  ) +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module") +
  ggplot2::theme_bw()  # Optional: Clean and minimal theme

dev.off()
# Plot expression profile
plot_expression_profiles(kuehne_expr_filtered, modules$modules)


# Enrichment (not work here)
enrichment <- bio_enrich(modules$modules,organism = "athaliana")
plot_enrichment(enrichment)


#### Part5: Phenotype association and visualization
phenotype_association <- associate_phenotype(
  modules$modules_eigengenes, 
  kuehne_traits[-1])
pdf("root_association.pdf", width = 5, height = 3)
plot_modules_phenotype(phenotype_association)
dev.off()

# Plot network of each module
# Change number in modules$modules$`1` to any module wanted to plot 
module_example <- modules$modules$`1`
graph <- build_graph_from_sq_mat(net$network[module_example, module_example])

layout_mod_1 <- plot_module(graph, upper_weight_th = 0.999, 
                            vertex.label.cex = 0, 
                            node_scaling_max = 13, 
                            legend_cex = 1,zoom=2)

# Plot the subclusters within selected module
net_mod_1 <- net$network[modules$modules$`1`, modules$modules$`1`] 
sub_clusters <- get_sub_clusters(net_mod_1)
pdf("leaf_network_m1.pdf", width = 15, height = 10)
layout_mod_1_sub_clust <- plot_module(graph, upper_weight_th = 0.9995,
                                      groups = sub_clusters,
                                      vertex.label.cex = 0, 
                                      node_scaling_max = 7, 
                                      legend_cex = 1,zoom=1)
dev.off()

# Save gene name for each module as csv
# Convert the modules to a list
module_list <- as.list(modules$modules)

# Convert the list into a data frame with each module as a column
module_df <- data.frame(do.call(cbind, lapply(module_list, function(x) {
  length_diff <- max(lengths(module_list)) - length(x)
  c(x, rep(NA, length_diff))  # Add NAs to make all columns the same length
})))

# Set column names to module names
colnames(module_df) <- names(module_list)

# Save the data frame as a CSV file
write.csv(module_df, "leaf_modules.csv", row.names = FALSE)


####============================================


### Part6: Venn diagram, Upset plot, heatmap, and PCA analysis
# Plot the Venn diagram
library(VennDiagram)
gene_list <- list(
  "leaf" = leaf_sig$gene,
  "stem" = stem_sig$gene,
  "root" = root_sig$gene
)
gene_list <- list(
  "leaf up" = leaf_sig_up$gene,
  "leaf down" = leaf_sig_down$gene,
  "stem up" = stem_sig_up$gene,
  "stem down" = stem_sig_down$gene,
  "root up" = root_sig_up$gene,
  "root down" = root_sig_down$gene
)
library(ggplot2)
library(ggVennDiagram)
pdf("venn_plot.pdf", width = 5, height = 5)
ggVennDiagram(gene_list, label_alpha = 0)+ scale_fill_gradient(low="grey90",high = "red")
dev.off()

# Plot the Upset plot
# Create a unique vector of all genes
all_genes <- unique(unlist(gene_list))
# Create a data frame where each column corresponds to a set
gene_matrix <- as.data.frame(sapply(gene_list, function(set) {
  as.integer(all_genes %in% set)
}))
rownames(gene_matrix) <- all_genes

# Convert rownames into a column to prepare for UpSetR
gene_matrix <- cbind(gene = rownames(gene_matrix), gene_matrix)
library(UpSetR)
pdf("upset_plot.pdf", width = 10, height = 7)
# Create the upset plot
upset(gene_matrix,
      sets = colnames(gene_matrix)[-1], # Exclude the 'gene' column
      sets.bar.color = "skyblue",         # Color for the set bars
      order.by = "freq",               # Order intersections by frequency
      main.bar.color = "salmon",          # Color for the intersection bars
      matrix.color = "black",          # Color for the intersection matrix
      text.scale = 1.5)                # Scale the text for readability
dev.off()





# Plot Heatmap with top 1,000 columns with the highest variance
library(pheatmap)

# Create an annotation for rows based on row names in kuehne_expr_filtered
row_annotation <- data.frame(
# Assign groups based on row names
Group <- ifelse(grepl("_L", rownames(kuehne_expr)), "leaf",
                  ifelse(grepl("_S", rownames(kuehne_expr)), "stem","root")))
rownames(row_annotation) <- rownames(kuehne_expr)  # Match annotation to row names in kuehne_expr_filtered
colnames(row_annotation)<- "Group"

# Calculate column variance
column_variance <- apply(kuehne_expr, 2, var)

# Select the top 1,000 columns with the highest variance
top_columns <- names(sort(column_variance, decreasing = TRUE))[1:1000]

# Subset the data frame based on the top columns
kuehne_expr_top <- kuehne_expr[, top_columns]

# Normalize each column to the range [0, 1]
kuehne_expr_filtered_normalized <- as.data.frame(
  apply(kuehne_expr_top, 2, function(x) (x - min(x)) / (max(x) - min(x)))
)

# Create the heatmap with pheatmap
pdf("heatmap_100variance.pdf", width = 10, height = 10)
pheatmap(
  kuehne_expr_filtered_normalized, 
  fontsize_row = 6,
  annotation_row = row_annotation,    # Add row annotations
  show_rownames = TRUE,               # Show row names
  show_colnames = FALSE,              # Hide column names
  cluster_rows = TRUE,                # Cluster rows
  cluster_cols = TRUE,                # Cluster columns
  main = ""  # Add a title
)
dev.off()

# PCA analysis 
# Perform PCA on the normalized expression data
pca_result <- prcomp(kuehne_expr_filtered_normalized, center = TRUE, scale. = TRUE)

# Create a data frame for PCA results
pca_data <- as.data.frame(pca_result$x)

# Add group information for coloring
pca_data$Group <- row_annotation$Group

# Load ggplot2 for PCA visualization
library(ggplot2)

# Plot the PCA
pdf("pca_top1000variancegene.pdf", width = 6, height = 4)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "",
      x = "PC1(33%)",
    y = "PC2(13%)"
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
dev.off()

# Perform PCA on the soil data
# Exclude the "tissue" column and select only numerical data
numerical_data <- kuehne_traits[c(1:20), -1]  # Remove the "tissue" column

# Perform PCA on the numerical columns
pca_result <- prcomp(numerical_data, center = TRUE, scale. = TRUE)

# Extract PCA scores (samples) and loadings (variables)
pca_scores <- as.data.frame(pca_result$x)  # Principal component scores for samples
pca_loadings <- as.data.frame(pca_result$rotation)  # Loadings (arrows)


# Load ggplot2 for visualization
library(ggplot2)

# PCA plot with arrows
pdf("pca_soil_with_arrows.pdf", width = 8, height = 6)
ggplot(data = pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.8, color = "gray") +  # Plot sample points in gray
  geom_text(aes(label = rownames(pca_scores)), size = 3, vjust = -0.5, hjust = 0.5) +  # Add sample names
  geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1 * 8, yend = PC2 * 8),
               arrow = arrow(length = unit(0.2, "cm")), color = "tomato", size = 0.5) +  # Plot arrows
  geom_text(data = pca_loadings, aes(x = PC1 * 8, y = PC2 * 8, label = rownames(pca_loadings)),
            color = "tomato", size = 5, vjust = -0.5) +  # Label arrows with variable names
  labs(
    title = "",
    x = "PC1(31%)",
    y = "Pc2(21%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",  # Remove legend (since points are a single color)
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
dev.off()

