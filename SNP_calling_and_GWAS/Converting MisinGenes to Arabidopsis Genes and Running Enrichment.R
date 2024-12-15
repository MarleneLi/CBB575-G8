gene_map <- read.table("new_gene.csv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gene_map) <- c("MisinGene","ArabidopsisGene")

traits <- c("Al", "Ca", "Fe", "K", "Cr", "Pb", "Mn", "As", "Cu", "Ni", "Zn", "S")

for (t in traits) {
  misingene_file <- paste0(t, "_misingenes.txt")
  if (!file.exists(misingene_file)) next
  
  misingenes <- readLines(misingene_file)
  misin_df <- data.frame(MisinGene=misingenes, stringsAsFactors=FALSE)
  
  merged <- merge(misin_df, gene_map, by="MisinGene", all.x=TRUE)
  final_arabidopsis_genes <- merged$ArabidopsisGene[merged$ArabidopsisGene != "None"]
  final_arabidopsis_genes <- na.omit(final_arabidopsis_genes)
  
  output_file <- paste0(t, "_arabidopsis_genes.txt")
  writeLines(final_arabidopsis_genes, output_file)
}





library(gprofiler2)
library(cowplot)
library(stringr)
# Define your traits array
traits <- c( "Ca", "Fe", "Cr", "Pb",  "As", "Cu", "Zn", "S")
#Al,"K","Mn","Ni" do not have enough genes
#install.packages("ggplot2")
library(ggplot2)
for (t in traits) {
  # Read the arabidopsis genes for this trait
  arab_genes_file <- paste0(t, "_arabidopsis_genes.txt")
  if (!file.exists(arab_genes_file)) {
    cat("Warning:", arab_genes_file, "not found. Skipping.\n")
    next
  }
  
  arab_genes <- readLines(arab_genes_file)
  
  # Run g:Profiler enrichment
  res <- gost(query = arab_genes, organism = "athaliana", user_threshold = 0.05)
  
  if(is.null(res$result) || nrow(res$result) == 0) {
    cat("No significant terms for", t, "\n")
    next
  }
  
  # Reorder by p_value ascending
  res_by_p <- res$result[order(res$result$p_value),]
  
  # Top 10 by p_value
  top10_p <- head(res_by_p, 10)
  top10_p$logP <- -log10(top10_p$p_value)
  
  # Reorder by term_size ascending
  res_by_size <- res$result[order(res$result$term_size),]
  
  # Top 10 by term_size
  top10_size <- head(res_by_size, 10)
  top10_size$logP <- -log10(top10_size$p_value)
  cat("Top 10 terms by smallest p_value:\n")
  print(head(res_by_p$term_name, 10))
  # res_by_size <- res$result[order(res$result$term_size), ]  # order by term_size ascending
  cat("\nTop 10 terms by smallest term_size:\n")
  print(head(res_by_size$term_name, 10))
  # Plot top 10 by p-value
  # Terms as factor in order of appearance (lowest p first)
  top10_p$term_name <- factor(top10_p$term_name, levels=rev(top10_p$term_name))
  
  p1 <- ggplot(top10_p, aes(x=logP, y=term_name)) +
    geom_bar(stat="identity", fill="#2C7BB6") +
    labs(title=paste(t, "- Top 10 Terms by p-value"),
         x="-log10(p-value)", y="") +
    theme_minimal()
  p1 <- p1 + theme(
    axis.text.y = element_text(size=10),
    plot.margin = unit(c(1,10,1,1), "lines")
  )
  # Save p-value-based plot
  p_file <- paste0(t, "_top_by_p.png")
  ggsave(p_file, p1, width=6, height=5)
  
  # Plot top 10 by term_size
  top10_size$term_name <- factor(top10_size$term_name, levels=rev(top10_size$term_name))
  
  p2 <- ggplot(top10_size, aes(x=logP, y=term_name)) +
    geom_bar(stat="identity", fill="#D7191C") +
    labs(title=paste(t, "- Top 10 Terms by size"),
         x="-log10(p-value)", y="") +
    theme_minimal()
  
  # Save term_size-based plot
  s_file <- paste0(t, "_top_by_size.png")
  ggsave(s_file, p2, width=6, height=5)
  p_combined <- plot_grid(p1, p2, nrow=1, labels = c("A", "B"))
  
  # Save the combined plot
  combined_file <- paste0(t, "_combined.png")
  ggsave(combined_file, p_combined, width=16, height=6)
  cat("Plots saved for", t, "\n")
}

