# ==============================================================================
# PRJNA794860 Complete RNA-Seq Analysis & Plotting Pipeline
# DESeq2, PCA, Volcano Plot, Heatmap, and Functional Enrichment (KEGG & GO)
# ==============================================================================

# 1. Load Required Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Set clean publication theme
theme_set(theme_classic(base_size = 12) + 
            theme(
              axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
            ))

# ==============================================================================
# 2. Load Count Matrix and Setup Metadata
# ==============================================================================
# Ensure 'Sample.csv' is in your working directory
count_data <- read.csv('Sample.csv', header = TRUE, row.names = 1)

# Filter low count genes across all 6 replicates
sample_cols <- c('PCOS_CT_R1', 'PCOS_CT_R2', 'PCOS_CT_R3', 'PCOS_DI_R1', 'PCOS_DI_R2', 'PCOS_DI_R3')
count_data <- count_data[rowSums(count_data[, sample_cols]) > 50, ]

# Setup sample metadata / condition factors
condition <- factor(c("PCOS_CT", "PCOS_CT", "PCOS_CT", "PCOS_DI", "PCOS_DI", "PCOS_DI"))
coldata <- data.frame(rownames = colnames(count_data), condition = condition)
rownames(coldata) <- coldata$rownames

# ==============================================================================
# 3. Run DESeq2 Differential Expression
# ==============================================================================
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("PCOS_CT", "PCOS_DI"))

dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result <- na.omit(deseq_result)

# Export complete DESeq2 results to CSV
deseq_df <- as.data.frame(deseq_result)
write.csv(deseq_df, file = "deseq_PCOS_CT_vs_DI_results.csv", row.names = TRUE)
cat("DESeq2 results successfully exported to 'deseq_PCOS_CT_vs_DI_results.csv'.\n")

# Variance Stabilizing Transformation for downstream plots
vsd <- vst(dds, blind = FALSE)

# ==============================================================================
# 4. Panel A: PCA Plot of Samples
# ==============================================================================
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("PCOS_CT" = "#E05A47", "PCOS_DI" = "#3F72AF")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  theme(legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14))

ggsave("Panel_A_PCA.png", pca_plot, width = 6, height = 5, dpi = 300)

# ==============================================================================
# 5. Panel B: Volcano Plot
# ==============================================================================
res_df <- as.data.frame(deseq_result)
res_df$diffexpressed <- "NS"
res_df$diffexpressed[res_df$log2FoldChange > 1.0 & res_df$padj < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1.0 & res_df$padj < 0.05] <- "DOWN"
res_df$diffexpressed <- factor(res_df$diffexpressed, levels = c("DOWN", "NS", "UP"))

top_genes <- res_df %>% 
  rownames_to_column("Gene") %>% 
  filter(diffexpressed != "NS") %>% 
  arrange(padj) %>% 
  head(15)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "#3F72AF", "NS" = "grey70", "UP" = "#E05A47")) +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed", linewidth = 0.6) +
  labs(title = "Volcano Plot", x = "log₂ Fold Change", y = "-log₁₀ Adjusted p-value") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14)) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, 
                  box.padding = 0.5, max.overlaps = Inf, show.legend = FALSE)

ggsave("Panel_B_Volcano.png", volcano_plot, width = 7, height = 5, dpi = 300)

# ==============================================================================
# 6. Panel C: Clustered Heatmap of Top DE Genes
# ==============================================================================
filtered_results <- deseq_result[deseq_result$padj < 0.05 & abs(deseq_result$log2FoldChange) > 1, ]
top10Upregulated <- head(filtered_results[order(filtered_results$log2FoldChange, decreasing = TRUE), ], 10)
top10Downregulated <- head(filtered_results[order(filtered_results$log2FoldChange, decreasing = FALSE), ], 10)
top20 <- rbind(top10Upregulated, top10Downregulated)

top20_vst <- assay(vsd)[rownames(top20), ]
col_palette <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(400)

group <- data.frame(condition = colData$condition)
rownames(group) <- colnames(top20_vst)
anno_colors <- list(condition = c("PCOS_CT" = "#E05A47", "PCOS_DI" = "#3F72AF"))

png("Panel_C_Heatmap.png", width = 2000, height = 2400, res = 300)
pheatmap(top20_vst, 
         scale = 'row',
         fontsize_row = 8,
         fontsize_col = 8,
         color = col_palette,
         annotation_col = group,
         annotation_colors = anno_colors,
         cluster_rows = TRUE,
         cutree_rows = 2,
         cluster_cols = FALSE,
         main = "Top 20 Differentially Expressed Genes Heatmap")
dev.off()

# ==============================================================================
# 7. Panels D–G: KEGG and GO Functional Enrichment (clusterProfiler)
# ==============================================================================
sig_genes <- rownames(filtered_results)
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if(nrow(gene_ids) > 0) {
  # D. KEGG Pathway Enrichment
  kegg_res <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
  if(!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    p_kegg <- barplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment") +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave("Panel_D_KEGG.png", p_kegg, width = 7, height = 6, dpi = 300)
  }
  
  # E. GO Biological Process (BP)
  go_bp <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
  if(!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    p_bp <- barplot(go_bp, showCategory = 10, title = "GO Enrichment - BP") +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave("Panel_E_GO_BP.png", p_bp, width = 7, height = 6, dpi = 300)
  }
  
  # F. GO Cellular Component (CC)
  go_cc <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05, readable = TRUE)
  if(!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
    p_cc <- barplot(go_cc, showCategory = 10, title = "GO Enrichment - CC") +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave("Panel_F_GO_CC.png", p_cc, width = 7, height = 6, dpi = 300)
  }
  
  # G. GO Molecular Function (MF)
  go_mf <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05, readable = TRUE)
  if(!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
    p_mf <- barplot(go_mf, showCategory = 10, title = "GO Enrichment - MF") +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave("Panel_G_GO_MF.png", p_mf, width = 7, height = 6, dpi = 300)
  }
}

cat("--- Complete PRJNA794860 Pipeline Finished Successfully! All Panels Generated. ---\n")