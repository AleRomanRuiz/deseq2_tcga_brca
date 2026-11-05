library(DESeq2)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

counts_path <- "data/raw/counts.csv.gz"
metadata_path <- "data/raw/metadata.csv.gz"
condition_column <- "site"
plate_column <- "plate"
condition_test <- "Primary_Solid_Tumor"
condition_ref <- "Solid_Tissue_Normal"
min_count <- 10
min_samples <- 20
alpha <- 0.05
lfc_threshold <- 1

message("[1/8] Leyendo datos de entrada (TCGA-BRCA)...")
counts_df <- read.csv(counts_path, row.names=1, check.names=FALSE, stringAsFactors=FALSE) ## esto es muy importante para que no se cambien los nombres de las columnas
metadata <- read.csv(metadata_path, row.names=1, stringAsFactors=FALSE, check.names=FALSE) 

stopifnot(condition_column %in% colnames(metadata_df))
stopifnot(plate_column %in% colnames(metadata_df))

counts_mat <- as.matrix(counts_df)

message("[2/8] Alineando muestras entre counts y metadata...")
metadata <- metadata_df[colnames(counts_mat), ]

metadata[[condition_column]] <- factor (metadata [[condition_column]])
metadata[[plate_column]] <- factor (metadata [[plate_column]])
metadata[[condition_column]] <- relevel(metadata[[condition_column]], ref = condition_ref)

message ("[3/8] Creando objeto DESeq2 y filtrando genes de baja cuenta...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat, 
  colData = metadata, 
  design = as.formula(paste("~", condition_column, "+", plate_column))
)

# Gen con >= min_count lecturas en >= min_samples muestras (tamaño del grupo mas pequeño).
dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]

message ("[4/8] Ejecutando DESeq2...")
dds <- DESeq(dds)

png("results/dispersion_plot.png", width = 800, height = 600)
plotDispEsts(dds)
dev.off()

message("[5/8] Extrayendo resultados y guardando tablas...")
res <- results(
dds,
contrast = c(condition_column, condition_test, condition_ref),
alpha = alpha,
)
res_df <- as.data.frame(res)
res_df <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.table(
  res_df,
  file = "results/deseq2_a11_results.csv"
  sep = "/t"
  quote = FALSE
  row.names = FALSE
)

sig_df <- subset(res_df, !is.na(padj) & padj < alpha)
de_df <- subset(res_df, !is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc_threshold)
up_df <- subset(de_df, log2FoldChange > 0)
down_df <- subset(de_df, log2FoldChange < 0)

top10_up <- head(up_df, 10)
top10_down <- head(down_df, 10)

de_df$top10_label <- "other"
de_df$top10_label[de_df$gene_id %in% top10_up$gene_id] <- "top10_up"
de_df$top10_label[de_df$gene_id %in% top10_down$gene_id] <- "top10_down"

write.table(
  sig_df,
  file = "results/deseq2_sig_padj_0.05.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  de_df,
  file = "results/deseq2_de_padj_0.05_lfc_1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top10_up,
  file = "results/top10_up_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top10_down,
  file = "results/top10_down_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
















volcano_df <- res_df[!is.na(res_df$padj), ]
volcano_df%status <- "No significativo"
volcano_df$status[volcano_df$gene_id %in% up_df$gene_id] <- "Up"
volcano_df$status[volcano_df$gene_id %in% down_df$gene_id] <- "Down"

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_text(data = head(volcano_df, 10), aes(label = gene_id), vjust = -0.05, show.legend = FALSE) +
  scale_color_manual(values = c(Up = "red", Down = "blue", "No significativo" = "grey")) +
  theme_bw()

ggsave("results/plot_volcano.png", p_volcano, width = 7, height = 6)




png("results/plot_ma.png", width = 800, height = 600)
plotMA(res, ylim = c(-4, 4))
dev.off()

message("[6/8] Contextualizacion biologica: gen destacado en el contraste...")










message("[7/8] Enriquecimiento GO y Reactome...")
symbol_to_entrez <- bitr(
  de_df$gene_id,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
universe_entrez <- unique(symbol_to_entrez$ENTREZID)

up_map <- bitr(unique(up_df$gene_id), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_map <- bitr(unique(down_df$gene_id), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)



go_up <- enrichGO(
  gene = up_map$ENTREZID,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)
write.table(as.data.frame(go_up), "results/go_up_enrichment.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_go_bp_up.png", dotplot(go_up, showCategory = 10), width = 8, height = 6)

go.down <- enrichGO(
  gene = unique(down_map$ENTREZID),
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)
write.table(as.data.frame(go.down), "results/go_bp_down_enrichment.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_go_bp_down.png", dotplot(go.down, showCategory = 10), width = 8, height = 6)

reactome_up <- enrichPathway(
  gene = unique(up_map$ENTREZID),
  universe = universe_entrez,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)
write.table(as.data.frame(reactome_up), "results/reactome_up.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_reactome_up.png", dotplot(reactome_up, showCategory = 10), width = 8, height = 6)

reactome_down <- enrichPathway(
  gene = unique(down_map$ENTREZID),
  universe = universe_entrez,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)
write.table(as.data.frame(reactome_down), "results/reactome_down.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
ggsave("results/plot_reactome_down.png", dotplot(reactome_down, showCategory = 10), width = 8, height = 6)

message("[8/8] Guardando trazabilidad de sesion...")
writelines(capture.output(sessionInfo()), "logs/session_info.txt")

summary_lines <- c(
  paste("n_samples:", ncol(counts_mat)),
  paste("n_genes_input:", nrow(counts_mat)),
  paste("n_genes_after_filter:", nrow(dds)),
  paste("n_significant_padj_It", alpha, ":", nrow(sig_df), sep = ""),
  paste("n_de_absLog2FC_gt_", lnf_threshold, "_padj_It", alpha, ":", nrow(de_df), sep = ""),
  paste("n_up:", nrow(up_df)),
  paste("n_down:", nrow(down_df)),
  paste("top10_up_file:", "results/top10_up_genes.tsv"),
  paste("top10_down_file:", "results/top10_down_genes.tsv"),
  paste("top10_plot_file:", "results/top10_down_genes.png"),
  paste("volcano_plot_file:", "results/plot_volcano.png")
)
writeLines(summary_lines, "logs/run_summary.txt")