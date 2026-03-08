# PROSES PENGOLAHAN DATA GEO DENGAN ANALISIS DEG
## Judul       = Data Expression in Alveolar Macrophages Induced by 
##               Lipopolysaccharide in Humans
## Sitasi      = Reynier F, de Vos AF, Hoogerwerf JJ, Bresser P et al. Gene 
##               expression profiles in alveolar macrophages induced by 
##               lipopolysaccharide in humans. Mol Med 2012 Dec 6;18(1):1303-11.
## Sampel Data = GSE40885

# TAHAP 0: Mempersiapkan Working Directory dan Package
setwd("/home/disspear/Documents/BRSP/Capstone Project")
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)

# TAHAP 1: Pengambilan Data Sampel dari GEO NCBI
gset <- getGEO("GSE40885", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]] # Download database lewat getGEO()
ex <- exprs(gset) # Matriks ekspresi gen
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)) # Menghitung nilai kuantil

# Logic penyesuaian rentang nilai berdasarkan Log2
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# TAHAP 2: Membuat Grup Data Sampel
group_info <- pData(gset)[["treatment:ch1"]] # Mengambil data berdasarkan kolom treatment
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)

# TAHAP 3: Membuat Desain Matrix Data
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Menentukan perbandingan biologis
grup_LPS <- nama_grup[1]
grup_saline <- nama_grup[2]
contrast_formula <- paste(grup_LPS, "-", grup_saline)
print(paste("Kontras yang dianalisis:", contrast_formula))

# TAHAP 4: Analisis Ekspresi Differensial (LIMMA)
fit <- lmFit(ex, design) # Membuat model linear
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design) # Membuat perbandingan grup
fit2 <- contrasts.fit(fit, contrast_matrix) # Membuat kontras pada model
fit2 <- eBayes(fit2) # Memakai Empirical Bayes

# Mengambil hasil akhir DEG
topTableResults <- topTable(fit2, adjust = "fdr", sort.by = "B", 
                            number = Inf, p.value = 0.01)
head(topTableResults)

# TAHAP 5: Pemberian Ulang Nama Gen
probe_ids <- rownames(topTableResults) # Mengambil ID Probe
gene_annotation <- AnnotationDbi::select(hgu133plus2.db, keys = probe_ids,
                                         columns = c("SYMBOL", "GENENAME"),
                                         keytype = "PROBEID") # Menggunakan database hgu133plus2.db
topTableResults$PROBEID <- rownames(topTableResults) # Combine dengan analisis LIMMA
topTableResults <- merge(topTableResults, gene_annotation, by = "PROBEID",
                         all.x = TRUE)
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

# TAHAP 6: Visualisasi Data DEG
# +) Pembuatan Boxplot
group_colors <- c("#fda002","#3076f3")
boxplot(ex, cex.axis = 0.7, col = group_colors, las = 2, outline = FALSE, 
        frame = FALSE, notch = TRUE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)")
legend("topright",legend = levels(gset$group),fill = unique(group_colors),
  cex = 0.8)

# +) Distribusi Nilai Ekspresi (Density Plot)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

# +) UMAP Plot
umap_input <- t(ex)
Total <- nrow(umap_input)

custom_umap_config <- umap.defaults
custom_umap_config$n_neighbors <- min(15, Total - 1) #Konfigurasi UMAP dengan n_neighbors (jika jumlah data lebih kecil dari jumlah input)
umap_result <- umap(umap_input, config = custom_umap_config)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

# +) Volcano Plot
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val <
                      0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val <
                      0.01] <- "DOWN"
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color =
                           status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Paru-paru Perokok")

# +) Heatmap
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)
mat_heatmap <- ex[top50$PROBEID, ]
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)
rownames(mat_heatmap) <- gene_label
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)
pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_col = 8,
  
  show_rownames = TRUE,
  fontsize_row = 6,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

# +) Menyimpan Hasil Pekerjaan
write.csv(topTableResults, "Hasil_GSE40885_DEG.csv")