---
title: "R Notebook"
output:
  html_document:
    df_print: paged
---

```{r}
# Set the knitr options
knitr::opts_chunk$set(
  echo = TRUE,
  error = TRUE,
  message = FALSE,
  warning = FALSE,
  results = "hold",
  fig.keep = "all",
  fig.show = "hold",
  comment = NA,
  dev = "svg"
)

library(limma)
library(Seurat)
library(dplyr)
library(RColorBrewer)
# library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(gghighlight)
# library(tibble)
library(ggh4x)
library(ComplexHeatmap)

# library(Matrix)
# library(NMF)
# library(EnhancedVolcano)
```

# Set path and load file
```{r}
# input_path <- r"(\\store.erasmusmc.nl\department\bmw\chien_data\Lab\Data_and_Analysis\Pin_Rui_Su\2 Data\PRS_220201_scRNAseq_Tripolar cells Plate 1-4\230223 process_from Chan mapped Seurat)"
snvObj <- readRDS(r"(\\store.erasmusmc.nl\department\bmw\chien_data\Lab\Data_and_Analysis\Pin_Rui_Su\2 Data\PRS_220201_scRNAseq_Tripolar cells Plate 1-4\230329 process_from Chan mapped SNV\combined_snv.rds)")
sObj <- readRDS(r"(\\store.erasmusmc.nl\department\bmw\chien_data\Lab\Data_and_Analysis\Pin_Rui_Su\2 Data\PRS_220201_scRNAseq_Tripolar cells Plate 1-4\220728 process_with SNV inferred factors\s143_s145_combined_seurat.rds)")
# write.csv(as.matrix(GetAssayData(sObj, assay = "RNA")), "RNA_143,145.csv")
```

# Output cell information
```{r}
snvObj$condition <- sObj$condition
snvObj$nmf_factors <- sObj$nmf_factors
snvObj$lineage <- sObj$lineage

sumTable <- sapply(unique(snvObj@meta.data$orig.ident), FUN = function(x) {
  dat <- snvObj@meta.data$condition[snvObj@meta.data$orig.ident==x]
  list(Bipolar = sum(dat=="Bipolar"), Tripolar = sum(dat=="Tripolar"), Control = sum(dat=="Control"))})
print("After QC")
print(t(sumTable))
```

# Run PCA
```{r}
snvObj <- NormalizeData(snvObj) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  JackStraw(dims = 50, num.replicate = 100) %>% 
  ScoreJackStraw(dims = 1:50)
JackStrawPlot(snvObj, dims = 1:50)
ElbowPlot(snvObj)
snvObj <- FindNeighbors(snvObj, dims = 1:10)
snvObj <- FindClusters(snvObj, resolution = 0.3)

DimPlot(snvObj, reduction = "pca", group.by = "orig.ident")
DimPlot(snvObj, reduction = "pca", group.by = "condition", label = TRUE)
DimPlot(snvObj, reduction = "pca", group.by = "nmf_factors", label = TRUE)

snvObj@meta.data$type_scluster <- str_c(snvObj@meta.data$condition, snvObj@meta.data$seurat_clusters)
snvObj@meta.data$type_scluster_NMF <- str_c(snvObj@meta.data$condition, "S", snvObj@meta.data$seurat_clusters, "N", snvObj@meta.data$nmf_factors)
```


# Run non-liner dimensional reduction
```{r}
snvObj <- RunUMAP(snvObj, dims = 1:50)
```

# Plot UMAP
```{r, fig.height = 5}
#| crop: true
# png(paste0(input_path, "\\s143s145_umap_batch.png"))
DimPlot(snvObj, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 2) + ggtitle("Batch") + theme(aspect.ratio=1) +
  scale_color_manual(values = alpha(c(brewer.pal(n = 5, name = "Dark2")), .5))
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_type.png"))
DimPlot(snvObj, reduction = "umap", group.by = "condition", label = TRUE, pt.size = 2) + ggtitle("Mitosis type") + theme(aspect.ratio=1) +
  scale_color_manual(values = alpha(c(brewer.pal(n = 5, name = "Dark2")), .5))
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_typefactor.png"))
DimPlot(snvObj, reduction = "umap", group.by = "lineage", label = TRUE, pt.size = 2) + ggtitle("Mitosis type + Factor") + theme(aspect.ratio=1) +
  scale_color_manual(values = alpha(c(brewer.pal(n = 6, name = "Dark2")), .5))
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_factor.png"))
DimPlot(snvObj, reduction = "umap", group.by = "nmf_factors", label = TRUE, pt.size = 2) + ggtitle("NMF Factor") + theme(aspect.ratio=1) +
  scale_color_manual(values = alpha(c(brewer.pal(n = 5, name = "Dark2")), .5))
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_scluster.png"))
DimPlot(snvObj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 2) + ggtitle("Seurat cluster") + theme(aspect.ratio=1) +
  scale_color_manual(values = alpha(c(brewer.pal(n = 5, name = "Dark2")), .5))
# dev.off()
```

# Plot heatmap
```{r}
mat <- as.matrix(GetAssayData(snvObj))

## scale the rows
mat <- t(scale(t(mat)))

cluster_anno <- snvObj@meta.data$nmf_factors %>% as.factor()

lt <- setNames(as.numeric(unique(snvObj@meta.data$nmf_factors)), brewer.pal(n = 3, name = "Dark2"))
heatmap(mat, ColSideColors = names(lt[cluster_anno]))

```
```{r}
Heatmap(mat, name = "Expression",  
        column_split = cluster_anno,
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(1, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        row_names_gp = gpar(fontsize = 2),
        top_annotation = HeatmapAnnotation(NMF = cluster_anno,
                                            annotation_name_gp = gpar(fontsize = 10)),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 0.01)
```



#########################################################################################

# Plot UMAP on specific clusters
```{r, fig.height=5}
itemList <- unique(sObj@meta.data$type_scluster_NMF)

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("BipolarS1N1", "BipolarS0N1")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("TripolarS1N1", "BipolarS1N1")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("BipolarS0N1", "BipolarS0N3")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("TripolarS1N2", "TripolarS1N1")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("TripolarS2N2", "TripolarS1N2")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("TripolarS1N2", "TripolarS1N1", "BipolarS1N2", "BipolarS1N1")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("BipolarS1N2", "BipolarS1N1")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("BipolarS1N2", "BipolarS1N1", "BipolarS0N2", "BipolarS0N1", "BipolarS0N3")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
hlList <- c("TripolarS2N2", "BipolarS2N2")
colorList <- setNames(rep("grey", length(itemList)), itemList)
colorList[hlList] <- brewer.pal(n = 8, name = "Dark2")[1:length(hlList)]
DimPlot(sObj, reduction = "umap", group.by = "type_scluster_NMF", label = TRUE, pt.size = 2) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust=0.5)) +
  ggtitle("Cell type", subtitle = paste(hlList, collapse = ', ')) + 
  scale_color_manual(values = alpha(colorList, .5)) +
  gghighlight(type_scluster_NMF %in% hlList, use_direct_label = FALSE)
# dev.off()

```

# Plot distribution
```{r}
# png(paste0(input_path, "\\dis_batch_factor.png"))
freq_table <- table(sObj@meta.data$nmf_factors, sObj@meta.data$orig.ident)
barplot(freq_table, beside = TRUE,
        main = "NMF factor",
        xlab = "Batch", ylab = "Counts",
        ylim = c(0,150), col = brewer.pal(n = nrow(freq_table), name = "Dark2"),
        legend = TRUE, args.legend = list(x = "top", ncol = 3, y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_type_factor.png"))
freq_table <- table(sObj@meta.data$nmf_factors, sObj@meta.data$condition)
barplot(freq_table, beside = TRUE,
        main = "NMF factor",
        xlab = "Mitosis type", ylab = "Counts",
        ylim = c(0,220), col = brewer.pal(n = nrow(freq_table), name = "Dark2"),
        legend = TRUE, args.legend = list(x = "top", ncol = 3, y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_factor_type.png"))
freq_table <- table(sObj@meta.data$condition, sObj@meta.data$nmf_factors)
barplot(freq_table, beside = TRUE,
        main = "Mitosis type",
        xlab = "NMF factor", ylab = "Counts",
        ylim = c(0,200), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_rel_factor_type.png"))
freq_table <- table(sObj@meta.data$condition, sObj@meta.data$nmf_factors)
freq_table2 <- freq_table / rowSums(freq_table) * 100
barplot(freq_table2, beside = TRUE,
        main = "Mitosis type",
        xlab = "NMF factor", ylab = "Percentage (%)",
        ylim = c(0,100), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_scluster_type.png"))
freq_table <- table(sObj@meta.data$condition, sObj@meta.data$seurat_clusters)
barplot(freq_table, beside = TRUE,
        main = "Mitosis type",
        xlab = "Seurat cluster", ylab = "Counts",
        ylim = c(0,250), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_rel_scluster_type.png"))
freq_table <- table(sObj@meta.data$condition, sObj@meta.data$seurat_clusters)
freq_table2 <- freq_table / rowSums(freq_table) * 100
barplot(freq_table2, beside = TRUE,
        main = "Mitosis type",
        xlab = "Seurat cluster", ylab = "Percentage (%)",
        ylim = c(0,100), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_factor_scluster.png"))
freq_table <- table(sObj@meta.data$seurat_clusters, sObj@meta.data$nmf_factors)
barplot(freq_table, beside = TRUE,
        main = "Seurat cluster",
        xlab = "NMF factor", ylab = "Counts",
        ylim = c(0,200), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()

# png(paste0(input_path, "\\dis_rel_factor_scluster.png"))
freq_table <- table(sObj@meta.data$seurat_clusters, sObj@meta.data$nmf_factors)
freq_table2 <- freq_table / rowSums(freq_table) * 100
barplot(freq_table2, beside = TRUE,
        main = "Seurat cluster",
        xlab = "NMF factor", ylab = "Percentage (%)",
        ylim = c(0,100), col = brewer.pal(n = nrow(freq_table), name = "Dark2")[1:nrow(freq_table)],
        legend = TRUE, args.legend = list(x = "top", ncol = nrow(freq_table), y.intersp = 0))
# dev.off()
```

# Run limma
# Function
```{r}
plotVolcano <- function(dgtable, x = "logFC", y = "P.Value", title = "Volcano plot", pCutoffCol = "adj.P.Val", CUTOFFfc = CUTOFFfc) {
  EnhancedVolcano::EnhancedVolcano(dgtable, x = x, y = y, pCutoffCol = pCutoffCol, lab = rownames(degtable),
                                   pCutoff = 0.05, FCcutoff = CUTOFFfc,
                                   xlim = c(-max(abs(dgtable[[x]])), max(abs(dgtable[[x]]))),
                                   ylim = c(-0.01, max(-log10(dgtable[[y]]))),
                                   title = title, subtitle = "", caption = "", legendPosition = "none",
                                   # selectLab = sigGenes[which(sigGenes %in% annotationVolcanos)],
                                   col = c("lightgray", "lightgray", "lightsteelblue2", "brown2"), colAlpha = 0.5,
                                   gridlines.major = F, gridlines.minor = F, drawConnectors = T, arrowheads = F, boxedLabels = F, max.overlaps = Inf,
                                   pointSize = 2, labSize = 3, axisLabSize = 12, titleLabSize = 12  )}
```

# Pairwise DEGs for Tripolar vs Bipolar
```{r dev='svg', fig.width=8, fig.height=8}
# The normalized expression
expressionData <- as.matrix(GetAssayData(sObj, assay = "RNA"))
print(paste("There are", sum(rowSums(expressionData) == 0), "genes have all zeros", sep = " "))
filteredData <- expressionData[which(rowSums(expressionData) != 0), ]

# Create model matrix
design <- model.matrix(~0 + condition + orig.ident, data = sObj@meta.data)
colnames(design) <- gsub("condition|orig.ident", "", colnames(design))
print(head(design))

# Fit the model
model <- lmFit(filteredData, design)
print(head(coef(model))) # means in each sample for each gene

# Identify DEGS
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(TrivsBi = Tripolar - Bipolar, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable_pooled_MPvsBP <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable_pooled_MPvsBP, "DEG_pooled_MPvsBP.csv")

# Output results
plotVolcano(degtable, title = "Multipolar vs Bipolar", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)
```

# DEGs correlated with Seurat cluster identifiers or NMF factors
```{r dev='svg', fig.width=8, fig.height=8}
# Correlation among Seurat cluster
design <- model.matrix(~as.numeric(sObj@meta.data$seurat_clusters) + sObj@meta.data$orig.ident) # Create model matrix
colnames(design) <- gsub("as.numeric(sObj@meta.data$seurat_clusters)", "scluster", colnames(design))
colnames(design) <- gsub("sObj@meta.data$orig.idents145", "145", colnames(design))
print(head(design))

model <- lmFit(filteredData, design) # Fit the model
print(head(coef(model))) # means in each sample for each gene

fittedmodel <- contrasts.fit(model, coef = 2) # test "seurat cluster" coefficient
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf)
# write.csv(degtable, "DEG_corr_seuratcluster.csv")
plotVolcano(degtable, title = "Correlation with Seurat cluster", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

# Correlation among NMF factor
design <- model.matrix(~as.numeric(sObj@meta.data$nmf_factors) + sObj@meta.data$orig.ident) # Create model matrix
colnames(design) <- gsub("as.numeric(sObj@meta.data$nmf_factors)", "nfactor", colnames(design))
colnames(design) <- gsub("sObj@meta.data$orig.idents145", "145", colnames(design))
print(head(design))

model <- lmFit(filteredData, design) # Fit the model
print(head(coef(model))) # means in each sample for each gene

fittedmodel <- contrasts.fit(model, coef = 2) # test "seurat cluster" coefficient
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf)
# write.csv(degtable, "DEG_MP2vsBP0.csv")
plotVolcano(degtable, title = "Correlation with NMF factor", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)
```

# Plot UMAP on highly-correlated genes
```{r, fig.height = 3, fig.width = 15}
# png(paste0(input_path, "\\s143s145_umap_batch.png"))
gene2plot <- c("PMEPA1", "HES1", "PSCA", "KRT17", "SERPINB2")
FeaturePlot(sObj, reduction = "umap", features = gene2plot, cols = alpha(c("lightgrey", "brown2"), .7), ncol = 5, pt.size = 1)
# dev.off()

# png(paste0(input_path, "\\s143s145_umap_batch.png"))
gene2plot <- c("ANLN", "AREG", "ACTB", "KRT81", "MTATP6P1")
FeaturePlot(sObj, reduction = "umap", features = gene2plot, cols = alpha(c("lightgrey", "brown2"), .7), ncol = 5, pt.size = 1)
# dev.off()
```


# Pairwise DEGs with Seurat cluster identifiers
```{r dev='svg', fig.width=8, fig.height=8}
design <- model.matrix(~0 + type_scluster + orig.ident, data = sObj@meta.data) # Create model matrix
colnames(design) <- gsub("type_scluster|orig.ident", "", colnames(design))
print(head(design))

model <- lmFit(filteredData, design) # Fit the model
print(head(coef(model))) # means in each sample for each gene

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar2 - Bipolar0, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsBP0.csv")
plotVolcano(degtable, title = "Multipolar2 vs Bipolar0", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar2 - Tripolar1, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsMP1.csv")
plotVolcano(degtable, title = "Multipolar2 vs Multipolar1", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar1 - Bipolar1, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP1vsBP1.csv")
plotVolcano(degtable, title = "Multipolar1 vs Bipolar1", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Bipolar1 - Bipolar0, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_BP1vsBP0.csv")
plotVolcano(degtable, title = "Bipolar1 vs Bipolar0", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)
```

# Individual DEGs for all samples with Seurat cluster identifiers
```{r dev='svg', fig.width=8, fig.height=8}
filteredData_rep <- filteredData # Create a matrix that replicates data for average comparison
colnames(filteredData_rep) <- str_c(colnames(filteredData_rep), "_r")
filteredData_av <- cbind(filteredData, filteredData_rep) 
metadata_rep <- sObj@meta.data
rownames(metadata_rep) <- str_c(rownames(metadata_rep), "_r")
metadata_rep$type_scluster <- "All"
metadata_av <- rbind(sObj@meta.data, metadata_rep)

design <- model.matrix(~0 + type_scluster + orig.ident, data = metadata_av) # Create model matrix
colnames(design) <- gsub("type_scluster|orig.ident", "", colnames(design))
print(head(design))

model <- lmFit(filteredData_av, design) # Fit the model
print(head(coef(model))) # means in each sample for each gene

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar2 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_MP2.csv")
plotVolcano(degtable, title = "Multipolar2 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar1 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_MP1.csv")
plotVolcano(degtable, title = "Multipolar1 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Tripolar0 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_MP0.csv")
plotVolcano(degtable, title = "Multipolar0 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Bipolar2 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_BP2.csv")
plotVolcano(degtable, title = "Bipolar2 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Bipolar1 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_BP1.csv")
plotVolcano(degtable, title = "Bipolar1 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(Bipolar0 - All, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_av_BP0.csv")
plotVolcano(degtable, title = "Bipolar0 vs Average", CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)
```

# DEGs with both Seurat clusters and NMF cluster identifiers
```{r dev='svg', fig.width=8, fig.height=8}
design <- model.matrix(~0 + type_scluster_NMF + orig.ident, data = sObj@meta.data) # Create model matrix
colnames(design) <- gsub("type_scluster_NMF|orig.ident", "", colnames(design))
print(head(design))

model <- lmFit(filteredData, design) # Fit the model
print(head(coef(model))) # means in each sample for each gene

contrastFml <- "BipolarS1N1 - BipolarS0N1"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_BPS1N1vsBPS0N1.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "TripolarS1N1 - BipolarS1N1"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsMP1.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "BipolarS0N1 - BipolarS0N3"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MPS0N1vsMPS0N3.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "TripolarS1N2 - TripolarS1N1"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsMP1.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "BipolarS1N2 - BipolarS1N1"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsMP1.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "TripolarS2N2 - TripolarS1N2"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
degtable_MPS2N2vsMPS1N2 <- degtable
# write.csv(degtable, "DEG_MPS2N2vsMPS1N2.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)

contrastFml <- "TripolarS2N2 - BipolarS2N2"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
degtable <- topTable(fittedmodel, number = Inf, adjust.method = "BH")
# write.csv(degtable, "DEG_MP2vsMP1.csv")
plotVolcano(degtable, title = gsub("-", "vs", contrastFml), CUTOFFfc = log2(1.5)) + theme(aspect.ratio=1)
```

# Test DEG correlation by correlation plot (not yet done)
```{r, fig.height = 8}
a <- degtable_pooled_MPvsBP %>% filter(adj.P.Val <= 0.05)
a$rowname <- row.names(a)
b <- degtable_MPS2N2vsMPS1N2 %>% filter(adj.P.Val <= 0.05)
b$rowname <- row.names(b)
c <- merge(a, b, by = "rowname")

ggplot(c, aes(x = logFC.x, y = logFC.y)) +
  geom_point(size = 3, color = "lightsteelblue2", alpha = .5) + 
  geom_text(aes(label = rowname), size = 3, check_overlap = TRUE) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1),
        aspect.ratio = 1) +
  scale_x_continuous(limits = c(-(max(abs(a$logFC))), (max(abs(a$logFC)))), breaks = scales::breaks_width(0.2)) +
  scale_y_continuous(limits = c(-(max(abs(b$logFC))), (max(abs(b$logFC)))), breaks = scales::breaks_width(0.2)) +
  labs(x = "DEG MP vs BP",
       y = "DEG MP-S2N2 vs MP-S1N2") +
  coord_axes_inside(labels_inside = TRUE)


```



# Test the DEG correlation by Venn diagram (not yet done)
```{r}
# Plot Venn diagram
## With fold change threshold
contrastFml <- "BipolarS1N1 - BipolarS0N1"
fittedmodel <- contrasts.fit(model, contrasts = makeContrasts(contrastFml, levels = colnames(design)))
fittedmodel <- eBayes(fittedmodel)
dt <- decideTests(fittedmodel, adjust.method = "BH", lfc = log2(1.5))
summary(dt)
vennDiagram(dt, include = c("up","down"), circle.col = c("darkorange1", "seagreen3", "dodgerblue1"))

## Without fold change threshold
dt <- decideTests(fittedmodel, adjust.method = "BH")
summary(dt)
vennDiagram(dt, include = c("up","down"), circle.col = c("darkorange1", "seagreen3", "dodgerblue1"))
```


```{r}
sessionInfo()
```

