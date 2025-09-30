########################################
############### PACKAGES ###############
########################################
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork)
library(reshape2)
library(hdf5r)

########################################
############### DATA IMPORT ############
########################################

# read the filtered Cell Ranger h5 matrices
spl_WT <- Read10X_h5("DATA/WT_SPsample_filtered_feature_bc_matrix.h5")
spl_KO <- Read10X_h5("DATA/KO_SP_sample_filtered_feature_bc_matrix.h5")
bm_WT  <- Read10X_h5("DATA/WT_BM_sample_filtered_feature_bc_matrix.h5")
bm_KO  <- Read10X_h5("DATA/KO_BM_sample_filtered_feature_bc_matrix.h5")

# create the appropriate Seurat
spl_WT_obj <- CreateSeuratObject(spl_WT, project = "spl_WT")
spl_KO_obj <- CreateSeuratObject(spl_KO, project = "spl_KO")
bm_WT_obj  <- CreateSeuratObject(bm_WT, project = "bm_WT")
bm_KO_obj  <- CreateSeuratObject(bm_KO, project = "bm_KO")

# let’s remove matrices to save memory
rm(spl_WT)
rm(spl_KO)
rm(bm_WT)
rm(bm_KO)


########################################
############### QUALITY ################
########################################

# note:
# this part is extended in QC_SR.R

#### SPLEEN WT

# mitochondrial genes and ribosomal proteins
spl_WT_obj[["percent.mt"]]  <- PercentageFeatureSet(spl_WT_obj, pattern = "^mt-")
spl_WT_obj[["percent.rb"]] <- PercentageFeatureSet(spl_WT_obj, pattern = "^Rps|^Rpl")
# doublets
doub_spl_WT <- read.csv("data/SCRUBLET/scrublet_WT_spl.csv",header = T,row.names = 1)
colnames(doub_spl_WT) <- c("Doublet_score","Is_doublet")
spl_WT_obj <- AddMetaData(spl_WT_obj,doub_spl_WT)

# filtering
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC == 'Pass','Low_nFeature',spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC != 'Pass' & spl_WT_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',spl_WT_obj@meta.data$QC,sep = ','),spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$percent.mt > 15 & spl_WT_obj@meta.data$QC == 'Pass','High_MT',spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC != 'Pass' & spl_WT_obj@meta.data$QC != 'High_MT',paste('High_MT',spl_WT_obj@meta.data$QC,sep = ','),spl_WT_obj@meta.data$QC)
table(spl_WT_obj[['QC']])




#### SPLEEN KO

## PrimPol KO spleen
spl_KO_obj[["percent.mt"]]  <- PercentageFeatureSet(spl_KO_obj, pattern = "^mt-")
spl_KO_obj[["percent.rb"]] <- PercentageFeatureSet(spl_KO_obj, pattern = "^Rps|^Rpl")
# doublets
doub_spl_KO <- read.csv("data/SCRUBLET/scrublet_KO_spl.csv",header = T,row.names = 1)
colnames(doub_spl_KO) <- c("Doublet_score","Is_doublet")
spl_KO_obj <- AddMetaData(spl_KO_obj,doub_spl_KO)


# filtering
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC == 'Pass','Low_nFeature',spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC != 'Pass' & spl_KO_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',spl_KO_obj@meta.data$QC,sep = ','),spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$percent.mt > 15 & spl_KO_obj@meta.data$QC == 'Pass','High_MT',spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC != 'Pass' & spl_KO_obj@meta.data$QC != 'High_MT',paste('High_MT',spl_KO_obj@meta.data$QC,sep = ','),spl_KO_obj@meta.data$QC)
table(spl_KO_obj[['QC']])


#### BONE MARROW WT

bm_WT_obj[["percent.mt"]]  <- PercentageFeatureSet(bm_WT_obj, pattern = "^mt-")
bm_WT_obj[["percent.rb"]] <- PercentageFeatureSet(bm_WT_obj, pattern = "^Rps|^Rpl")
# doublets
doub_bm_WT <- read.csv("data/SCRUBLET/scrublet_WT_BM.csv",header = T,row.names = 1)
colnames(doub_bm_WT) <- c("Doublet_score","Is_doublet")
bm_WT_obj <- AddMetaData(bm_WT_obj,doub_bm_WT)

# filtering
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC == 'Pass','Low_nFeature',bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC != 'Pass' & bm_WT_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',bm_WT_obj@meta.data$QC,sep = ','),bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$percent.mt > 15 & bm_WT_obj@meta.data$QC == 'Pass','High_MT',bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC != 'Pass' & bm_WT_obj@meta.data$QC != 'High_MT',paste('High_MT',bm_WT_obj@meta.data$QC,sep = ','),bm_WT_obj@meta.data$QC)
table(bm_WT_obj[['QC']])

#### BONE MARROW PRIMPOL KO

bm_KO_obj[["percent.mt"]]  <- PercentageFeatureSet(bm_KO_obj, pattern = "^mt-")
bm_KO_obj[["percent.rb"]] <- PercentageFeatureSet(bm_KO_obj, pattern = "^Rps|^Rpl")
# doublets
doub_bm_KO <- read.csv("data/SCRUBLET/scrublet_KO_BM.csv",header = T,row.names = 1)
colnames(doub_bm_KO) <- c("Doublet_score","Is_doublet")
bm_KO_obj <- AddMetaData(bm_KO_obj,doub_bm_KO)

# filtering
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC == 'Pass','Low_nFeature',bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC != 'Pass' & bm_KO_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',bm_KO_obj@meta.data$QC,sep = ','),bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$percent.mt > 15 & bm_KO_obj@meta.data$QC == 'Pass','High_MT',bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC != 'Pass' & bm_KO_obj@meta.data$QC != 'High_MT',paste('High_MT',bm_KO_obj@meta.data$QC,sep = ','),bm_KO_obj@meta.data$QC)
table(bm_KO_obj[['QC']])


########################################
############### SPLEEN #################
########################################

# filter
spl_WT_obj_pass <- subset(spl_WT_obj, subset = QC == "Pass")
spl_KO_obj_pass <- subset(spl_KO_obj, subset = QC == "Pass")

# integration
# creat a list
spl_list <- list()
spl_list[["spl_WT"]] <- spl_WT_obj_pass
spl_list[["spl_KO"]] <- spl_KO_obj_pass

# normalize/find HVG
for (i in 1:length(spl_list)) {
  spl_list[[i]] <- NormalizeData(spl_list[[i]], verbose = F)
  spl_list[[i]] <- FindVariableFeatures(spl_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}
# to find integration anchors and actually perform integration
# it takes time...
spl_anchors <- FindIntegrationAnchors(object.list = spl_list, dims = 1:30)
spl_seurat <- IntegrateData(anchorset = spl_anchors, dims = 1:30)
# to save ram
rm(spl_list)
rm(spl_anchors)

# before integration
DefaultAssay(spl_seurat) <- "RNA"
spl_seurat <- NormalizeData(spl_seurat, verbose = F)
spl_seurat <- FindVariableFeatures(spl_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
spl_seurat <- ScaleData(spl_seurat, verbose = F)
spl_seurat <- RunPCA(spl_seurat, npcs = 30, verbose = F)
spl_seurat <- RunUMAP(spl_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(spl_seurat,reduction = "umap") +
  plot_annotation(title = "Spleen, before integration")

# the integrated assay (it’s already normalized and HVGs are selected)
DefaultAssay(spl_seurat) <- "integrated"
spl_seurat <- ScaleData(spl_seurat, verbose = F)
spl_seurat <- RunPCA(spl_seurat, npcs = 30, verbose = F)
spl_seurat <- RunUMAP(spl_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(spl_seurat, reduction = "umap") + 
  plot_annotation(title = "Spleen, after integration (Seurat 3)")

# split
DimPlot(spl_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

# let's cluster
spl_seurat <- FindNeighbors(spl_seurat, dims = 1:30, k.param = 10, verbose = F)
spl_seurat <- FindClusters(spl_seurat, verbose = F)
# spl_seurat <- FindClusters(spl_seurat, verbose = F, resolution = 0.1)
DimPlot(spl_seurat,label = T) + NoLegend()



# to test several clusters
for (res in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)) {
  spl_seurat <- FindClusters(spl_seurat, resolution = res, verbose = FALSE)
  print(DimPlot(spl_seurat, label = TRUE) + ggtitle(paste("Resolution", res)))
}

# table with the cluster
count_table <- table(spl_seurat@meta.data$seurat_clusters, spl_seurat@meta.data$orig.ident)
count_table

### annotation

# Extraer las matrices de cada muestra
expr_WT <- GetAssayData(spl_seurat, assay = "RNA", layer = "data.spl_WT")
expr_KO <- GetAssayData(spl_seurat, assay = "RNA", layer = "data.spl_KO")

# Combina columnas (células) de ambos
expr <- cbind(expr_WT, expr_KO)


spl_sce <- SingleCellExperiment(
  assays = list(logcounts = expr),
  colData = spl_seurat@meta.data
)


# mouse inmune reference
ref <- celldex::ImmGenData()
# otro puede ser:
# MouseRNAseqData()
# ImmGenData

# Corre SingleR para anotar
pred_spl <- SingleR(test = spl_sce, ref = ref, labels = ref$label.fine)

# Añade la anotación a tu objeto Seurat
spl_seurat$SingleR_labels <- pred_spl$labels

# Acorta los labels (esto es solo un ejemplo, puedes adaptarlo)
spl_seurat$SingleR_clean <- gsub("\\s*\\([^\\)]+\\)", "", spl_seurat$SingleR_labels)

# Visualiza con nombres más cortos
DimPlot(spl_seurat, group.by = "SingleR_clean", label = TRUE, repel = TRUE) + NoLegend()


# para anotar los clusteres

# 1. Extraer metadata
meta <- spl_seurat@meta.data

# 2. Crear tabla con la anotación más frecuente por clúster
annot_df <- meta %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster) %>%
  count(SingleR_clean) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. Crear un vector de anotación por clúster
annot_vec <- setNames(annot_df$SingleR_clean, annot_df$cluster)

# 4. Asignar anotaciones a cada célula usando su clúster
meta$cluster_annotation <- annot_vec[as.character(meta$seurat_clusters)]

# 5. Asegurarse de que se guarda en el objeto Seurat
spl_seurat@meta.data <- meta

DimPlot(spl_seurat, group.by = "cluster_annotation", label = TRUE, repel = TRUE) + NoLegend()




# if i want to check specific genes in cluster
# volcano plot 
VlnPlot(
  spl_seurat,
  features = "H2ax",
  group.by = "cluster_annotation",   # usa la anotación en lugar del número
  split.by = "orig.ident"            # WT vs KO
) + 
  ggtitle("H2ax expression by annotated clusters (WT vs KO)")


FeaturePlot(
  spl_seurat,
  features = "H2ax",
  reduction = "umap",
  split.by = "orig.ident",   # genera un panel para WT y otro para KO
  cols = c("lightgrey", "red")
)



# para buscar marcadores de cada cluster:

DefaultAssay(spl_seurat) <- "RNA"
spl_seurat <- JoinLayers(spl_seurat)
# findallmarkers
markers_spl <- FindAllMarkers(spl_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers_spl)
# write.csv(markers_spl, "markers_spleen.csv")

# top genes per cluster
top10 <- markers_spl %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
# Opcional: heatmap con top genes
DoHeatmap(spl_seurat, features = top10$gene) + NoLegend()

# para plotear todo

# Paso 1: Define las anotaciones para los 16 clústeres
annotations <- c(
  "B cells",                    # cluster 0
  "T cells",                    # cluster 1
  "T cell activation",          # cluster 2
  "Neutrophils",                # cluster 3
  "T cells",                    # cluster 4
  "Vesicle-mediated transport", # cluster 5
  "Antigen presentation",       # cluster 6
  "T cell activation",          # cluster 7
  "ER-stress",                  # cluster 8
  "Chromosome seggregation",    # cluster 9
  "B cell activation",          # cluster 10
  "Neutrophil activation",      # cluster 11
  "NK cells",                   # cluster 12
  "T cells/ribosomes",          # cluster 13
  "Myeloid cell differentiation", # cluster 14
  "Leukocyte/granulocyte"       # cluster 15
)

# Asegurar que las anotaciones estén bien nombradas
names(annotations) <- as.character(0:15)

# Paso 2: Obtener el vector de clúster por celda
cluster_ids <- spl_seurat$seurat_clusters
names(cluster_ids) <- colnames(spl_seurat)

# Paso 3: Mapear anotaciones por celda
cluster_annotations <- annotations[as.character(cluster_ids)]

# Paso 4: Añadir al metadata
spl_seurat@meta.data$cluster_annotation_string <- cluster_annotations

# Visualizar
DimPlot(spl_seurat, group.by = "cluster_annotation_string", label = TRUE, repel = TRUE) + NoLegend()





# para chequear marcadores:
FeaturePlot(spl_seurat, features = c("Adgre1", "Cd68", "Hba-a1", "Ly6g"), label = TRUE)
VlnPlot(spl_seurat, features = c("Adgre1", "Ly6g", "Klf1"), group.by = "cluster_annotation")

# usando find markers
markers <- FindMarkers(bm_seurat, ident.1 = "Macrophages", group.by = "cluster_annotation")
head(markers)




########################################
############ BONE MARROW ###############
########################################

# integration
# creat a list
bm_list <- list()
bm_list[["bm_WT"]] <- bm_WT_obj
bm_list[["bm_KO"]] <- bm_KO_obj

# normalize/find HVG
for (i in 1:length(bm_list)) {
  bm_list[[i]] <- NormalizeData(bm_list[[i]], verbose = F)
  bm_list[[i]] <- FindVariableFeatures(bm_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}
# to find integration anchors and actually perform integration
# it takes time...
bm_anchors <- FindIntegrationAnchors(object.list = bm_list, dims = 1:30)
bm_seurat <- IntegrateData(anchorset = bm_anchors, dims = 1:30)
# to save ram
rm(bm_list)
rm(bm_anchors)

# before integration
DefaultAssay(bm_seurat) <- "RNA"
bm_seurat <- NormalizeData(bm_seurat, verbose = F)
bm_seurat <- FindVariableFeatures(bm_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
bm_seurat <- ScaleData(bm_seurat, verbose = F)
bm_seurat <- RunPCA(bm_seurat, npcs = 30, verbose = F)
bm_seurat <- RunUMAP(bm_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(bm_seurat,reduction = "umap") +
  plot_annotation(title = "Bone marrow, before integration")

# the integrated assay (it’s already normalized and HVGs are selected)
DefaultAssay(bm_seurat) <- "integrated"
bm_seurat <- ScaleData(bm_seurat, verbose = F)
bm_seurat <- RunPCA(bm_seurat, npcs = 10, verbose = F)
bm_seurat <- RunUMAP(bm_seurat, reduction = "pca", dims = 1:10, verbose = F)

DimPlot(bm_seurat, reduction = "umap") + 
  plot_annotation(title = "Bone marrow, after integration (Seurat 3)")

# split
DimPlot(bm_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

# let's cluster
bm_seurat <- FindNeighbors(bm_seurat, dims = 1:10, k.param = 10, verbose = F)
# less resolution, less clusters
bm_seurat <- FindClusters(bm_seurat, verbose = F, resolution = 0.1)
DimPlot(bm_seurat,label = T) + NoLegend()

plot_integrated_clusters(bm_seurat)

# number of cells for each cluster
count_table_bm <- table(bm_seurat@meta.data$seurat_clusters, bm_seurat@meta.data$orig.ident)
count_table_bm


### annotation

# Extraer las matrices de cada muestra
expr_WT_bm <- GetAssayData(bm_seurat, assay = "RNA", layer = "data.bm_WT")
expr_KO_bm <- GetAssayData(bm_seurat, assay = "RNA", layer = "data.bm_KO")

# Combina columnas (células) de ambos
expr_bm <- cbind(expr_WT_bm, expr_KO_bm)

bm_sce <- SingleCellExperiment(
  assays = list(logcounts = expr_bm),
  colData = bm_seurat@meta.data
)

# annotation

# mouse inmune reference
ref <- celldex::ImmGenData()

# Corre SingleR para anotar
pred_bm <- SingleR(test = bm_sce, ref = ref, labels = ref$label.fine)

# Añade la anotación a tu objeto Seurat
bm_seurat$SingleR_labels <- pred_bm$labels
# limpiamos
bm_seurat$SingleR_clean <- gsub("\\s*\\([^\\)]+\\)", "", bm_seurat$SingleR_labels)

# Visualiza con nombres más cortos
DimPlot(bm_seurat, group.by = "SingleR_clean", label = TRUE, repel = TRUE) + NoLegend()


# para anotar los clusteres

# 1. Extraer metadata
meta <- bm_seurat@meta.data

# 2. Crear tabla con la anotación más frecuente por clúster
annot_df <- meta %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster) %>%
  count(SingleR_clean) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. Crear un vector de anotación por clúster
annot_vec <- setNames(annot_df$SingleR_clean, annot_df$cluster)

# 4. Asignar anotaciones a cada célula usando su clúster
meta$cluster_annotation <- annot_vec[as.character(meta$seurat_clusters)]

# 5. Asegurarse de que se guarda en el objeto Seurat
bm_seurat@meta.data <- meta

DimPlot(bm_seurat, group.by = "cluster_annotation", label = TRUE, repel = TRUE) + NoLegend()

# ploteamos los porcentajes
plot_integrated_clusters_annotated(bm_seurat)

# volcano plot 
VlnPlot(
  bm_seurat,
  features = "H2ax",
  group.by = "cluster_annotation",   # usa la anotación en lugar del número
  split.by = "orig.ident"            # WT vs KO
) + 
  ggtitle("H2ax expression by annotated clusters (WT vs KO)")






# to annotate manually

# para buscar marcadores de cada cluster:

DefaultAssay(bm_seurat) <- "RNA"
bm_seurat <- JoinLayers(bm_seurat)
# findallmarkers
markers_bm <- FindAllMarkers(bm_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers_bm)
#write.csv(markers_bm, "OUTPUT/markers_bonemarrow.csv")

# para plotear todo


# Paso 1: Define las anotaciones para los 13 clústeres
annotations <- c(
  "Ribosomal-B cells",          # cluster 0
  "T cells",                    # cluster 1
  "T cell activation",          # cluster 2
  "Neutrophils/Macrophages",    # cluster 3
  "T cells",                    # cluster 4
  "Vesicle-mediated transport", # cluster 5
  "Antigen presentation/Macrophages",# cluster 6
  "T cell activation",          # cluster 7
  "ER-stress/Pancreas",         # cluster 8
  "Chromosome seggregation",    # cluster 9
  "B cell activation",          # cluster 10
  "Neutrophil activation",      # cluster 11
  "NK cells",                   # cluster 12
  "T cell",                     # cluster 13
  "Myeloid cell differentiation",      # cluster 14
  "Granulocyte",                   # cluster 15
)

# Asegurar que las anotaciones estén bien nombradas
names(annotations) <- as.character(0:15)

# Paso 2: Obtener el vector de clúster por celda
cluster_ids <- bm_seurat$seurat_clusters
names(cluster_ids) <- colnames(bm_seurat)

# Paso 3: Mapear anotaciones por celda
cluster_annotations <- annotations[as.character(cluster_ids)]

# Paso 4: Añadir al metadata
bm_seurat@meta.data$cluster_annotation_string <- cluster_annotations

# Visualizar
DimPlot(bm_seurat, group.by = "cluster_annotation_string", label = TRUE, repel = TRUE) + NoLegend()

