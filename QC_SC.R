# Librerías necesarias
library(Seurat)


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

# mitochondrial genes and ribosomal proteins
spl_WT_obj[["percent.mt"]]  <- PercentageFeatureSet(spl_WT_obj, pattern = "^mt-")
spl_WT_obj[["percent.rb"]] <- PercentageFeatureSet(spl_WT_obj, pattern = "^Rps|^Rpl")
# doublets
doub_spl_WT <- read.csv("data/SCRUBLET/scrublet_WT_spl.csv",header = T,row.names = 1)
colnames(doub_spl_WT) <- c("Doublet_score","Is_doublet")
spl_WT_obj <- AddMetaData(spl_WT_obj,doub_spl_WT)

# graphs
VlnPlot(spl_WT_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(spl_WT_obj, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(spl_WT_obj, feature1 = "nFeature_RNA", feature2 = "Doublet_score")



# filtering
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "percent.mt")
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "percent.rb")
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Is_doublet")
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC == 'Pass','Low_nFeature',spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC != 'Pass' & spl_WT_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',spl_WT_obj@meta.data$QC,sep = ','),spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$percent.mt > 15 & spl_WT_obj@meta.data$QC == 'Pass','High_MT',spl_WT_obj@meta.data$QC)
spl_WT_obj[['QC']] <- ifelse(spl_WT_obj@meta.data$nFeature_RNA < 500 & spl_WT_obj@meta.data$QC != 'Pass' & spl_WT_obj@meta.data$QC != 'High_MT',paste('High_MT',spl_WT_obj@meta.data$QC,sep = ','),spl_WT_obj@meta.data$QC)
table(spl_WT_obj[['QC']])

# to plot after filtering:
# Por QC
FeatureScatter(spl_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QC")



## PrimPol KO spleen

spl_KO_obj[["percent.mt"]]  <- PercentageFeatureSet(spl_KO_obj, pattern = "^mt-")
spl_KO_obj[["percent.rb"]] <- PercentageFeatureSet(spl_KO_obj, pattern = "^Rps|^Rpl")
# doublets
doub_spl_KO <- read.csv("data/SCRUBLET/scrublet_KO_spl.csv",header = T,row.names = 1)
colnames(doub_spl_KO) <- c("Doublet_score","Is_doublet")
spl_KO_obj <- AddMetaData(spl_KO_obj,doub_spl_KO)

# graphs
VlnPlot(spl_KO_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
FeatureScatter(spl_KO_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(spl_KO_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(spl_KO_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(spl_KO_obj, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(spl_KO_obj, feature1 = "nFeature_RNA", feature2 = "Doublet_score")

# filtering
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC == 'Pass','Low_nFeature',spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC != 'Pass' & spl_KO_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',spl_KO_obj@meta.data$QC,sep = ','),spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$percent.mt > 15 & spl_KO_obj@meta.data$QC == 'Pass','High_MT',spl_KO_obj@meta.data$QC)
spl_KO_obj[['QC']] <- ifelse(spl_KO_obj@meta.data$nFeature_RNA < 500 & spl_KO_obj@meta.data$QC != 'Pass' & spl_KO_obj@meta.data$QC != 'High_MT',paste('High_MT',spl_KO_obj@meta.data$QC,sep = ','),spl_KO_obj@meta.data$QC)
table(spl_KO_obj[['QC']])

# to plot after filtering:
# Por QC
FeatureScatter(spl_KO_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QC")


# WT bone marrow
bm_WT_obj[["percent.mt"]]  <- PercentageFeatureSet(bm_WT_obj, pattern = "^mt-")
bm_WT_obj[["percent.rb"]] <- PercentageFeatureSet(bm_WT_obj, pattern = "^Rps|^Rpl")
# doublets
doub_bm_WT <- read.csv("data/SCRUBLET/scrublet_WT_BM.csv",header = T,row.names = 1)
colnames(doub_bm_WT) <- c("Doublet_score","Is_doublet")
bm_WT_obj <- AddMetaData(bm_WT_obj,doub_bm_WT)

# graphs
VlnPlot(bm_WT_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
FeatureScatter(bm_WT_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(bm_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(bm_WT_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(bm_WT_obj, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(bm_WT_obj, feature1 = "nFeature_RNA", feature2 = "Doublet_score")

# filtering
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC == 'Pass','Low_nFeature',bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC != 'Pass' & bm_WT_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',bm_WT_obj@meta.data$QC,sep = ','),bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$percent.mt > 15 & bm_WT_obj@meta.data$QC == 'Pass','High_MT',bm_WT_obj@meta.data$QC)
bm_WT_obj[['QC']] <- ifelse(bm_WT_obj@meta.data$nFeature_RNA < 500 & bm_WT_obj@meta.data$QC != 'Pass' & bm_WT_obj@meta.data$QC != 'High_MT',paste('High_MT',bm_WT_obj@meta.data$QC,sep = ','),bm_WT_obj@meta.data$QC)
table(bm_WT_obj[['QC']])

# to plot after filtering:
# Por QC
FeatureScatter(bm_WT_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QC")

# PrimPol KO bone marrow

bm_KO_obj[["percent.mt"]]  <- PercentageFeatureSet(bm_KO_obj, pattern = "^mt-")
bm_KO_obj[["percent.rb"]] <- PercentageFeatureSet(bm_KO_obj, pattern = "^Rps|^Rpl")
# doublets
doub_bm_KO <- read.csv("data/SCRUBLET/scrublet_KO_BM.csv",header = T,row.names = 1)
colnames(doub_bm_KO) <- c("Doublet_score","Is_doublet")
bm_KO_obj <- AddMetaData(bm_KO_obj,doub_bm_KO)

# graphs
VlnPlot(bm_KO_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
FeatureScatter(bm_KO_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(bm_KO_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(bm_KO_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(bm_KO_obj, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(bm_KO_obj, feature1 = "nFeature_RNA", feature2 = "Doublet_score")

# filtering
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$Is_doublet == 'True','Doublet','Pass')
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC == 'Pass','Low_nFeature',bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC != 'Pass' & bm_KO_obj@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',bm_KO_obj@meta.data$QC,sep = ','),bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$percent.mt > 15 & bm_KO_obj@meta.data$QC == 'Pass','High_MT',bm_KO_obj@meta.data$QC)
bm_KO_obj[['QC']] <- ifelse(bm_KO_obj@meta.data$nFeature_RNA < 500 & bm_KO_obj@meta.data$QC != 'Pass' & bm_KO_obj@meta.data$QC != 'High_MT',paste('High_MT',bm_KO_obj@meta.data$QC,sep = ','),bm_KO_obj@meta.data$QC)
table(bm_KO_obj[['QC']])

# to plot after filtering:
# Por QC
FeatureScatter(bm_KO_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QC")