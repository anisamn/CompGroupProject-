#set up environment and load packages
#library(Seurat)
#library(patchwork)


library(dplyr)
library(Seurat)
library(patchwork)


#run this line once before running code below 
#reticulate::py_install(packages ='umap-learn')

### SETUP THE SEURAT OBJECT ###

#get full directory to mouse heart GEO data 
current_path<-getwd()

#run in Rstudio (add repository to paths at all locations in code): 
#zoneI
#SAN_path<-paste(current_path, "/scRNA-Seq-Variation-Pipeline/mouse_heart_GEO_data/SAN_GEO", sep="")    #path to SAN GEO data 

#zoneII
#AVN_path<-paste(current_path, "/scRNA-Seq-Variation-Pipeline/mouse_heart_GEO_data/AVN_GEO", sep="")    #path to AVN GEO data 

#zoneIII
#LPF_path<-paste(current_path, "/scRNA-Seq-Variation-Pipeline/mouse_heart_GEO_data/LPF_GEO", sep="")    #path to LPF GEO data 
#RPF_path<-paste(current_path, "/scRNA-Seq-Variation-Pipeline/mouse_heart_GEO_data/RPF_GEO", sep="")    #path to RPF GEO data 



#run in terminal: 
SAN_path<-paste(current_path, "/mouse_heart_GEO_data/SAN_GEO", sep="")    #path to SAN GEO data 
AVN_path<-paste(current_path, "/mouse_heart_GEO_data/AVN_GEO", sep="")    #path to AVN GEO data 
LPF_path<-paste(current_path, "/mouse_heart_GEO_data/LPF_GEO", sep="")    #path to LPF GEO data 
RPF_path<-paste(current_path, "/mouse_heart_GEO_data/RPF_GEO", sep="")    #path to RPF GEO data 


#load data form Cell Ranger from 3 different heart zones
SAN.data<-Read10X(SAN_path)
AVN.data<-Read10X(AVN_path)
LPF.data<-Read10X(LPF_path)
RPF.data<-Read10X(RPF_path)


#initialize the Seurat objects with the raw (non-normalized data)
zoneI<-CreateSeuratObject(counts = SAN.data, project = "zone I", min.cells = 3, min.features = 200)

zoneII<-CreateSeuratObject(counts = AVN.data, project = "zone II", min.cells = 3, min.features = 200)

zoneIIILPF<-CreateSeuratObject(counts = LPF.data, project = "zone IIILPF",min.cells = 3, min.features = 200)
zoneIIIRPF<-CreateSeuratObject(counts = RPF.data, project = "zone IIIRPF", min.cells = 3, min.features = 200)
zoneIII.combined <- merge(zoneIIILPF, y = zoneIIIRPF, add.cell.ids = c("LPF", "RPF"), project = "zone III")



### STANDARD PRE-PROCESSING WORKFLOW ###

#QC and selecting cells for further analysis
zoneI[["percent.mt"]] <- PercentageFeatureSet(zoneI, pattern = "^MT-")
zoneII[["percent.mt"]] <- PercentageFeatureSet(zoneII, pattern = "^MT-")
zoneIII.combined[["percent.mt"]] <- PercentageFeatureSet(zoneIII.combined, pattern = "^MT-")

#visualize QC metrics as a violin plot
VlnPlot(zoneI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zoneII, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zoneIII.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#use featureScatter to visualize feature-feature relationships, also can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(zoneI, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(zoneII, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneII, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(zoneIII.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneIII.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

zoneI <- subset(zoneI, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
zoneII <- subset(zoneI, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
zoneIII.combined <- subset(zoneIII.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



### NORMALIZING THE DATA ###

#normalize the data by employing a global-scaling normalization method “LogNormalize”
zoneI <- NormalizeData(zoneI, normalization.method = "LogNormalize", scale.factor = 10000)
zoneII <- NormalizeData(zoneII, normalization.method = "LogNormalize", scale.factor = 10000)
zoneIII.combined <- NormalizeData(zoneIII.combined, normalization.method = "LogNormalize", scale.factor = 10000)


### IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION) ###

#calculate a subset of features that exhibit high cell-to-cell variation in the dataset by directly modeling the mean-variance relationship inherent in single-cell data 
zoneI <- FindVariableFeatures(zoneI, selection.method = "vst", nfeatures = 2000)
zoneII <- FindVariableFeatures(zoneII, selection.method = "vst", nfeatures = 2000)
zoneIII.combined <- FindVariableFeatures(zoneIII.combined, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10zoneI <- head(VariableFeatures(zoneI), 10)
top10zoneII <- head(VariableFeatures(zoneII), 10)
top10zoneIII <- head(VariableFeatures(zoneIII.combined), 10)

#plot variable features with and without labels 
plot1 <- VariableFeaturePlot(zoneI)
plot2 <- LabelPoints(plot = plot1, points = top10zoneI, repel = FALSE)
plot1 + plot2

plot1 <- VariableFeaturePlot(zoneII)
plot2 <- LabelPoints(plot = plot1, points = top10zoneII, repel = FALSE)
plot1 + plot2

plot1 <- VariableFeaturePlot(zoneIII.combined)
plot2 <- LabelPoints(plot = plot1, points = top10zoneIII, repel = FALSE)
plot1 + plot2



###SCALING THE DATA###

#apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0, Scales the expression of each gene, so that the variance across cells is 1
#all.genes <- rownames(zoneI)
#all.genes <- rownames(zoneII)
#all.genes <- rownames(zoneIII.combined)


#include the 3 heart zones
all.genes <- rownames(zoneI)
zoneI <- ScaleData(zoneI, features = all.genes)

all.genes <- rownames(zoneII)
zoneII <- ScaleData(zoneII, features = all.genes)

all.genes <- rownames(zoneIII.combined)
zoneIII.combined <- ScaleData(zoneIII.combined, features = all.genes)


###  PERFORM LINEAR DIMENSIONAL REDUCTION ###

#perform PCA on the scaled data
zoneI <- RunPCA(zoneI, features = VariableFeatures(object = zoneI))
zoneII <- RunPCA(zoneII, features = VariableFeatures(object = zoneII))
zoneIII.combined <- RunPCA(zoneIII.combined, features = VariableFeatures(object = zoneIII.combined))

### DETERMINE THE 'DIMENSIONALITY' OF THE DATASET ###

#NOTE: This process can take a long time for big datasets, comment out for expediency. More
#approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time

zoneI <- JackStraw(zoneI, num.replicate = 100)
zoneI <- ScoreJackStraw(zoneI, dims = 1:15)
JackStrawPlot(zoneI, dims = 1:15)

zoneII <- JackStraw(zoneII, num.replicate = 100)
zoneII <- ScoreJackStraw(zoneII, dims = 1:15)
JackStrawPlot(zoneII, dims = 1:15)

zoneIII.combined <- JackStraw(zoneIII.combined, num.replicate = 100)
zoneIII.combined <- ScoreJackStraw(zoneIII.combined, dims = 1:15)
JackStrawPlot(zoneIII.combined, dims = 1:15)

#heuristic alternative to JackStraw  
ElbowPlot(zoneI)
ElbowPlot(zoneII)
ElbowPlot(zoneIII.combined)



###CLUSTER CELLS###
zoneI <- FindNeighbors(zoneI, dims = 1:14)
zoneI <- FindClusters(zoneI, resolution = 0.5)
tabI<-table(Idents(zoneI))


zoneII <- FindNeighbors(zoneII, dims = 1:14)
zoneII <- FindClusters(zoneII, resolution = 0.5)
tabII<-table(Idents(zoneII))


zoneIII.combined <- FindNeighbors(zoneIII.combined, dims = 1:14)
zoneIII.combined <- FindClusters(zoneIII.combined, resolution = 0.5)
tabIII<-table(Idents(zoneIII.combined))

#look at cluster IDs of the first 5 cells
head(Idents(zoneI), 5)
head(Idents(zoneII), 5)
head(Idents(zoneIII.combined), 5)


### RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/t-SNE) ###
#Run t-distributed Stochastic Neighbor Embedding
zoneI_tsne_path<-paste(current_path, "/seurat_output/zoneI_tsne.jpeg", sep="")    #path to zone I tsne
zoneII_tsne_path<-paste(current_path, "/seurat_output/zoneII_tsne.jpeg", sep="")    #path to zone II tsne
zoneIII_tsne_path<-paste(current_path, "/seurat_output/zoneIII_tsne.jpeg", sep="")    #path to zone III tsne


zoneI <- RunTSNE(zoneI,dims.use = 1:15,reduction.use = "pca")
jpeg(file = zoneI_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneI, reduction = "tsne")
dev.off()


zoneII <- RunTSNE(zoneII,dims.use = 1:15, reduction.use = "pca")
jpeg(file = zoneII_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneII, reduction = "tsne")
dev.off()

zoneIII.combined <- RunTSNE(zoneIII.combined,dims.use = 1:15, reduction.use = "pca")
jpeg(file = zoneIII_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneIII.combined, reduction = "tsne")
dev.off()


###FINDING DIFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS))###

#cluster numbers described in Goodyer et al paper
zoneI_cluster_path<-paste(current_path, "/seurat_output/zoneI_C9_HF.csv", sep="")    #path to cluster 9 variable ft 
zoneII_cluster_path<-paste(current_path, "/seurat_output/zoneII_C4_HF.csv", sep="")    #path to cluster 4 variable ft
zoneIII_cluster_path<-paste(current_path, "/seurat_output/zoneIII_C13_HF.csv", sep="")    #path to cluster 13 variable ft

#find all markers of Cluster 9 in zone I
cluster9.markers <- FindMarkers(zoneI, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)
write.csv(cluster9.markers, zoneI_cluster_path) #write deferentially expressed features to csv 


#find all markers of Cluster 4 in zone II
cluster4.markers <- FindMarkers(zoneII, ident.1 = 9, min.pct = 0.25)
head(cluster4.markers, n = 5)
write.csv(cluster4.markers, zoneII_cluster_path) #write deferentially expressed features to csv


#find all markers of Cluster 13 in Zone III
cluster13.markers <- FindMarkers(zoneIII.combined, ident.1 = 9, min.pct = 0.25)
head(cluster13.markers, n = 5)
write.csv(cluster13.markers, zoneIII_cluster_path) #write differentials expressed features to csv


