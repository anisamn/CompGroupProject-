#load Seurat 
library(Seurat)
library(DropletUtils)

#Load gene of interest
gene_int <- commandArgs(trailingOnly = TRUE)

##Creating SAN data cell matrix##

#load the SAN dataset 


san.data <- Read10X(data.dir = 'mouse_heart_GEO_data/SAN_GEO')
dim(san.data) 

#Initialize the Seurat object with raw data 
#remove genes with less than three cells 
#remove cells with less than 200 genes expressed 
san <- CreateSeuratObject(counts = san.data, min.cells = 3, min.features = 200) 
san 


#create a subset with cells that have the Mybphl gene
####CHANGE Mybphl to gene_int
san.ourgene = subset(x = san, subset = Mybphl > 0)
san.ourgene 

#write out to CSV 
###CHANGE cellsWithMYBPHL_san to cellsWithGene_san
###Add paste as name_san
#name_san <- paste("cellsWith", gene_int, "_SAN.csv", sep="")
cellsWithGene_san <- GetAssayData(object = san.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithGene_san, "cellsWithMYBPHL_SAN.csv") 

#Add path as path_san
#path_san <- paste("cellsWith", gene_int, "_SAN.path", sep="")
#write10xCounts(x = san.ourgene@assays$RNA@counts, path = path_san)

#LPF data 
lpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/LPF_GEO') 
dim(lpf.data) 

lpf <- CreateSeuratObject(counts = lpf.data, min.cells = 3, min.features = 200)

lpf

### CHANGE Mybphl to gene_int
lpf.ourgene = subset(x = lpf, subset = Mybphl > 0) 
lpf.ourgene 

### CHANGE cellsWithMYBPHL_lpf to cellWithGene_lpf
### ADD paste as name_lpf
name_lpf <- paste("cellsWith", gene_int, "_LPF.csv", sep="")
cellsWithGene_lpf <- GetAssayData(object = lpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithGene_lpf, "cellsWithMYBPHL_LPF.csv") 

#Add path as path_lpf
#path_lpf <- paste("cellsWith", gene_int, "_LPF.path", sep="")
#write10xCounts(x = lpf.ourgene@assays$RNA@counts, path = path_lpf)

#RPF data 
rpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/RPF_GEO') 
dim(rpf.data) 

rpf <- CreateSeuratObject(counts = rpf.data, min.cells = 3, min.features = 200) 
rpf 

### CHANGE Mybphl to gene_int
rpf.ourgene = subset(x = rpf, subset = Mybphl > 0) 
rpf.ourgene 

### CHANGE cellsWithMYBPHL_rpf to cellWithGene_rpf
### ADD paste as name_rpf
name_rpf <- paste("cellsWith", gene_int, "_RPF.csv", sep="")
cellsWithGene_rpf <- GetAssayData(object = rpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithGene_rpf, "cellsWithMYBPHL_RPF.csv") 
#Add path as path_rpf
#path_rpf <- paste("cellsWith", gene_int, "_RPF.path", sep="")
#write10xCounts(x = rpf.ourgene@assays$RNA@counts, path = path_rpf)

#AVN data 
avn.data <- Read10X(data.dir = 'mouse_heart_GEO_data/AVN_GEO') 
dim(avn.data) 

avn <- CreateSeuratObject(counts = avn.data, min.cells = 3, min.featuers = 200)
avn 

### CHANGE Mybphl to gene_int
avn.ourgene = subset(x = avn, subset = Mybphl > 0) 
avn.ourgene 

### CHANGE cellsWithMYBPHL_avn to cellWithGene_avn
### ADD paste as name_avn
name_avn <- paste("cellsWith", gene_int, "_AVN.csv", sep="")
cellsWithGene_avn <- GetAssayData(object = avn.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithGene_avn, "cellsWithMYBPHL_AVN.csv") 
#Add path as path_avn
#path_avn <- paste("cellsWith", gene_int, "_AVN.path", sep="")
#write10xCounts(x = avn.ourgene@assays$RNA@counts, path = path_avn)

print('finished') 

