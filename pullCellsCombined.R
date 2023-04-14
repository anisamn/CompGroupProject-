#load Seurat 
library(Seurat)
library(DropletUtils)

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
san.ourgene = subset(x = san, subset = Mybphl > 0)
san.ourgene 

#write out to CSV 
cellsWithMYBPHL_san <- GetAssayData(object = san.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_san, 'cellswithMYBPHL_SAN.csv') 

write10xCounts(x = san.ourgene@assays$RNA@counts, path = 'cellswithMYBPHL_SAN.path')

#LPF data 
lpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/LPF_GEO') 
dim(lpf.data) 

lpf <- CreateSeuratObject(counts = lpf.data, min.cells = 3, min.features = 200)

lpf

lpf.ourgene = subset(x = lpf, subset = Mybphl > 0) 
lpf.ourgene 

cellsWithMYBPHL_lpf <- GetAssayData(object = lpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_lpf, 'cellswithMYBPHL_LPF.csv') 
write10xCounts(x = lpf.ourgene@assays$RNA@counts, path = 'cellswithMYBPHL_LPF.path')
#RPF data 
rpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/RPF_GEO') 
dim(rpf.data) 

rpf <- CreateSeuratObject(counts = rpf.data, min.cells = 3, min.features = 200) 
rpf 

rpf.ourgene = subset(x = rpf, subset = Mybphl > 0) 
rpf.ourgene 

cellsWithMYBPHL_rpf <- GetAssayData(object = rpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_rpf, 'cellswithMYBPHL_RPF.csv') 
write10xCounts(x = rpf.ourgene@assays$RNA@counts, path = 'cellswithMYBPHL_RPF.path')
#AVN data 
avn.data <- Read10X(data.dir = 'mouse_heart_GEO_data/AVN_GEO') 
dim(avn.data) 

avn <- CreateSeuratObject(counts = avn.data, min.cells = 3, min.featuers = 200)
avn 

avn.ourgene = subset(x = avn, subset = Mybphl > 0) 
avn.ourgene 

cellsWithMYBPHL_avn <- GetAssayData(object = avn.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_avn, 'cellswithMYBPHL_AVN.csv') 
write10xCounts(x = avn.ourgene@assays$RNA@counts, path = 'cellswithMYBPHL_AVN.path')
print('finished') 


