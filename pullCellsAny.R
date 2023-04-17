#load Seurat 
#library(Seurat)

print("I have been called")
gene_int <- commandArgs(trailingOnly = TRUE)
test <- paste("test", gene_int, sep="")
print(test)

.f = function() {
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
san.ourgene = subset(x = san, subset = gene_int > 0)
san.ourgene 

#write out to CSV 
###CHANGE cellsWithMYBPHL_san to cellsWithGene_san
###Add paste
cellsWithGene_san <- paste("cellsWith", gene_int, "_san", sep="")
cellsWithMYBPHL_san <- GetAssayData(object = san.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_san, 'cellswithMYBPHL_SAN.csv') 

#LPF data 
lpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/LPF_GEO') 
dim(lpf.data) 

lpf <- CreateSeuratObject(counts = lpf.data, min.cells = 3, min.features = 200)

lpf

lpf.ourgene = subset(x = lpf, subset = Mybphl > 0) 
lpf.ourgene 

cellsWithMYBPHL_lpf <- GetAssayData(object = lpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_lpf, 'cellswithMYBPHL_LPF.csv') 

#RPF data 
rpf.data <- Read10X(data.dir = 'mouse_heart_GEO_data/RPF_GEO') 
dim(rpf.data) 

rpf <- CreateSeuratObject(counts = rpf.data, min.cells = 3, min.features = 200) 
rpf 

rpf.ourgene = subset(x = rpf, subset = Mybphl > 0) 
rpf.ourgene 

cellsWithMYBPHL_rpf <- GetAssayData(object = rpf.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_rpf, 'cellswithMYBPHL_RPF.csv') 

#AVN data 
avn.data <- Read10X(data.dir = 'mouse_heart_GEO_data/AVN_GEO') 
dim(avn.data) 

avn <- CreateSeuratObject(counts = avn.data, min.cells = 3, min.featuers = 200)
avn 

avn.ourgene = subset(x = avn, subset = Mybphl > 0) 
avn.ourgene 

cellsWithMYBPHL_avn <- GetAssayData(object = avn.ourgene, assay = 'RNA', slot = 'data') 

write.csv(cellsWithMYBPHL_avn, 'cellswithMYBPHL_AVN.csv') 

print('finished') 

}
print("This is the end of the Rscript")
