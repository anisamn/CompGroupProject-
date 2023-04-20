# ScRNA-seq mouse heart clustering - Group Project 

## About the project
![COMPgpWorkflow](https://user-images.githubusercontent.com/125702969/227393937-d72af62e-aec8-4467-913a-bcf7bc7c1e8c.png)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Single-cell RNA sequencing (scRNA-seq) is used for global gene expression analysis at a single cell resolution and to improve research on cell-type heterogeneity within the same tissue type. An understanding of cell heterogeneity is significant for developing targeted drugs and for the creation of new therapeutic approaches. It can also allow for a deeper understanding of the tissues and systems that use these cells. To identify the cells present, and map them to clusters, the pipeline scRNA-Seq Variation in Mouse Heart was created. The tool was created to retrieve fastq files from NCBI and map the reads pulled to a heart model created by CellRanger. The CellRanger output, combined with analysis by Seurat, can then be used to make plots and tables of significant clusters and differentially expressed features.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To add to the original scRNA-Seq Variation in Mouse Heart tool, a new pipeline will be added after the clustering stage of pipeline two that will create a list of highly expressed MYBPHL positive cells, identify a list of highly expressed genes from these cells, then re-cluster the list of highly expressed genes. Details on this can be seen in the flowchart below, and in the design document under the wiki tab.  

### References 
Goodyer, W. R., Beyersdorf, B. M., Paik, D. T., Tian, L., Li, G., Buikema, J. W., Chirikian, O., Choi, S., Venkatraman, S., Adams, E. L., Tessier-Lavigne, M., Wu, J. C., &amp; Wu, S. M. (2019). Transcriptomic profiling of the developing cardiac conduction system at single-cell resolution. Circulation Research, 125(4), 379–397. https://doi.org/10.1161/circresaha.118.314578 

## Workflow 

### Step one: Run pipeline one

This step pulls mouse heart scRNA-seq data from NCBI, organizes the FastQ files into a folder called mouse_heart_SRA_data, and stores genomic data into a folder called mouse_genome. A final output stores barcodes, features, and matrix files into a folder called cellranger_output. 

Scripts: 
1. retrieveDataSRA.py pulls mouse heart fastq data from NCBI, and a pre-computed reference genome from Cell Ranger. 
2. reformatDataSRA.py converts the fastq files to a Cell Ranger readable form 
3. mapReads.py maps sequencing reads to the reference genome using Cell Ranger 

Input files: 
- twelve SRA files, each associated with different regions of the heart
  - 3 AVN: Atrioventricular node/His bundle 
  - 3 LPF: Left purkinje fiber 
  - 3 RPF: Right purkinje fiber 
  - 3 SAN: Sinoatrial node 

Output files: 

- Three folders with mouse heart genomic data 

### Step two: Run pipeline two

This step pulls mouse heart GEO data from NCBI, and stores the data into a folder called mouse_heart_GEO_data. The new information is then reformatted into Seurat readable form, and outputs are stored in seurat_output. This is where Clustering.R is run to pass the reformatted data through Seurat. 

Scripts: 
1. retrieveDataGEO.py pulls mouse heart GEO data from NCBI containing associated barcodes, features, and matrix files. 
2. reformatDataGEO.py converts the barcodes, genes and matrix files to Seurat readable forms, and creates heart region associated output files. 
3. Clustering.R runs Seurat to cluster the cell data by performing QC, filtration, normalization, and PCA and dimensional reduction. 

Input files: 
- four folders from each region of the heart, each containing barcode, feature, and matrix information 

Output files: 
- seurat_output folder contianing matrices of cell and gene prevelance based on different heart zones, and three t-SNE plots mapping cells to the three zones of the heart. 

### Step three: Run pipeline three 
Pipeline three is the new parameter that can recluster cells based on a gene specified by the user. This outputs six files; three CSV files containing cells associated with the gene, and three t-SNE plots of the cells that contain the gene input by the user. 

Scripts: 
1. pipeline_3.py
2. cellpull.py 
3. pullCellsAny.R 

Input files: 
The output mapped CSV files in the seruat_output folder. 
- SAN_GEO
- LPF_GEO 
- RPF_GEO 
- AVN_GEO

## Example Run  
1. Clone the repository  
```
https://github.com/anisamn/CompGroupProject-
```
2. Run pipeline_2.py to download the required Mouse Heart scRNA-seq data. 
```
nohup python pipeline_2.py > pipeline2_run.out
```
3. Run pipeline three to see clusters of specified cells -- try entering 'MYBPHL to test 
```
python pipeline_3.py 
```
4. Outputs can be found in the file -- 


## Languages and Packages Required 

### Installation for CellRanger
CellRanger: analysis tool that processes single cell data, generates barcode matrices, and clustering

CellRanger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in

```
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1641366506&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDEzNjY1MDZ9fX1dfQ__&Signature=kaV8~ZabHhyDykUhbN~F78PDQfNZ64IamgsGc1nOSghFKPr0fbZ3WJk-2eWYh7IEt-KupenYP89W1zHi4lrxF~ZBbuP4NTaKEAa-G6ILJoX-VdyFnktkXFYDHgzEJ8ABq-NM6RWn20WD3a9BITNHTIWPtxjM-NaXAuR5uc5PuAEgjSDaQ2QBAQr~1q4aSM-~vJt~ia5e8acTz9RlM24EluLqfO59VCtAorP-5iJRwvLw9DjfrTlDtWfy3M2LSXp5OGmVJH1WUQReLK~0iZX2e8~vrHlAYpuxMa0Lgil6oHQ5s6vc~Dod3Aqpjb9sM~wuVo80zi4EqJ5nq0LU8SNbiQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

### Installation for Seurat

Seurat: An R package that performs QC and analysis of scRNA-seq data

Seurat: https://satijalab.org/seurat/articles/install.html

```
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
```

### Package Installation for R: 

```
library(Seurat) 
library(patchwork) 
library(dplyr) 
```

### Package Installation for Python: 

```
import os
import retrieveDataSRA 
import reformatDataSRA
import mapReads 
import retrieveDataGEO
import reformatDataGEO
import clusteringAnalysis

```
