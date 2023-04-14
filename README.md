# CompGroupProject-##Languages and Packages
## Languages and Packages

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

# Flags for input files: 

### 
