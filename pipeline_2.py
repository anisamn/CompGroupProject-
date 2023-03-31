#set up environment
import retrieveDataGEO
import reformatDataGEO
import clusteringAnalysis

#GEO Accession  -- CCS Component -- SRA run (raw data)
#GSM3885058     -- Zone I: SAN region -- SRR9290711
#GSM3885059	-- Zone II: AVN/His region -- SRR9290713
#GSM3885060	-- Zone III (Left): BB/PF region (Left Ventricle) -- SRR9290715
#GSM3885061	-- Zone III (Right): BB/PF region (Right Ventricle) -- SRR9290717


#retrieve mouse heart data from NCBI's GEO database *synonymous with Cell Ranger output*
retrieveDataGEO.getGEOdata()

#move barcode, genes, matrix files of each CCS component into designated folder (SAN,AVN,LPF,RPF)
reformatDataGEO.moveGEO()

#perform clustering and other statistical analyses using Seurat
clusteringAnalysis.performClustering()
