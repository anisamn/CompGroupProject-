#set up environment to retrieve, reformat, and map SRA reads
import retrieveDataSRA
import reformatDataSRA
import mapReads

#From NCBI GEO database
#GEO Accession -- CCS Component -- SRA run (raw data)
#GSM3885058 -- Zone I: SAN region -- SRR9290711
#GSM3885059	-- Zone II: AVN/His region -- SRR9290713
#GSM3885060	-- Zone III (Left): BB/PF region (Left Ventricle) -- SRR9290715
#GSM3885061	-- Zone III (Right): BB/PF region (Right Ventricle) -- SRR9290717

#created list of SRA sample runs (one of each)

SRRs = ['SRR9290711', 'SRR9290713', 'SRR9290715', 'SRR9290717']
print(SRRs) 


#retrieve mouse heart data from NCBI's SRA + retrieve mouse reference genome
retrieveDataSRA.getSRAdata(SRRs)
print('finished data retrieve') 


retrieveDataSRA.getRefGenome()
print('finished ref genome') 

#rename fastq files to Cell Ranger compatabile format
reformatDataSRA.renameFastqs(SRRs)
print('finished reformatting') 

#map sample reads using Cell Ranger and create Cell Ranger output folder
mapReads.runCellRanger(SRRs)
print('cell ranger run finished') 


