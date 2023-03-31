import os


#function to retrieve mouse heart data from NCBI's SRA database
def getSRAdata(SRRs):
    current_path = os.getcwd()   #create and change to SRA data folder directory
    folder = "/mouse_heart_SRA_data"
    os.mkdir(current_path + folder)
    os.chdir(current_path + folder)

    for SRR in SRRs:
        if(SRR[-1] == '3'):
            SRR_data_address = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/' + SRR + '/' + SRR + '.1'
            wget_SRR = 'wget' + ' ' + SRR_data_address
            fastq_dump_SRR = 'fastq-dump -I --split-files' + ' ' + SRR + '.1'
            os.system(wget_SRR)  #retrieve SRA data using wget command
            os.system(fastq_dump_SRR) #uncompress data and convert to paired-end fastq files using fastq-dump command


        else:
            SRR_data_address = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/' + SRR + '/' + SRR + '.1'
            wget_SRR = 'wget' + ' ' + SRR_data_address
            fastq_dump_SRR = 'fastq-dump -I --split-files' + ' ' + SRR + '.1'
            os.system(wget_SRR)   #retrieve SRA data using wget command
            os.system(fastq_dump_SRR)  #uncompress data and convert to paired-end fastq files using fastq-dump command

    os.chdir(current_path) #change to current directory


#function to retrieve mouse reference genome index
def getRefGenome():
    current_path = os.getcwd()  #create and change to mouse genome folder directory
    folder = "/mouse_genome"
    print('current path: ',  current_path + '/mouse_genome') 
    
    os.mkdir(current_path + folder)
    
    os.chdir(current_path + folder)

    wget_genome = 'wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz'  #retrieve cell ranger mouse genome
    decompress_genome = 'tar -zxvf refdata-gex-mm10-2020-A.tar.gz'
    os.system(wget_genome)
    os.system(decompress_genome)

    os.chdir(current_path) #change to current directory



