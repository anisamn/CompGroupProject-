import os

#function to reame fastq files to original format
def renameFastqs(SRRs):
    current_path = os.getcwd()                                               #change into SRA data folder directory
    folder = "/mouse_heart_SRA_data"
    os.chdir(current_path + folder)


    for i in range(0,len(SRRs)):
        #rename fastq files to cellranger compatitable format
        if(SRRs[i][-1] == '1'):
            rename_fastq1 = 'mv ' + SRRs[i] + '.1_1.fastq ' + 'SAN' + '_S' + str(i + 1) + '_L005' + '_R1' + '_001.fastq'
            rename_fastq2 = 'mv ' + SRRs[i] + '.1_2.fastq ' + 'SAN' + '_S' + str(i + 1) + '_L005' + '_R2' + '_001.fastq'
            rename_fastq3 = 'mv ' + SRRs[i] + '.1_3.fastq ' + 'SAN' + '_S' + str(i + 1) + '_L005' + '_I1' + '_001.fastq'

        elif(SRRs[i][-1] == '3'):
            rename_fastq1 = 'mv ' + SRRs[i] + '.1_1.fastq ' + 'AVN' + '_S' + str(i + 1) + '_L005' + '_R1' + '_001.fastq'
            rename_fastq2 = 'mv ' + SRRs[i] + '.1_2.fastq ' + 'AVN' + '_S' + str(i + 1) + '_L005' + '_R2' + '_001.fastq'
            rename_fastq3 = 'mv ' + SRRs[i] + '.1_3.fastq ' + 'AVN' + '_S' + str(i + 1) + '_L005' + '_I1' + '_001.fastq'

        elif(SRRs[i][-1] == '5'):
            rename_fastq1 = 'mv ' + SRRs[i] + '.1_1.fastq ' + 'LPF' + '_S' + str(i + 1) + '_L005' + '_R1' + '_001.fastq'
            rename_fastq2 = 'mv ' + SRRs[i] + '.1_2.fastq ' + 'LPF' + '_S' + str(i + 1) + '_L005' + '_R2' + '_001.fastq'
            rename_fastq3 = 'mv ' + SRRs[i] + '.1_3.fastq ' + 'LPF' + '_S' + str(i + 1) + '_L005' + '_I1' + '_001.fastq'

        elif(SRRs[i][-1] == '7'):
            rename_fastq1 = 'mv ' + SRRs[i] + '.1_1.fastq ' + 'RPF' + '_S' + str(i + 1) + '_L005' + '_R1' + '_001.fastq'
            rename_fastq2 = 'mv ' + SRRs[i] + '.1_2.fastq ' + 'RPF' + '_S' + str(i + 1) + '_L005' + '_R2' + '_001.fastq'
            rename_fastq3 = 'mv ' + SRRs[i] + '.1_3.fastq ' + 'RPF' + '_S' + str(i + 1) + '_L005' + '_I1' + '_001.fastq'

        os.system(rename_fastq1)
        os.system(rename_fastq2)
        os.system(rename_fastq3)

    os.chdir(current_path)                                                     #change to current directory

