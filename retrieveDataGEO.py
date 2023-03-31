import os


#function to retrieve mouse heart data from NCBI's GEO database
#*synonymous with Cell Ranger output*

def getGEOdata():
    current_path = os.getcwd()                                               #create and change to GEO data folder directory
    folder = "/mouse_heart_GEO_data"

    os.mkdir(current_path + folder)
    os.chdir(current_path + folder)

    wget_GEO = "wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132658/suppl/GSE132658_RAW.tar'"
    decompress_GEO = "tar -xvf GSE132658_RAW.tar"
    os.system(wget_GEO)
    os.system(decompress_GEO)

    os.chdir(current_path)                                                     #change to current directory
