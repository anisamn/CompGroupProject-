import os

#run Clustering. to normalize and cluster data

def performClustering():
    current_path = os.getcwd()                                  #create seurat results folder
    results_folder = "/seurat_output_GOI"
    os.mkdir(current_path + results_folder)


    seurat_cmd = 'Rscript Clustering_MYBPHL.R'
    os.system(seurat_cmd)

    print("Pipeline is finished running! Open 'seurat_ouput' folder to view significant results!")
    print(' _  _')
    print('(o)(o)--. ')
    print(' \../ (  )')
    print(' m\/m--m--')
