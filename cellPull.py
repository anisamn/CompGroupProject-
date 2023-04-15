#This is the python function to get user input and run pullCellsAny.r
import os

#Set up pullCells function using the previously inputted gene of interest
def pullCells(interest_gene):
	#create and call command for Rscript and our gene of interest
	pull_cmd = 'Rscript pullCellsAny.R ' + interest_gene
	os.system(pull_cmd)
