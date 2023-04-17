#This is the outline for Pipeline 3
#Set up environment
import cellPull
#Ask for user input
interest_gene = input("Enter gene of interest: ")
#Pull all cells of indicated the indicated gene
cellPull.pullCells(interest_gene)
