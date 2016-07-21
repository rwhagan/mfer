#!/usr/local/bin/python3

'''This script formats and transposes OTU tables, merges them with metadata files, and performs categorical analyses. Usage: mfer.py -t otu_table -m metadata_file -n categorical_variable'''
__author__ = "Richard Hagan, David Jacobson, Allison Mann, and Krithivasan Sankaranarayanan"
__license__ = "GPL"
__version__ = "1.0"

import pandas as pd
from pandas import DataFrame
import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='This script formats and transposes OTU tables, merges them with metadata files, and performs categorical analyses. Usage: mfer.py -t otu_table -m metadata_file -n categorical_variable')

#Required arguments
parser.add_argument('-t', '--table', help="The table you wish to format.", required=True)
parser.add_argument('-m', '--metadata', help="The metadata file you wish to merge", required=True) 

#Optional arguments
parser.add_argument('-n', '--categoryName', help='Metadata category used to compare taxa significance. Default is all categories')
parser.add_argument('-c', '--collapsedOut', help='Desired collapsed OTU table name, default is collapsed_otu_table.txt', default='collapsed_otu_table.txt')
parser.add_argument('-o', '--compareOut', help='Output file name for cateogry compare step, default is compareOut.txt', default='compareOut.txt')
parser.add_argument('-v', '--mergedOut', help='Output file name for merged table step, default is mergedOut.txt', default='mergedOut.txt')
parser.add_argument('-r', '--transpose', help='Output file name for transpose table step, default is tranposedTable.txt', default='transposedTable.txt')
parser.add_argument('--fdr', help='Optional. Choose fdr value for significance, default is 0.1',default=0.1)
parser.add_argument('--pval', help='Optional. Choose p-value for significant, default is 0.05',default=0.05)
parser.add_argument('-s', '--significanceOutput', help = 'the output from the R script, filtered for significant pval and fdr' , default = 'significanceOutput.txt')

args = parser.parse_args()

class colors:
	RUN = '\033[94m'
	COMPLETE = '\033[92m'
	ENDC = '\033[0m'
	FAIL = '\033[91m'

def prelimCheck():
	assert os.path.exists(args.table), colors.FAIL + 'ERROR: File does not exist: %s. Is the path correct?' % args.table
	assert os.path.exists(args.metadata), colors.FAIL + 'ERROR: File does not exist: %s. Is the path correct?' % args.metadata
	if os.path.isfile('Compare.R'):
		print('Running...')
	else:
		print(colors.FAIL + 'ERROR: R script not found, is it in your current working directory?' + colors.ENDC)


def tableCheck():
	global cleanedTable
	print(colors.RUN + 'Checking OTU table format...' + colors.ENDC)
	with open(args.table, 'r') as f:
		first = f.readline() #pull out only first line of file
		if '# ' in first: #check to see if leading header is there
			inTable = pd.read_csv(args.table, sep='\t', skiprows=1) #if so, skip
		else:
			inTable = pd.read_csv(args.table, sep='\t') #else read in normally

	inTable.replace(to_replace=' *', value='', regex=True, inplace=True) #remove spaces from taxa strings
	inTable.columns = [x.strip().replace(' ', '_') for x in inTable.columns] #replace any spaces from column names and change to underscores

	cleanedTable = inTable.groupby('#OTU_ID').sum() #group by taxonomy string, get sum across all rows
	print(colors.COMPLETE + "Number of duplicate taxa strings to be collapsed: %i" % (inTable.shape[0] - cleanedTable.shape[0]+1) + colors.ENDC) #number of taxa collapsed, added one to account for differences in shape	
	
	with open(args.collapsedOut, 'w') as outfile:
		cleanedTable.to_csv(outfile, sep='\t')
	outfile.close()

	print(colors.RUN + "Table checking complete. Collapsed OTU table written to: %s..." % args.collapsedOut + colors.ENDC)


def transpose(cleanedTable): #the input to this function should be Allie's cleaned table
	global transposedTable
	transposedTable = cleanedTable.transpose()
	transposedTable.to_csv(args.transpose, sep='\t',)
	print(colors.RUN + "Table transposed. Written to %s..." % args.transpose + colors.ENDC)


def merging():
	global metadataHeaders
	print(colors.RUN + 'Merging %s and %s...' % (args.metadata, args.transpose) + colors.ENDC)
	with open(args.metadata, 'r') as m:
		metadataHeaders = list(m.readline().split())[1:]
	metaDF = pd.read_csv(args.metadata, sep='\t')
	taxaDF = pd.read_csv(args.transpose, sep = '\t')
	mergedDF = pd.merge(metaDF, taxaDF, left_index = True, right_index = True, how = 'inner') ## merge the metadata and transposed files using the first header in each file as the criteria
	### should we do merge inner or outer? 
	del mergedDF['Unnamed: 0'] ### delete the column that has sample names from taxa table
	mergedDF.columns = [x.strip().replace('#', '') for x in mergedDF.columns] #remove hash symbol	
	mergedDF.to_csv(args.mergedOut, sep = '\t', index = False)
	
def Rcomp():
	if args.categoryName:
		print(colors.RUN + 'Comparing groups, split by %s..' % args.categoryName + colors.ENDC)
		subprocess.call(args=['Rscript Compare.R %s %s %s' %(args.mergedOut, args.categoryName, args.compareOut)], shell=True)
		print(colors.COMPLETE + 'Complete! Output written to %s' % args.compareOut + colors.ENDC)
	else:
		for category in metadataHeaders:
			print(colors.RUN + 'Comparing groups, split by %s..' % category + colors.ENDC)
			subprocess.call(args=['Rscript Compare.R %s %s %s' %(args.mergedOut, category, category + "compared.txt")], shell=True)
			print(colors.COMPLETE + 'Complete! Output written to %s' % (category + "compared.txt") + colors.ENDC)	
	

def filterR():
	df =  pd.read_csv(args.compareOut, sep = '\t')
	userPval = float(args.pval)
	userFDR = float(args.fdr)
	newDF = df.ix[(df['pval']<=userPval) & (df['fdr']<=userFDR)]
	newDF.to_csv(args.significanceOutput, sep = '\t', index = False)


prelimCheck()
tableCheck()
transpose(cleanedTable)
merging()
Rcomp()
filterR()




