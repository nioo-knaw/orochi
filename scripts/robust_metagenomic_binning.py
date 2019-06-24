'''
Generates a dictionary with the id of each contig
The dictionary holds the linkage score between the contigs.
Parse the different bins to record which contigs they contain
According to a scoring system for the different binning tools a different value is placed in correspondence of two contigs when they are found to be interacting
A soft threshold is used to remove the weak linkages
k-mer clustering is used to generate bins
Run checkm on the new formed bins
map the raw reads back to the contigs
pool the reads from a bin in new fastq files
rerun spades
rerun binning tools
rerun checkm
evaluate results
'''

import os
import numpy as np
import pickle
import argparse
import subprocess
import scipy.sparse

def commandLineParser():
	""" Argparse settings (Command-line parameters and parser)
	Returns:
		dict: Dictionary containing all the parameter options and its supplied or default value.
	"""
	parser = argparse.ArgumentParser(description='generates stats of a binner')
	
	# Files/directories/general options
	parser.add_argument('-i','--assembly', type=str, required=True, default='/home/NIOO/vittoriot/data/concoct/contigs/contigs.fasta', help='location of the spades output')
	parser.add_argument('-l','--logfile', type=str, required=True, default='/home/NIOO/vittoriot/output/maxbin2/checkm/log.txt', help='Path to the checkm log file.')
	parser.add_argument('-o','--outputDir', type=str, required=True, default='/home/NIOO/vittoriot/output/robust_binning/', help='Directory where all output files will be stored.')
	parser.add_argument('-v','--verbose', required=False , default=False, action='store_true', help='Display additional output/information in STDOUT.')
	parser.add_argument('-b','--bins', type=str, required=True, default='/home/NIOO/vittoriot/output/maxbin2/,/home/NIOO/vittoriot/data/concoct/bins/,/data/ngs/ME/raaijmakers_group/victorc/endophytic_meta/analysis/miplushi/binning/groopm/spades/,/data/ngs/ME/raaijmakers_group/victorc/endophytic_meta/analysis/miplushi/binning/metabat/spades/', help='Comma "," separated path to the binning tools output folders')
	parser.add_argument('-x','--extension', type=str, required=False, default='.fasta,.fa,.fna', help='Comma "," separated extensions of the fasta formats')
	return vars(parser.parse_args())

def get_files(argOptions):
	file_path = []
	folders = argOptions['bins'].split('|')
	for F in folders:
		f = os.listdir(F)
		f = [F+x for x in f if x.endswith(('.fasta','.fa','.fna'))]
		file_path.extend(f)

def create_2d_table(argOptions):
	F = open(argOptions['assembly']);	f = F.read()
	id_list = np.zeros(f.count('>'), dtype=np.int);	F.close()
	n = -1;	table = {}
	with open(argOptions['assembly']) as assembly:
		for line in assembly:
			if line.startswith('>'):
				n+=1
				id_list[n] = line.split('_')[-1]
				table[id_list[n]] = []
	return table

if __name__ == '__main__':
	argOptions = commandLineParser()
	get_files(argOptions)
	tab = create_2d_table(argOptions)
	
