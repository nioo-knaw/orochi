import os
import argparse
import subprocess
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def commandLineParser():
	""" Argparse settings (Command-line parameters and parser)
	Returns:
		dict: Dictionary containing all the parameter options and its supplied or default value.
	"""
	parser = argparse.ArgumentParser(description='generates stats of a binner')
	
	# Files/directories/general options
	parser.add_argument('-l','--checkmlog', type=str, required=True, help='Path to the checkm log file.')
	parser.add_argument('-o','--outputDir', type=str, required=True, help='Path to the directory to store all files')
	parser.add_argument('-v','--verbose', required=False , default=False, action='store_true', help='Display additional output/information in STDOUT.')
	parser.add_argument('-i','--input', type=str, required=True, help='Folder where the bins are stored.')
	parser.add_argument('-x','--extension', type=str, required=False, default='.fasta', help='extension of the bins format, ".fasta" is default')
	
	return parser.parse_args()

class checkm:
	'''Parse the Checkm log file'''
	def __init__(self, argOptions):
		self.log_file = open(str(argOptions.checkmlog)).read()
		self.stats_dict = {}

	def stats(self):
		'''split file per row and column, store relevant information for each bin (row)'''
		rows = self.log_file.split('\n')[3:-2]
		self.bin_stats = [r.split(' ') for r in rows]
		for n in xrange(len(rows)):
			self.bin_stats[n] = [self.bin_stats[n][x] for x in xrange(len(self.bin_stats[n])) if self.bin_stats[n][x]]
			self.bin_stats[n][0] = ''.join(c for c in self.bin_stats[n][0] if c.isdigit())
			relevant = [1,2,12,13,14]
			self.stats_dict[self.bin_stats[n][0]] = [self.bin_stats[n][x] for x in relevant]
		print ('Base stats:')
		checkm.summary(self.stats_dict)

	def summary(self, d):
		'''Prints number of bins after filtering and other informative stats'''
		self.summary_stats = []
		for score in xrange(2,5):
			s = 0
			for k in d.keys():
				s+= float(d[k][score])
			self.summary_stats.append(s/len(d))
		print ('number of bins %d:' %(len(d)))
		print ('avg. completeness %d, avg contamination %d, avg heterogeneicity %d\n' %(self.summary_stats[0], self.summary_stats[1], self.summary_stats[2]))

	def filtering(self):
		soft = [50.00, 30.00];	self.s_dict = {}
		hard = [90.00, 5.00];	self.h_dict = {}
		for k in self.stats_dict.keys():
			if float(self.stats_dict[k][2]) > soft[0] and float(self.stats_dict[k][3]) < soft[1]:
				self.s_dict[k] = self.stats_dict[k]
			if float(self.stats_dict[k][2]) > hard[0] and float(self.stats_dict[k][3]) < hard[1]:
				self.h_dict[k] = self.stats_dict[k]
		
class bins:
	'''checks number of bins, bin size, number of contigs and n50 for each bin''' 
	def __init__(self,argOptions):
		self.folder = argOptions.input
		self.ext = argOptions.extension
		self.bin_stats = {}

	def base_stats(self):
		bins = os.listdir(self.folder)
		bins = [x for x in bins if x.endswith(self.ext)]
		for x in bins:
			k = ''.join(c for c in x if c.isdigit())
			stat = (os.stat(self.folder+x))
			self.bin_stats[k] = [stat.st_size/float(1000000)]
			f = open(self.folder+x).read()
			self.bin_stats[k].append(f.count('>'))
			self.bin_stats[k].append(round(self.bin_stats[k][0]*1000/f.count('>'),3))
		print (self.bin_stats)
		print 'number of binned contigs:'
		print (np.sum([self.bin_stats[x][1] for x in self.bin_stats.keys()]))
		print 'metagenome size [in Mb]:'
		print (np.sum([self.bin_stats[x][0] for x in self.bin_stats.keys()]))
	def plot(self, checkm):
		
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		xs = [self.bin_stats[x][0] for x in self.bin_stats.keys()]
		ys = [self.bin_stats[x][1] for x in self.bin_stats.keys()]
		zs = [self.bin_stats[x][2] for x in self.bin_stats.keys()]
		ax.scatter(xs, ys, zs, c='b', marker='o')

		ax.set_xlabel('Bin size')
		ax.set_ylabel('# of contigs')
		ax.set_zlabel('average contig size')

		plt.show()
		
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		xs = [float(checkm.stats_dict[x][2]) for x in checkm.stats_dict.keys()]
		ys = [float(checkm.stats_dict[x][3]) for x in checkm.stats_dict.keys()]
		zs = [float(checkm.stats_dict[x][4]) for x in checkm.stats_dict.keys()]
		ax.scatter(xs, ys, zs, c='b', marker='o')

		ax.set_xlabel('completeness')
		ax.set_ylabel('contamination')
		ax.set_zlabel('heterogeneicity')

		plt.show()

if __name__ == '__main__':
	argOptions = commandLineParser()
	checkm = checkm(argOptions)
	checkm.stats()
	checkm.filtering()
	print ('soft filter stats:')
	checkm.summary(checkm.s_dict)
	print ('hard filter stats:')
	checkm.summary(checkm.h_dict)
	bins = bins(argOptions)
	bins.base_stats()
	bins.plot(checkm)
