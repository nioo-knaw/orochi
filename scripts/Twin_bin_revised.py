'''Initial imports and print settings for numpy'''

from __future__ import division
import os
import argparse
import subprocess
import pickle
import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=200)
np.set_printoptions(threshold=np.nan)


def commandLineParser():
        """ Argparse settings (Command-line parameters and parser)
        Returns:
                dict: Dictionary containing all the parameter options and its supplied or default value.
        """
        parser = argparse.ArgumentParser(description='Parses diamond taxonomy and creates identifies bins')
	
	parser.add_argument('-i', '--input', type=str, required=True, help='Path to the folder that contains all binning tools output.')
	parser.add_argument('-dt', '--diamondTaxonomy', type=str, required=True, help='Path to diamond nr taxonomy output.')
        parser.add_argument('-o','--output', type=str, required=True, help='Directory where all output files will be stored.')
        parser.add_argument('-v','--verbose', required=False , default=False, action='store_true', help='Display additional output/information in STDOUT.')
        parser.add_argument('-T', '--tools', type=str, required=False, default='metabat,maxbin,concoct,groopm', help='comma separated name of the tools')
        return vars(parser.parse_args())


def megan(argOptions):
	'''Read diamond.nr-taxonomy.tsv file and parse its content to classify the bins
	each read classification is split in the different taxonomic levels, 
	all reads classifications for a given taxonomic level are grouped for each bin.
	If the percentage of reads belonging to a given taxonomic group is above 10%, 
	we include it in the classification group'''
	txt = '';	dic_pic = {}
	path = argOptions['diamondTaxonomy']
	In = open(path).read()
	In = In.split('\n')
	C = {};	O = {}; F = {}
	G = {};	S = {};	SP = {}
	print p
	for line in In:
		line = line.split('[')
		line = [x for x in line if x]
		c = [x for x in line if x.startswith('Class]')]
		o = [x for x in line if x.startswith('Order]')]
		f = [x for x in line if x.startswith('Family]')]
		g = [x for x in line if x.startswith('Genus]')]
		s = [x for x in line if x.startswith('Species]')]
		sp = [x for x in line if x.startswith('Species+]')]
		for x in c:
			x = x.split(';')
			if len(x) > 1:
				if x[0][7:] in C:
					C[x[0][7:]] += int(x[1])
				elif len(x)>=2:
					C[x[0][7:]] = int(x[1])
		for x in o:
			x = x.split(';')
			if len(x) > 1:
				if x[0][7:] in O:
					O[x[0][7:]] += int(x[1])
				elif len(x)>=2:
					O[x[0][7:]] = int(x[1])
		for x in f:
			x = x.split(';')
			if len(x) > 1:
				if x[0][8:] in F:
					F[x[0][8:]] += int(x[1])
				elif len(x)>=2:
					F[x[0][8:]] = int(x[1])
		for x in g:
			x = x.split(';')
			if len(x) > 1:
				if x[0][7:] in G:
					G[x[0][7:]] += int(x[1])
				elif len(x)>=2:
					G[x[0][7:]] = int(x[1])
		for x in s:
			x = x.split(';')
			if len(x) > 1:
				if x[0][9:] in S:
					S[x[0][9:]] += int(x[1])
				elif len(x)>=2:
					S[x[0][9:]] = int(x[1])
		for x in sp:
			x = x.split(';')
			if len(x) > 1:
				if x[0][10:] in SP:
					SP[x[0][10:]] += int(x[1])
				elif len(x)>=2:
					SP[x[0][10:]] = int(x[1])
	s = len(In)*100
#	s = sum([C[x] for x in C.keys()])
	for c in C.keys():
		if C[c]/s >= 0.1:
#			print 'Class : ' + c
#			print round(C[c]/s, 2)
			txt+= 'Class : ' + c + '\n' + str(round(C[c]/s, 2)) + '\n'
			if 'Class' not in dic_pic.keys():
				dic_pic['Class'] = [(c, round(C[c]/s, 2))]
			else:
				dic_pic['Class'].append((c, round(C[c]/s, 2)))

#	s = sum([O[x] for x in O.keys()])
	for o in O.keys():
		if O[o]/s >= 0.1:
#			print 'Order : ' + o
#			print round(O[o]/s, 2)
			txt+= 'Order : ' + o + '\n' + str(round(O[o]/s, 2)) + '\n'
			if 'Order' not in dic_pic.keys():
				dic_pic['Order'] = [(o, round(O[o]/s, 2))]
			else:
				dic_pic['Order'].append((o, round(O[o]/s, 2)))

#	s = sum([F[x] for x in F.keys()])
	for f in F.keys():
		if F[f]/s >= 0.1:
#			print 'Family : ' + f
#			print round(F[f]/s, 2)
			txt+= 'Family : ' + f + '\n' + str(round(F[f]/s, 2)) + '\n'
			if 'Family' not in dic_pic.keys():
				dic_pic['Family'] = [(f, round(F[f]/s, 2))]
			else:
				dic_pic['Family'].append((f, round(F[f]/s, 2)))

#	s = sum([G[x] for x in G.keys()])
	g_class = 'False'
	for g in G.keys():
		if G[g]/s >= 0.1:
			g_class = 'True'
#			print 'Genus : ' + g
#			print G[g]/s
			txt+= 'Genus : ' + g + '\n' + str(round(G[g]/s, 2)) + '\n'
			if 'Genus' not in dic_pic.keys():
				dic_pic['Genus'] = [(g, round(G[g]/s, 2))]
			else:
				dic_pic['Genus'].append((g, round(G[g]/s, 2)))

	if g_class == 'False':
		txt+= 'Genus : ' + 'Unclassified' + '\n' + str(round(G[g]/s, 2)) + '\n'
		if 'Genus' not in dic_pic.keys():
			dic_pic['Genus'] = [('Unclassified', round(G[g]/s, 2))]
		else:
			dic_pic['Genus'].append(('Unclassified', round(G[g]/s, 2)))

#	s = sum([S[x] for x in S.keys()])
	for s1 in S.keys():
		if S[s1]/s >= 0.1:
#			print 'Species : ' + s1
#			print S[s1]/s
			txt+= 'Species : ' + s1 + '\n' + str(round(S[s1]/s, 2)) + '\n'
			if 'Species' not in dic_pic.keys():
				dic_pic['Species'] = [(s1, round(S[s1]/s, 2))]
			else:
				dic_pic['Species'].append((s1, round(S[s1]/s, 2)))

#	s = sum([SP[x] for x in SP.keys()])
	for sp in SP.keys():
		if SP[sp]/s >= 0.1:
#			print 'Strain : ' + sp
#			print SP[sp]/s
			txt+= 'Strain : ' + sp + '\n' + str(round(SP[sp]/s, 2)) + '\n'
			if 'Strain' not in dic_pic.keys():
				dic_pic['Strain'] = [(sp, round(SP[sp]/s, 2))]
			else:
				dic_pic['Strain'].append((sp, round(SP[sp]/s, 2)))

	file_out = open('%s/classification.txt' %(argOptions['output']), 'w')
	file_out.write(txt)
	file_out.close()
	return dic_pic

def run_megan(B):
	'''megan command if the file is not found. [requires diamond output]'''
	cmd = "java -Xmx32G -Djava.awt.headless=true -Duser.language=en -Duser.region=US -cp '/data/tools/MEGAN/5.10.0/jars/MEGAN.jar:/data/tools/MEGAN/5.10.0/jars/data.jar' megan.tools.TaxonPathClassifier -i %s/diamond.tsv -f Detect -ms 50 -me 0.01 -tp 50 -g2t /data/db/megan/gi_taxid_prot-4March2015.bin -o %s/diamond.nr-taxonomy.tsv" %(B, B)
#	print cmd
	if not os.path.exists(argOptions['diamondTaxonomy']):
		subprocess.check_call(cmd, shell='True')

def classify(gp, bc):
	'''Unhash the print statements for additional information on the bins and their classification'''
	groups = {}
	for genus in set(gp):
		gb = []
		for k in bc.keys():
			if genus in bc[k]:
				gb.append(k)
		groups[genus] = gb
#		print genus
#		print len(gb)
#		print gb
	return groups

def build_matrix(g, argOptions):
	''''''
	mat = {}
	for genus in g.keys():
		mat[genus] = np.zeros((len(g[genus]), len(g[genus])))
		genus_contigs = set([])
		genome_contigs = {}
		for genome in g[genus]:
			genome_contigs[genome] = []
			bin_id = genome.split('/')[-1]
			tool_name = genome.split('/')[-2]
			F = open('%s/%s/%s' %(argOptions['input'],tool_name, bin_id)).read()
			F = F.split('>')
			for contig in F:
				contig = contig.split('\n')
				genus_contigs = genus_contigs | set([contig[0]])
				if contig[0]:
					genome_contigs[genome].append(contig[0])
		genus_contigs = [x for x in genus_contigs if x]
		x = -1;	y = -1
		for genome in g[genus]:
			x += 1
			for contig in genome_contigs[genome]:
				y = -1
				for genome2 in g[genus]:
					y += 1
					if contig in genome_contigs[genome2]:
						if int(contig.split('_')[3]) > 5000:
							mat[genus][x, y] += int(contig.split('_')[3])
	return mat


def percentage_matrix(mat):
	'''claculates the percentage identity between two bins'''
	for genus in mat.keys():
#		print genus
#		print mat[genus]
		for x in xrange(len(mat[genus])):
			for y in xrange(len(mat[genus][x])):
				if x != y:
					mat[genus][x][y] = mat[genus][x][y]/(min(mat[genus][x][x], mat[genus][y][y]))
			mat[genus][x] = np.around(mat[genus][x], decimals = 3)
#		print mat[genus]
	return mat

def find_twins(pmat, g):
	'''identifies twin bins [checks how many overlap more than 80% of their size]'''
	draft_genomes = {};	 draft_count = 0
	for genus in pmat.keys():
		print genus
		print pmat[genus]
		genus_set = set([])
		UC = 0
		for x in xrange(len(pmat[genus])):
			if len([n for n in pmat[genus][x] if n>0.8]) > 2:
				print pmat[genus][x]
				UC += 1
				print g[genus][x]
		print UC

def extract_draft(pmat, g):
	''''''
	D = []
	for genus in pmat.keys():
		print genus
		print g[genus]
		p_draft = set([x for x in xrange(len(pmat[genus])) if pmat[genus][x][x] > 750000 and pmat[genus][x][x] < 15000000])
		if len(pmat[genus])>1:
			for x in p_draft:
				draft = [X for X in p_draft if pmat[genus][x][X] > 0.5]
				if draft:
					p_draft = p_draft - set(draft)
					D.append([g[genus][n] for n in draft])
		else:
			print g[genus]
			D.append(g[genus])
			print len(D)
	return D
		
def organize_draft(D, argOptions):
	'''stores twin bins together'''
	counter = 0
	for twin in D:
		l = ''
		for x in twin:
			x = x.split('/')
			path = '%s/%s/%s' %(argOptions['input'],x[-2],x[-1])
			l+=path + ' '
		counter += 1
		if not os.path.exists('%sdraft_%s/' %(argOptions['output'], counter)):
			subprocess.check_call('mkdir %sdraft_%s/' %(argOptions['output'], counter), shell=True)
			cmd = 'cp %s %sdraft_%s/' %(l, argOptions['output'], counter)
			subprocess.check_call(cmd, shell=True)

if __name__ == '__main__':
	argOptions = commandLineParser()
	tools = argOptions['tools'].split(',')
	genus_pool = [];	bin_content = {}
	for f in tools:
		for b in os.listdir(f):
			if os.path.exists(f+b+'/diamond.tsv'):
				B = f+b
				run_megan(argOptions)
				if not os.path.exists('%s/classification.p' %(B)):
					dic_pic = megan(B, argOptions)
					pickle.dump(dic_pic, open('%s/classification.p' %(B), 'w'))
				pic = pickle.load(open('%s/classification.p' %(B)))
				if 'Genus' in pic.keys():
					genus_pool.extend([x[0] for x in pic['Genus']])
					bin_content[B] = [x[0] for x in pic['Genus']]
	g = classify(genus_pool, bin_content)
	mat = build_matrix(g, argOptions)
	pmat = percentage_matrix(mat)
	find_twins(mat, g)
	D = extract_draft(pmat, g)
	organize_draft(D, argOptions)
