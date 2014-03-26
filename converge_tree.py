#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Test for convergent signals in bipartitions of trees.

Requires: 
   A data file with taxon<tab>value
   A tree file with matching taxon names.

Also: dendropy, which can be installed with:
   sudo pip install dendropy

Version 0.1 - First port of PERL version by Joseph F. Ryan
"""

import sys
from optparse import OptionParser


def main():
	Options,Args = get_options()

	if not(Options.DataFileName and Options.TreeFileName):
		sys.stderr.write("You must specify a data file [-d] (taxon<tab>trait) and Newick tree file [-t].\n")
	else:

		# print rh_depths
	
		rh_depths = get_depths(Options.DataFileName)
		ra_taxa, ra_clades = get_taxa_and_clades(Options.TreeFileName,rh_depths)
		ra_max_scores,depthlist = max_scores(rh_depths,ra_taxa)
	
		scorelist,ra_scores = get_score(ra_taxa,ra_clades,rh_depths)
		# for j in ra_scores:
		# 	print j
			
		diff, best_score, best_count = get_best_score_and_count(scorelist,ra_scores)
		max_possible = ra_max_scores[best_count]
		standev = stdev(depthlist)
	
		if standev:
			num_stdevs = round((max_possible - best_score) / standev,1);
		else:
			num_stdevs = 0
	
		Header = "TREEFILE\tSCORE\tMAX_POSSIBLE\tNUM_STDEVS\tSTDEV\tCOUNT_DIFF_IN_TWO_CLADES\tTAXA_IN_CLADE\tTOTAL_TAXA"
		if Options.PrintHead:
			print Header
		print "{0}\t{1:.1f}\t{2:.1f}\t{3:.1f}\t{4:.2f}\t{5}\t{6}\t{7}".format(Options.TreeFileName,
		        best_score,max_possible,num_stdevs,standev,diff,best_count,len(ra_taxa))


def mean(a):
	return float(sum(a))/len(a)
	
def stdev(a):
	# note: added -1 to the divisor
	s2 = sum( (x - mean(a))**2  for x in a )/ (len(a)-1)
	return s2**0.5
	
def get_options():
	parser = OptionParser(usage = __doc__)
#	parser.add_option("-v", "--version", action="store_true", dest="version", default=False, help="Print program version")
	parser.add_option("-d", "--datafile", action="store", type="string", dest="DataFileName", help="Path to data file")
	parser.add_option("-p", "--printheader", action="store_true", dest="PrintHead", default=False, help="Print header line before output")
	parser.add_option("-t", "--treefile", action="store", type="string", dest="TreeFileName", help="Path to tree file")
	(options, args) = parser.parse_args()
	return options, args


def get_depths(Fname):
	"""Read in dictionary of depths"""
	DD = {}
	with open(Fname,'rU') as Dfile:
		for Line in Dfile:
			if Line:
				T,V = Line.rstrip("\n").split()
				if V.replace('.','').isdigit():
					DD[T] = float(V)
	return DD
	
def get_taxa_and_clades(Tname,dep):
	import re
	"""Parse tree file and depths"""
	taxa = []
	clades = []
	clearmatch = r':[e\d\-\.]+'
	nodematch = r'[,\(\)]'
	NodeList = sorted(get_tree_lines(Tname))
	for N in NodeList:
		Node = re.sub(clearmatch,'',N) 
		# look for parentheses-free lines (single taxa)
		if not re.search(nodematch,Node):
			# if the node is a key in the depth dictionary
			if Node in dep:
				taxa.append(Node)
			else:
				sys.exit("No data found for {}".format(Node))
		# pull clade names out of multi-element lists
		else:
			# print "Node:",Node
			Node = Node.replace('(',"").replace(')','')
			Fields = Node.split(',')
			# print "Fields:",Fields
			clades.append(Fields)
	return taxa,clades
	
def get_tree_lines(Tname):
	stringlist =[]
	from dendropy import Tree
	tree = Tree.get_from_path(Tname,"newick")
	for nd in tree.postorder_internal_node_iter():
	    for child in nd.child_nodes():
	        stringlist.append(child.as_newick_string())
	return (stringlist)
	
	
def max_scores(dep,taxa):
	mscore = {}
	dv = [dep[d] for d in dep.keys() if d in taxa]
	sdepth = sorted(dv)[::-1]
	for i in range(1,len(sdepth)+1):
		mscore[i] = sum(sdepth[:i])
	return mscore, dv
	

def get_score(tax,clades,dep):
	low_diff = 1e9
	fullscore=[]
	scores = []
	labels = ['score','s_a','s_b','c_a','c_b','diff']
	for clad in clades:
		score_a = [ dep[ta] for ta in clad ]
		score_b = [ dep[tb] for tb in tax if tb not in clad]
		count_a = len(score_a)
		count_b = len(score_b)
		sum_a   = sum(score_a)
		sum_b   = sum(score_b)
		c_diff  = abs(count_a - count_b)
		d_diff  = abs(sum_a - sum_b)
		if (c_diff <= low_diff):
			low_diff = c_diff
			fullscore.append(dict(zip(labels,
			          (d_diff, sum_a, sum_b, count_a, count_b, c_diff))))
			scores.append(c_diff)
	return scores,fullscore

def get_best_score_and_count(scorelist,ra_scores):
	# print scorelist
	sortsc = [x for (y,x) in sorted(zip(scorelist,ra_scores))]
	# print sortsc
	diff = sortsc[0]['diff']
	score = 0
	for sc in sortsc:
		if sc['diff'] == diff:
			if sc['score'] > score:
				# use logical 0 or 1 to index the value
				AorB = ['b','a'][(sc['s_a'] > sc['s_b'])]
				best_sc = sc['s_'+ AorB] 
				best_ct = sc['c_'+ AorB] 
		else:
			break
	return sortsc[0]['diff'],best_sc, best_ct


def get_best(scorelist,ra_scores):
	position = scorelist.index(max(scorelist))
	# use logical 0 or 1 to index the value
	AorB = ['b','a'][(ra_scores[position]['s_a'] > ra_scores[position]['s_b'])]
	
	best_sc    = ra_scores[position]['s_'+ AorB] 
	best_ct = ra_scores[position]['c_'+ AorB] 
	return ra_scores[position]['diff'],best_sc, best_ct
	

if __name__ == "__main__":
	main()
