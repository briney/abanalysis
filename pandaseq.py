#!/usr/bin/python
# filename: pandaseq.py

###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author:	Bryan Briney
# @props: 	PANDAseq team (github.com/neufeld/pandaseq)
# @cite: 	Andre P Masella, Andrea K Bartram, Jakub M Truszkowski, Daniel G Brown and Josh D Neufeld
#		 	PANDAseq: paired-end assembler for illumina sequences. 
#		 	BMC Bioinformatics 2012, 13:31.
#		 	www.biomedcentral.com/1471-2105/13/31
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################

import os
import glob
import subprocess as sp
from multiprocessing import cpu_count

def list_files(d):  
	return sorted([f for f in glob.glob(d + '/*') if os.path.isfile(f)])

def pair_files(files, nextseq):
	pairs = {}
	for f in files:
		if nextseq:
			f_prefix = '_'.join(os.path.basename(f).split('_')[:3])
		else:
			f_prefix = '_'.join(os.path.basename(f).split('_')[:2])
		if f_prefix in pairs:
			pairs[f_prefix].append(f)
		else:
			pairs[f_prefix] = [f,]
	return pairs

def batch_pandaseq(f, r, o):
	cmd = 'pandaseq -f {0} -r {1} -d rbfkms -T {3} -w {2}'.format(f, r, o, cpu_count())
	sp.Popen(cmd, shell=True, stderr=sp.STDOUT, stdout=sp.PIPE).communicate()

def merge_reads(files, output, nextseq, i):
	files.sort()
	f = files[0]
	r = files[1]
	if nextseq:
		lane = os.path.basename(f).split('_')[-3]
		sample_id = os.path.basename(f).split('_')[0]
		sample = sample_id + '_' + lane
	else:
		sample = os.path.basename(f).split('_')[0]
	print_sample_info(i, sample)
	o = os.path.join(output, '{}.fasta'.format(sample))
	batch_pandaseq(f, r, o)
	# print_sample_end()


def print_start_info():
	print ''
	print ''
	print '========================================'
	print 'Merging reads with PANDAseq'
	print '========================================'
	print ''

def print_input_info(files):
	print 'The input directory contains {} pair(s) of files to be merged.\n'.format(len(files) / 2)

def print_sample_info(i, sample):
	print '[ {} ]  Processing sample {}'.format(str(i), sample)

def print_sample_end():
	print 'Done.\n'


def run(input, output, nextseq):
	print_start_info()
	files = list_files(input)
	print_input_info(files)
	pairs = pair_files(files, nextseq)
	for i, pair in enumerate(sorted(pairs.keys())):
		if len(pairs[pair]) == 2: 
			merge_reads(pairs[pair], output, nextseq, i)

