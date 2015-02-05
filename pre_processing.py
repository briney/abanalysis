#!/usr/bin/python
# filename: uaid.python

###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import os
from Bio import SeqIO
from multiprocessing import Pool, cpu_count


def process_raw_seqs(f, out_dir):
	out_file = os.path.join(out_dir, os.path.basename(f))
	open(out_file, 'w').write('')
	out_handle = open(out_file, 'a')
	for s in SeqIO.parse(open(f, 'r'), 'fasta'):
		out_handle.write('>{0}_{1}_\n{1}\n'.format(s.id, str(s.seq)))

def print_start_info():
	print ''
	print ''
	print '========================================'
	print "Pre-processing samples"
	print '========================================'
	print ''

def print_input_info(files):
	print 'Pre-processing sequences from {0} files to retain raw sequence data after IgBLAST...'.format(len(files))

def print_sample_end():
	print 'Done.\n'

def run(files, output):
	print_start_info()
	print_input_info(files)
	p = Pool(processes=cpu_count())
	for f in files:
		p.apply_async(process_raw_seqs, args=(f, output))
	p.close()
	p.join()


