#!/usr/bin/python
# filename: batch_merge.py

###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @props: IgBLAST team (http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi)
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import os
import glob
import shutil
import argparse

import pandaseq

parser = argparse.ArgumentParser("Batch merging of paired-end reads with PANDAseq")
parser.add_argument('-i', '--in', dest='input', required=True, 
					help="The input directory, containing paired FASTQ files (uncompressed or gzip compressed). Required.")
parser.add_argument('-o', '--out', dest='output', required=True, 
					help="The output directory, will contain merged FASTA files. Required.")
parser.add_argument('-n', '--nextseq', dest='nextseq', default=False, action='store_true', 
					help="Use flag if run was performed on a NextSeq sequencer.")
args = parser.parse_args()


def make_direc(d):
	if not os.path.exists(d):
		os.mkdir(d)

def remove_direc(d):
	shutil.rmtree(d)

def list_files(d):  
	return sorted([f for f in glob.glob(d + '/*') if os.path.isfile(f)])

def bin_files(files):
	file_bins = {}
	for f in files:
		f_pre = '_'.join(os.path.basename(f).split('_')[:-1])
		if f_pre in file_bins:
			file_bins[f_pre].append(f)
		else:
			file_bins[f_pre] = [f,]
	return file_bins

def concat(d):
	files = list_files(d)
	file_bins = bin_files(files)
	for bin in file_bins:
		outfile = os.path.join(args.output, '{}.fasta'.format(bin))
		with open(outfile, 'w') as o:
			for f in file_bins[bin]:
				with open(f) as i:
					for line in i:
						o.write(line)

def main():
	make_direc(args.output)
	if args.nextseq:
		temp = os.path.join(args.output, 'temp')
		make_direc(temp)
		o = temp
	else:
		o = args.output
	pandaseq.run(args.input, o, args.nextseq)
	if args.nextseq:
		print ''
		print 'Concatenating NextSeq lane files for each sample...'
		concat(o)
		remove_direc(o)
		print 'Done.'
		print ''



if __name__ == '__main__':
	main()


