#!/usr/bin/python
# filename: ab_analysis.py

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
import time
import math
import glob
import platform
import argparse
import threading
from subprocess import Popen, PIPE
from multiprocessing import Pool, cpu_count

from Bio import SeqIO

from blast_parse import BlastParse


parser = argparse.ArgumentParser("Antibody annotation with IgBLAST.")
parser.add_argument('-i', '--in', dest='input', required=True, 
					help="The input file, to be split and processed in parallel. \
					If a directory is given, all files in the directory will be iteratively processed.")
parser.add_argument('-o', '--out', dest='output', required=True, 
					help="The output directory, which will contain JSON or tab-delimited output files.")
parser.add_argument('-l', '--log', dest='log', default='', 
					help="The log file, to which the BlastParse log info will be written. \
					Default is stdout.")
parser.add_argument('-t', '--temp', dest='temp_dir', default='', 
					help="The directory in which temp files will be stored.  \
					If the directory doesn't exist, it will be created.  \
					Defaults to './temp_files'.")
parser.add_argument('-p', '--threads', dest='num_threads', default=0, type=int, 
					help="Number of parallel igblastn instances to spawn. \
					Defaults to max available processors.")
parser.add_argument('-v', '--tsv', dest="tsv_out", action='store_true', default=False, 
					help="NOT YET IMPLEMENTED. If set, the final output (from BlastParse) will be in tab-delimited format.  \
					Defaults to JSON output.")
parser.add_argument('-m', '--merge', dest="merge", action='store_true', default=False, 
					help="Use if the input files are paired-end FASTQs (either gzip compressed or uncompressed) \
					from Illumina platforms. Prior to running IgBLAST, reads will be merged with pandaseq. \
					Requires that pandaseq is installed.")
parser.add_argument('-n', '--next_seq', dest="next_seq", action='store_true', default=False, 
					help="Use if the run was performed on a NextSeq sequencer. \
					Multiple lane files for the same sample will be merged.")
parser.add_argument('-u', '--uaid', dest="uaid", type=int, default=None, 
					help="Use if the input files contain unique antibody identifiers (UAIDs). \
					UAIDs will be identified and incorporated into the output JSON file.")
parser.add_argument('-b', '--basespace', dest="use_basespace", default=False, action='store_true', 
					help="NOT YET IMPLEMENTED. Use flag if files should be downloaded directly from BaseSpace. \
					Files will be downloaded into the directory provided with the '-i' flag, which should be empty.")
parser.add_argument('-d', '--debug', dest="debug", action='store_true', default=False, 
					help="If set, will write all failed/exception sequences to file and give more informative errors.")
parser.add_argument('-s', '--species', dest='species', default='human', 
					choices=['human', 'macaque', 'mouse'])
args = parser.parse_args()


class launch_thread(threading.Thread):

	def __init__(self, in_file, out_file):
		threading.Thread.__init__(self)
		self.in_file = in_file
		self.out_file = out_file
		binary = './igblastn_' + platform.system().lower()
		self.cmd = '{3} -germline_db_V database/{0}_gl_V -germline_db_J database/{0}_gl_J -germline_db_D database/{0}_gl_D ' + \
				   '-organism {0} -domain_system imgt -auxiliary_data optional_file/{0}_gl.aux ' + \
				   '-show_translation -outfmt 3 -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 ' + \
				   '-query {1} -out {2}'.format(args.species, self.in_file, self.out_file, binary)

	def run(self):
		p = Popen(self.cmd, shell=True, stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
		


#####################################################################
# 
#                        FILES AND DIRECTORIES
#
#####################################################################

def build_temp_dirs():
	if args.temp_dir != '': 
		temp_directory = args.temp_dir
	else: 
		temp_directory = "./temp_files"
	temp_out_directory = temp_directory + "/temp_output"
	if not os.path.exists(temp_directory): os.mkdir(temp_directory)
	if not os.path.exists(temp_out_directory): os.mkdir(temp_out_directory)
	return temp_directory, temp_out_directory

def build_output_dir():
	output_dir = args.output
	if not os.path.exists(output_dir): os.mkdir(output_dir)
	return output_dir

def list_files(d): 
	if os.path.isdir(d):
		expanded_dir = os.path.expanduser(d)
		return sorted(glob.glob(expanded_dir + '/*'))
	else:
		return [d,]

def file_length(file):
	c = 0
	with open(file) as f:
		for i, l in enumerate(f):
			if l[0] == '>': c += 1
	return c

def num_procs():
	if args.num_threads > 0: 
		return args.num_threads
	return cpu_count()

def clean_up(temp_files, temp_out_files, temp_directory, temp_out_directory):
	for file in temp_out_files:
		os.remove(file)
	os.rmdir(temp_out_directory)
	for file in temp_files:
		os.remove(file)
	os.rmdir(temp_directory)

def file_splitter(file, splitlen, num_seqs, temp_directory):
	counter = 1
	files_list = []
	lines = open(file, 'r').read().replace(' ', '_').split('>')
	for line in range(0, num_seqs+1, splitlen):
		output = lines[line:line+splitlen]
		temp_filename = temp_directory + "/tempfile_" + str(counter)
		files_list.append(temp_filename)
		open(temp_filename, "w").write("")
		temp_file = open(temp_directory + "/tempfile_" + str(counter), "a")
		if counter == 1:
			temp_file.write('>'.join(output))
			counter += 1
		else:
			temp_file.write('>' + '>'.join(output))
			counter += 1
		temp_file.close()
	return files_list
		


#####################################################################
# 
#                            PRINTING
#
#####################################################################

def print_input_info(i):
	print ''
	print ''
	print '========================================'
	print 'Parallel IgBLAST'
	print '========================================'
	print ''	
	if len(i) > 1:
		print 'Input is a directory of {} files.\n\n'.format(len(i))
	else:
		print 'Input is a single file.\n\n'

def print_infile(i):
	b = os.path.basename(i)
	print '-'*len(b)
	print b
	print '-'*len(b)

def print_summary_output(g, e, f, blast_time, parse_time):
	total_seqs = g+e+f
	print ''
	print 'Out of {} total sequences:'.format(total_seqs)
	print '{} sequences processed normally'.format(g)
	print '{} sequences passed sanity checks, but could not be processed'.format(e)
	print '{} sequences failed sanity checks are were not processed'.format(f)
	print ''
	print 'IgBLAST took {0} seconds ({1} sequences per second)'.format(blast_time, total_seqs/blast_time)
	print 'parsing took {0} seconds ({1} sequences per second)'.format(parse_time, total_seqs/parse_time)
	print ''



#####################################################################
# 
#                              PARSING
#
#####################################################################

def line_generator(blast_file):
	f = open(blast_file, 'r')
	for line in f:
		yield line

def block_generator(blast_file):
	l = line_generator(blast_file)
	line = next(l)
	while line.find('Query= ') == -1: line = next(l)
	block = line.replace('Query= ', '')
	while True:
		try:
			line = next(l)
			while line.find('Query= ') == -1: 
				block += line
				line = next(l)
			yield block
			block = line.replace('Query= ', '')
		except StopIteration:
			yield block
			break
	raise StopIteration

def do_parse(blastout):
	out_file = blastout + '.json'
	result = []
	pool = Pool(processes=cpu_count())
	for i in block_generator(blastout):
		try:
			if args.debug:
				result.append(parser(i))
			else:
				result.append(pool.apply_async(parser, (i,)))
		except StopIteration:
			break
	pool.close()
	pool.join()
	good, exceptions, failed = process_parse_data(result, out_file)
	result = []
	return good, exceptions, failed

def parser(i):
	bp = BlastParse(i, species=args.species, tsv=args.tsv_out, log=args.log, debug=args.debug, uaid=args.uaid)
	if bp.sanity_checks() < 1:
		output = bp.parse()
		return output
	else:
		return ['', '', i]

def process_parse_data(results, out_file):
	good = 0
	exceptions = 0
	failed = 0
	r_handle = build_result_handle(out_file)
	if args.debug:
		e_handle = build_exception_handle(out_file)
		f_handle = build_failed_handle(out_file)
	for result in results:
		if args.debug:
			if result[0] != '':
				r_handle.write(result[0])
				good += 1
			elif result[1] != '':
				e_handle.write(result[1])
				exceptions += 1
			elif result[2] != '':
				f_handle.write(result[2])
				failed += 1
		else:
			r = result.get()
			if r[0] != '':
				r_handle.write(r[0])
				good += 1
			elif r[1] != '': exceptions += 1
			elif r[2] != '': failed += 1
	return good, exceptions, failed

def build_result_handle(out_file):
	open(out_file, 'w').write('')
	return open(out_file, 'a')

def build_exception_handle(out_file):
	e_file = out_file.split('.')[0] + '_exceptions'
	open(e_file, 'w').write('')
	return open(e_file, 'a')

def build_failed_handle(out_file):
	f_file = out_file.split('.')[0] + '_failed'
	open(f_file, 'w').write('')
	return open(f_file, 'a')



#####################################################################
# 
#                       INPUT PROCESSING
#
#####################################################################

def check_input(input_list):
	format = format_check(input_list[0])
	if format == 'fasta':
		return input_list
	else:
		return convert_to_fasta(input_list)

def format_check(in_file):
	with open(in_file) as f:
		line = f.next()
		while line == '':
			line = f.next()
		if line.startswith('>'):
			return 'fasta'
		elif line.startswith('@'):
			return 'fastq'
		else:
			raise RuntimeError('Input files must be in either FASTA or FASTQ format.')

def convert_to_fasta(input_list):
	fasta_dir = args.input + 'fastas/'
	if not os.path.exists(fasta_dir): 
		os.mkdir(fasta_dir)
	for f in input_list:
		out_file = os.path.join(fasta_dir, os.path.basename(f).split('.')[0])
		open(out_file, 'w').write('')
		out_handle = open(out_file, 'a')
		for s in SeqIO.parse(f, 'fastq'):
			out_handle.write('>{0}\n{1}\n'.format(s.id, str(s.seq)))
	return list_files(fasta_dir)

def merge_reads():
	import pandaseq
	merge_dir = args.input + 'merged_reads/'
	if not os.path.exists(merge_dir): 
		os.mkdir(merge_dir)
	pandaseq.run(args.input, merge_dir, nextseq=args.next_seq)
	return list_files(merge_dir)

def preprocess(files):
	import pre_processing
	processed_dir = args.input + 'processed/'
	if not os.path.exists(processed_dir): 
		os.mkdir(processed_dir)
	pre_processing.run(files, processed_dir)
	return list_files(processed_dir)

def download_files():
	from basespace import BaseSpace
	bs = BaseSpace()
	bs.download(args.input)
	args.merge = True



#####################################################################
# 
#                            IgBLAST
#
#####################################################################

def do_igblast(i, out_dir):
	o_prefix = os.path.basename(i).split('.')[0]
	o = os.path.join(out_dir, o_prefix + '_blastout')

	# parallel IgBLASTn
	blast_start = time.time()
	blastout = parallel_igblast(i,o)
	blast_end = time.time()
	blast_time = blast_end - blast_start

	# parse the IgBLASTn output
	parse_start = time.time()
	good_seqs, exc_seqs, failed_seqs = do_parse(blastout)
	parse_end = time.time()
	parse_time = parse_end - parse_start
	print_summary_output(good_seqs, exc_seqs, failed_seqs, blast_time, parse_time)


def parallel_igblast(in_file, out_file):	
	num_seqs = file_length(in_file)
	threads = num_procs()
	split_length = int(math.ceil(float(num_seqs) / threads))
	temp_directory, temp_out_directory = build_temp_dirs()
	split_files = file_splitter(in_file, split_length, num_seqs, temp_directory)
	thread_list = []
	blastout_list = []

	# run IgBLASTn in parallel
	for f in split_files:
		temp_out_file = os.path.join(temp_out_directory, os.path.basename(f).split('.')[0] + "_blastout")
		t = launch_thread(f, temp_out_file)
		t.start()
		thread_list.append(t)
		blastout_list.append(temp_out_file)	
	for thread in thread_list:
		thread.join()

	# combine all blastout files into a single output file
	open(out_file, 'w').write('')
	with open(out_file, 'w') as out_handle:
		for f in blastout_list:
			with open(f) as in_handle:
				for line in in_handle:
					out_handle.write(line)
	clean_up(split_files, blastout_list, temp_directory, temp_out_directory)
	return out_file


def main():
	# if args.use_basespace:
	# 	download_files()
	if args.merge:
		input_list = merge_reads()
	else:
		input_list = list_files(args.input)
		input_list = check_input(input_list)
	input_list = preprocess(input_list)
	print_input_info(input_list)
	output_dir = build_output_dir()
	for i in input_list:
		if os.path.isfile(i):
			print_infile(i)
			do_igblast(i, output_dir)


if __name__ == '__main__':
	main()
		
		
		
		
		