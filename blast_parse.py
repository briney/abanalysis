#!/usr/bin/python
# filename: blast_parse.py

###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 0.3.0
# @author: Bryan Briney
# @props: IgBLAST team (http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi)
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import re
import sys
import math
import json
import traceback
import collections
import warnings

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo as matlist


class BlastParse(object):


	def __init__(self, input, species='human', log='', tsv=False, debug=False, uaid=None):

		self.input = input.split('\n')
		self.species = species
		self.debug = debug
		self.tsv = tsv
		self.uaid_len = uaid
		self.log = self._get_log(log)


	def __str__(self):
		return '\n'.join(i for i in self.input)


	def get_id(self):
		return self.input[0]




######################################################
#
#                      PARSE
#
######################################################


	def parse(self):

		warnings.filterwarnings("ignore")
		try:
			# break the input into region-defined chunks
			self._get_chunks()

			# parse each chunk
			self._parse_chunks()

			# assemble the vdj sequences
			self._assemble_vdj()

			# regions
			self._var_regions()

			# junction
			self._find_junction()

			# identities
			self._get_identities()

			# output
			return self._build_output()

		except Exception, e:
			e = traceback.format_exc()
			return self._exception_output(e)




######################################################
#
#                  SANITY CHECKS
#
######################################################


	def sanity_checks(self):

		try:
			if self._rearrangement_alignment_check() > 0: return 1
			probs = self._gene_check()
			if probs > 0: return 1	
			if self._region_check() > 0: return 1
			return 0
		except Exception:
			return 1


	def _rearrangement_alignment_check(self):

		i = iter(self.input)
		row = next(i)
		try:
			while row.find('V-(D)-J rearrangement summary') == -1: row = next(i)
			while row.find('Alignment summary') == -1: row = next(i)
			return 0
		except StopIteration: 
			if self.debug: print "Rearrangement check failed."
			return 1


	def _gene_check(self):

		i = iter(self.input)
		row = next(i)

		# get to the rearrangement line, grab the v-gene
		while row.find('V-(D)-J rearrangement summary') == -1: row = next(i)
		line = next(i).split()
		v_gene = line[0]

		# heavy chains
		if v_gene[:3] == 'IGH':
			if line[3] != 'VH' or line[2][:2] != 'IG':
				if self.debug: print 'Gene check failed.'
				return 1

		# light chains
		elif v_gene[:3] in ('IGK', 'IGL'):
			if line[2] not in ('VK', 'VL') or line[1][:2] != 'IG':
				if self.debug: print 'Gene check failed.'
				return 1

		return  0


	def _region_check(self):

		regions = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT']
		i = iter(self.input)
		row = next(i)

		# get to the start of the regions block
		while row.find('Alignment summary') == -1: row = next(i)
		self.first_region = regions.index(next(i).split()[0])
		countdown = regions.index('FR3-IMGT') - self.first_region

		# check the FWR3 line
		if countdown == 0: return 0
		while countdown > 0:
			row = next(i)
			countdown -= 1
		if row.split()[0] != 'FR3-IMGT':
			if self.debug: print "Region check failed. Didn't find FR3-IMGT region. Found '{}' instead.".format(row.split()[0])
			return 1
		elif row.split()[3] == 'N/A' or int(row.split()[3]) < 25:
			if self.debug: print 'Region check failed.'
			return 1

		# check the Total line
		while row.find('Total') == -1: row = next(i)
		if row.split()[3] == 'N/A' or int(row.split()[3]) < 25:
			if self.debug: print 'Region check failed.'
			return 1

		return 0




######################################################
#
#                 SPLIT INTO CHUNKS
#
######################################################


	def _get_chunks(self):

		self.chunks = {}
		self.iter = iter(self.input)
		row = self._query_chunk()
		row = self._score_chunk(row)
		row = self._classification_chunk(row)
		row = self._rearrangement_chunk(row)
		row = self._junction_chunk(row)
		row = self._alignment_summary_check(row)
		self._alignment_chunk(row)


	def _query_chunk(self):

		row = next(self.iter)
		self.chunks['query'] = []
		while row != '':
			self.chunks['query'].append(row.strip())
			row = next(self.iter)
		return next(self.iter)


	def _score_chunk(self, row):

		self.chunks['score'] = []
		while row.find('Sequences producing significant alignments:') == -1: 
			row = next(self.iter)
		row = next(self.iter)
		row = next(self.iter)
		while row != '': 
			self.chunks['score'].append(row.rstrip('\n').split())
			row = next(self.iter)
		return row


	def _classification_chunk(self, row):

		while row.find('Domain classification') == -1: row = next(self.iter)
		row = next(self.iter)
		self.chunks['domain_class'] = row.rstrip('\n').split()
		return next(self.iter)


	def _rearrangement_chunk(self, row):

		while row.find('V-(D)-J rearrangement summary') == -1: row = next(self.iter)
		row = next(self.iter)
		self.chunks['rearrangement_summary'] = row.rstrip('\n').split()
		return next(self.iter)


	def _junction_chunk(self, row):

		while row.find('V-(D)-J junction details') == -1: row = next(self.iter)
		row = next(self.iter)
		self.chunks['junction_details'] = row.rstrip('\n').split()
		return next(self.iter)


	def _alignment_summary_check(self, row):

		self.chunks['alignment_summary'] = []
		while row.find('Alignment summary') == -1: row = next(self.iter)
		row = next(self.iter)
		while row != '':
			self.chunks['alignment_summary'].append(row)
			row = next(self.iter)
		return next(self.iter)


	def _alignment_chunk(self, row):

		self.chunks['alignment'] = []
		while row.find('Alignments') == -1: row = next(self.iter)
		while row.find('Lambda') == -1: 
			self.chunks['alignment'].append(row.rstrip('\n'))
			row = next(self.iter)
		while self.chunks['alignment'][0] == '': del self.chunks['alignment'][0]
		#while self.chunks['alignment'][-1] == '': del self.chunks['alignment'][-1]




######################################################
#
#                   PARSE CHUNKS
#
######################################################


	def _parse_chunks(self):

		self._parse_query_chunk(self.chunks['query'])
#		self._parse_classification_chunk(self.chunks['domain_class'])
		self._parse_rearrangement_chunk(self.chunks['rearrangement_summary'])
		self._parse_score_chunk(self.chunks['score'])
		self._parse_alignment_chunk(self.chunks['alignment'])


	def _parse_query_chunk(self, chunk):
		c_string = ''.join(chunk)
		c_string = c_string.replace('\n', '')
		self.seq_id = '_'.join(c_string.split('_')[:-2])
		self.raw_input = c_string.split('_')[-2]
		# self.uaid = ''
		# if self.uaid_len:
		# 	self.uaid = self.raw_input[:self.uaid_len]
		self.uaid = self.raw_input[:self.uaid_len] if self.uaid_len else ''
			# self.seq_id = self.seq_id.split('_')[0]
		if self.debug: print '\n\n' + self.seq_id


	def _parse_classification_chunk(self, chunk):

		pass
		#TODO


	def _parse_rearrangement_chunk(self, chunk):

		chains = {	'IGHV': 'heavy',
					'VH':   'heavy',
					'IGKV': 'kappa',
					'IGLV': 'lambda'}

		self.chain = chains[re.split('[1-9]', chunk[0])[0]]	

		self.var_gene = chunk[0]
		self.join_gene = chunk[-6]
		if self.chain == 'heavy':
			if 'IGHD' in chunk[1]: self.div_gene = chunk[1]
			else: self._hc_without_div()
		else: self.has_div = False
		self._gene_ids()

		self.productivity = chunk[-2]


	def _parse_score_chunk(self, chunk):

		self.var_bitscore = float(chunk[0][1])
		self.var_evalue = float(chunk[0][2])

		self.join_bitscore = float(chunk[-1][1])
		self.join_evalue = float(chunk[-1][2])

		if self.chain == 'heavy' and self.div_gene != '':
			self.div_bitscore = float(chunk[1][1])
			self.div_evalue = float(chunk[1][2])


	def _parse_alignment_chunk(self, chunk):

		self.blocks = self._build_alignment_blocks(chunk)
		self._get_pretty_alignment(self.blocks)

		self._v_alignment_parser(self.blocks)
		self._junction_sequence_parser(self.chunks['junction_details'])
		if self.has_div: self._d_alignment_parser(self.blocks)
		self._j_alignment_parser(self.blocks)




######################################################
#
#                    GENE IDs
#
######################################################


	def _gene_ids(self):

		self._v_gene_id()
		self._j_gene_id()
		if self.chain == 'heavy' and self.div_gene != '': self._d_gene_id()


	def _v_gene_id(self):

		v_split = re.split('[-\*]', self.var_gene)
		self.var_gene_fam = re.sub('VH|IGHV|VK|IGKV|VL|IGLV', '', v_split[0])
		if len(v_split) > 3: self.var_gene_gene = '-'.join(v_split[1:3])
		else: self.var_gene_gene = v_split[1]
		self.var_gene_allele = v_split[-1]


	def _j_gene_id(self):

		j_split = re.split('[-\*]', self.join_gene)
		self.join_gene_gene = re.sub('IGHJ|IGKJ|IGLJ', '', j_split[0])
		self.join_gene_allele = j_split[1]


	def _d_gene_id(self):

		d_split = re.split('[-\*]', self.div_gene)
		self.div_gene_fam = re.sub('IGHD|IGKD|IGLD', '', d_split[0])
		self.div_gene_gene = d_split[1]
		self.div_gene_allele = d_split[2]


	def _hc_without_div(self):

		self.div_gene = ''
		self.has_div = False
		self.div_gene_fam = ''
		self.div_gene_gene = ''
		self.div_gene_allele = ''
		self.div_bitscore = ''
		self.div_evalue = ''
		self.div_nt_identity = ''
		self.div_muts_nt = ''
		self.div_muts_nt_count = ''
		self.div_ins = ''
		self.div_del = ''




######################################################
#
#                ALIGNMENT PARSING
#
######################################################


	def _build_alignment_blocks(self, chunk):

		# define the start point for each alignment block
		starts = []
		for i, line in enumerate(chunk):
			if 'Query_' in line: starts.append(i-2)

		# build a list of the alignment blocks
		blocks = []
		for s in starts:
			block = []
			i = iter(chunk[s:])
			row = next(i)
			try:
				while row.split() == '': row = next(i)
				while row != '':
					block.append(row)
					row = next(i)
				blocks.append(block)
			except StopIteration: pass
		return blocks


	def _get_pretty_alignment(self, blocks):

		self.pretty_alignment = '\n\n'.join(['\n'.join(b) for b in blocks])


	def _get_vdj_alignment_blocks(self, blocks, gene):

		a_blocks = []
		for block in blocks:
			for line in block:
				if line[0] == gene: 
					a_blocks.append(block)
		return a_blocks


	def _align_block_position(self, block):

		position = re.search('<', block[0]).start()
		length = len(block[0].strip())

		return position, length



#-----------------
#   VARIABLE
#-----------------


	def _v_alignment_parser(self, blocks):

		# get the raw data out of the alignment blocks
		v_blocks = self._get_vdj_alignment_blocks(blocks, 'V')
		self.v_region_delims = self._var_region_delimiters(v_blocks)
		position, length = self._align_block_position(v_blocks[0])
		self._var_germ_aa_sequence(v_blocks, position, length)
		self._var_aa_sequence(v_blocks, position, length)
		self._var_nt_alignment(v_blocks)
		self._var_nt_sequence(v_blocks)
		self.var_readframe = self._var_reading_frame()
		self._var_germ_offsets(v_blocks)

		# fix ambigs and indels
		if 'N' in self.var_nt_seq: 
			self._fix_v_ambigs()
		self.var_ins = []
		self.var_del = []
		if self._v_indel_check(): 
			self._fix_v_indels()
		# self.var_aa_seq = self._v_retranslate()

		# mutations
		self._var_nt_mutations()
		self._var_aa_mutations()


	def _var_region_delimiters(self, blocks):

		delims = ''
		for block in blocks:
			delims += block[0].strip()
		return delims


	def _var_germ_aa_sequence(self, blocks, p, l):

		self.var_germ_aa_seq_raw = ''
		for block in blocks:
			self.var_germ_aa_seq_raw += block[4][p:p+l]
		self.var_germ_aa_seq = ''.join(self.var_germ_aa_seq_raw.split())


	def _var_aa_sequence(self, blocks, p, l):

		self.var_aa_seq_raw = ''
		for block in blocks:
			self.var_aa_seq_raw += block[1][p:p+l]
		self.var_aa_seq = ''.join(self.var_aa_seq_raw.split())[:len(self.var_germ_aa_seq)]


	def _var_nt_alignment(self, blocks):

		align = ''
		for block in blocks:
			align += block[3].split()[5]
		while align[-1] == '-': align = align[:-1]
		self.var_nt_alignment = align


	def _var_nt_sequence(self, blocks):

		nt_seq = ''
		for block in blocks:
			nt_seq += block[2].split()[2]
		self.var_nt_seq = nt_seq[:len(self.var_nt_alignment)]


	def _var_reading_frame(self):

		alignment_start = re.search("\S", self.var_aa_seq_raw).start()
		return alignment_start - 1


	def _var_germ_offsets(self, blocks):

		self.var_germ_offset_nt = int(blocks[0][3].split()[4])
		self.var_germ_offset_aa = int((self.var_germ_offset_nt - 1 + self.var_readframe) / 3)


	def _var_nt_mutations(self):

		self.var_muts_nt = []
		for i in range(len(self.var_nt_seq)):
			if self.var_nt_alignment[i] != '.':
				m = i + self.var_germ_offset_nt
				mut = '{0}:{1}>{2}'.format(m, self.var_nt_alignment[i], self.var_nt_seq[i])
				self.var_muts_nt.append(mut)

		self.var_muts_nt_count = len(self.var_muts_nt)


	def _var_aa_mutations(self):

		self.var_muts_aa = []
		[(a1, a2, score, begin, end)] = pairwise2.align.globalxx(self.var_germ_aa_seq, self.var_aa_seq, one_alignment_only=True)
		for i in range(len(a1)):
			if a1[i] == '-' or a2[i] == '-':
				continue
			if a1[i] != a2[i]:
				m = i + self.var_germ_offset_aa
				mut = '{0}:{1}>{2}'.format(m, a1[i], a2[i])
				self.var_muts_aa.append(mut)

		self.var_muts_aa_count = len(self.var_muts_aa)



#-----------------
#   DIVERSITY
#-----------------


	def _d_alignment_parser(self, blocks):

		# get the raw data out of the alignment blocks
		d_blocks = self._get_vdj_alignment_blocks(blocks, 'D')
		d_start, d_length = self._div_nt_alignment(d_blocks)
		self._div_nt_seq(d_blocks, d_start, d_length)

		# fix ambigs and indels
		if 'N' in self.div_nt_seq: self._fix_d_ambigs
		self.div_ins = []
		self.div_del = []
		if self._d_indel_check():
			self._fix_d_indels()

		# mutations
		self._div_nt_mutations()


	def _div_nt_alignment(self, blocks):

		self.div_nt_alignment = ''
		for block in blocks:
			i = iter(block)
			row = next(i)
			while row[0] != 'D': row = next(i)
			self.div_nt_alignment += row.split()[5]
		return self._trim_div_nt_alignment()


	def _trim_div_nt_alignment(self):
		
		untrimmed_length = len(self.div_nt_alignment)
		while self.div_nt_alignment[0] == '-': self.div_nt_alignment = self.div_nt_alignment[1:]
		start = untrimmed_length - len(self.div_nt_alignment) + self.d_trunc_start
		while self.div_nt_alignment[-1] == '-': self.div_nt_alignment = self.div_nt_alignment[:-1]
		if self.d_trunc_end > 0:
			self.div_nt_alignment[:-self.d_trunc_end]
		return start, len(self.div_nt_alignment)


	def _div_nt_seq(self, blocks, s, l):

		self.div_nt_seq = ''
		for block in blocks:
			i = iter(block)
			row = next(i)
			while 'Query_' not in row: row = next(i)
			self.div_nt_seq += row.split()[2]
		self.div_nt_seq = self.div_nt_seq[s:s+l]


	def _div_nt_mutations(self):

		self.div_muts_nt = []
		for i in range(len(self.div_nt_seq)):
			if self.div_nt_alignment[i] != '.':
				m = i + self.var_germ_offset_nt + (len(self.var_nt_seq + self.n1_nt))
				mut = '{0}:{1}>{2}'.format(m, self.div_nt_alignment[i], self.div_nt_seq[i])
				self.div_muts_nt.append(mut)

		self.div_muts_nt_count = len(self.div_muts_nt)



#-----------------
#   JOINING
#-----------------


	def _j_alignment_parser(self, blocks):

		# get the raw data out of the alignment blocks
		j_blocks = self._get_vdj_alignment_blocks(blocks, 'J')
		j_length = self._join_nt_alignment(j_blocks)
		self._join_nt_seq(j_blocks, j_length)

		# fix ambigs and indels
		if 'N' in self.join_nt_seq: self._fix_j_ambigs()
		self.join_ins = []
		self.join_del = []
		if self._j_indel_check():
			self._fix_j_indels()

		# mutations
		self._join_nt_mutations()


	def _join_nt_alignment(self, blocks):

		self.join_nt_alignment = ''
		for block in blocks:
			i = iter(block)
			row = next(i)
			while row[0] != 'J': row = next(i)
			try:
				self.join_nt_alignment += row.split()[5]
			except IndexError:
				if row.split()[4].endswith('-.'): self.join_nt_alignment += row.split()[4]
				else: raise Exception
		while self.join_nt_alignment[0] == '-': self.join_nt_alignment = self.join_nt_alignment[1:]
		return len(self.join_nt_alignment)


	def _join_nt_seq(self, blocks, l):

		self.join_nt_seq = ''
		for block in blocks:
			i = iter(block)
			row = next(i)
			while 'Query_' not in row: row = next(i)
			self.join_nt_seq += row.split()[2]
		start = len(self.join_nt_seq) - l
		self.join_nt_seq = self.join_nt_seq[start:]


	def _join_nt_mutations(self):

		self.join_muts_nt = []
		for i in range(len(self.join_nt_seq)):
			if self.join_nt_alignment[i] != '.':
				m = i + self.join_start
				mut = '{0}:{1}>{2}'.format(m, self.join_nt_alignment[i], self.join_nt_seq[i])
				self.join_muts_nt.append(mut)

		self.join_muts_nt_count = len(self.join_muts_nt)




######################################################
#
#                CDR3 AND JUNCTION
#
######################################################


	def _junction_sequence_parser(self, j):

		for i in range(len(j)):
			if j[i] == 'N/A': j[i] = ''

		if self.chain == 'heavy':
			self.junction_var = j[0]
			self.junction_join = j[4]
			self._n1_parser(j[1])
			self._n2_parser(j[3])
			self._div_nt_parser(j[2])
		else:
			self.junction_var = j[0]
			self.junction_join = j[2]
			self._n_parser(j[1])
		self._get_join_start_position()


	def _n1_parser(self, seq):

		if seq.startswith('('):
			self.d_trunc_start = len(seq) - 2
			self.n1_nt = ''
			self.junction_var = self.junction_var + re.sub('[()]', '', seq)
		else:
			self.d_trunc_start = 0
			self.n1_nt = seq


	def _n2_parser(self, seq):

		if seq.startswith('('):
			self.d_trunc_end = len(seq) - 2
			self.n2_nt = ''
			self.junction_join = re.sub('[()]', '', seq) + self.junction_join
		else:
			self.d_trunc_end = 0
			self.n2_nt = seq


	def _div_nt_parser(self, seq):

		if seq == '':
			self.has_div = False
			self.div_nt = ''
		elif seq.startswith('('):
			self.has_div = False
			self.div_nt = ''
			if self.n1_nt == '':
				self.junction_var = self.junction_var + re.sub('[()]', '', seq)
			else:
				self.junction_join = re.sub('[()]', '', seq) + self.junction_join
		else:
			self.has_div = True
			self.div_nt = seq

		self.junction_nt = self.junction_var + self.n1_nt + self.div_nt + self.n2_nt + self.junction_join

	
	def _n_parser(self, seq):

		if seq.startswith('('):
			self.n_nt = ''
			self.junction_var = self.junction_var + re.sub('[()]', '', seq)
			self._j_overlap_offset = len(re.sub('[()]', '', seq))
		else:
			self.n_nt = seq
			self._j_overlap_offset = 0

		self.junction_nt = self.junction_var + self.n_nt + self.junction_join


	def _get_join_start_position(self):

		if self.chain == 'heavy':
			self.join_start = self.var_germ_offset_nt + len(self.var_nt_seq + self.n1_nt + self.div_nt + self.n2_nt)
		else:
			self.join_start = self.var_germ_offset_nt + len(self.var_nt_seq + self.n_nt) + self._j_overlap_offset


	def _find_junction(self):

		start = self._find_junction_start()
		end = self._find_junction_end(start)

		self.junction_aa = self.vdj_aa[start:end]
		self.CDR3_aa = self.vdj_aa[start+1:end-1]


	def _find_junction_start(self):

		scores = []
		matrix = matlist.blosum62

		fr3_start = len(self.FR1_aa + self.CDR1_aa + self.FR2_aa + self.CDR2_aa)
		germ = self._get_v_chunk(self.var_gene)
		for i in range(fr3_start, len(self.vdj_aa)-5):
			chunk = self.vdj_aa[i:i+6]
			if '*' in chunk:
				scores.append(0.0)
				continue
			scores.append(pairwise2.align.globalds(germ, chunk, matrix, -100, -20, one_alignment_only=1, score_only=1))

		max_score = max(scores)
		max_score_positions = [i+4 for i,j in enumerate(scores) if j == max_score]
		return (max_score_positions[0] + fr3_start)


	def _find_junction_end(self, start):

		scores = []
		matrix = matlist.blosum62

		germ = self._get_j_chunk(self.join_gene)
		for i in range(start, len(self.vdj_aa)-1):
			if i < len(self.vdj_aa)-4:
				chunk = self.vdj_aa[i:i+5]
			else:
				chunk = self.vdj_aa[i:]
			if '*' in chunk:
				scores.append(0.0)
				continue
			scores.append(pairwise2.align.globalds(germ[:len(chunk)], chunk, matrix, -100, -20, one_alignment_only=1, score_only=1))

		max_score = max(scores)
		max_score_positions = [i+1 for i,j in enumerate(scores) if j == max_score]
		return (max_score_positions[-1] + start)


	def _get_v_chunk(self, gene):

		trunc_gene = gene[:5]
		germs = {
			'human':	{	'IGHV1' : 'AVYYCA',
							'IGHV2' : 'ATYYCA',
							'IGHV3' : 'AVYYCA',
							'IGHV4' : 'AVYYCA',
							'IGHV5' : 'AMYYCA',
							'IGHV6' : 'AVYYCA',
							'IGHV7' : 'AVYYCA',
							'IGKV1' : 'ATYYCQ',
							'IGKV2' : 'GVYYCM',
							'IGKV3' : 'AVYYCQ',
							'IGKV4' : 'AVYYCQ',
							'IGLV1' : 'ADYYCQ',
							'IGLV2' : 'ADYYCS',
							'IGLV3' : 'ADYYCQ',
							'IGLV4' : 'ADYYCQ',
							'IGLV5' : 'ADYYCM',
							'IGLV6' : 'ADYYCQ',
							'IGLV7' : 'AEYYCL',
							'IGLV8' : 'SDYYCV',
							'IGLV9' : 'SDYHCG'},

			'macaque':	{	'IGHV1' : 'AVYYCA',
							'IGHV2' : 'ATYYCA',
							'IGHV3' : 'AVYYCA',
							'IGHV4' : 'AVYYCA',
							'IGHV5' : 'ATYYCA',
							'IGHV6' : 'AVYYCA',
							'IGHV7' : 'AVYYCA',
							'IGKV1' : 'ATYYCQ',
							'IGKV2' : 'GVYYCM',
							'IGKV3' : 'AVYYCQ',
							'IGKV4' : 'AVYYCQ',
							'IGKV5' : 'AYYFCQ',
							'IGKV6' : 'ATYYCQ',
							'IGKV7' : 'ADYYCL',
							'IGLV1' : 'ADYYCQ',
							'IGLV2' : 'ADYYCS',
							'IGLV3' : 'ADYYCQ',
							'IGLV4' : 'ADYYCQ',
							'IGLV5' : 'ADYYCM',
							'IGLV6' : 'ADYYCQ',
							'IGLV7' : 'AEYYCW',
							'IGLV8' : 'SDYYCT',
							'IGLV9' : 'SDYHCG'
							},
			
			'mouse':	{	'IGHV1' : 'AVYYCA',
							'IGHV2' : 'AIYYCA',
							'IGHV3' : 'ATYYCA',
							'IGHV4' : 'ALYYCA',
							'IGHV5' : 'AMYYCA',
							'IGHV6' : 'GIYYCT',
							'IGHV7' : 'ATYYCA',
							'IGHV8' : 'ATYYCA',
							'IGHV9' : 'ATYFCA',
							'IGKV1' : 'ATYYCQ',
							'IGKV2' : 'GVYYCA',
							'IGKV3' : 'ATYYCQ',
							'IGKV4' : 'ATYYCQ',
							'IGKV5' : 'GVYYCQ',
							'IGKV6' : 'AVYFCQ',
							'IGKV7' : 'THYYCA',
							'IGKV8' : 'AVYYCQ',
							'IGKV9' : 'ADYYCL',
							'IGLV1' : 'AIYFCA',
							'IGLV2' : 'AMYFCA',
							'IGLV3' : 'AIYICG',
							'IGLV4' : 'AIYFCA',
							'IGLV5' : 'AIYFCA',
							'IGLV6' : 'AIYFCA',
							'IGLV7' : 'AIYFCA',
							'IGLV8' : 'AIYFCA'},
			
			'rabbit':	{}}

		return germs[self.species][trunc_gene]


	def _get_j_chunk(self, gene):

		trunc_gene = gene[:5]
		if self.debug:
			print self.species
			print trunc_gene
		germs = {
			'human':	{	'IGHJ1' : 'WGQGT', 
							'IGHJ2' : 'WGRGT', 
							'IGHJ3' : 'WGQGT', 
							'IGHJ4' : 'WGQGT',
							'IGHJ5' : 'WGQGT',
							'IGHJ6' : 'WGQGT',
							'IGKJ1' : 'FGQGT',
							'IGKJ2' : 'FGQGT',
							'IGKJ3' : 'FGPGT',
							'IGKJ4' : 'FGGGT',
							'IGKJ5' : 'FGQGT',
							'IGLJ1' : 'FGTGT',
							'IGLJ2' : 'FGGGT',
							'IGLJ3' : 'FGGGT',
							'IGLJ4' : 'FGGGT',
							'IGLJ5' : 'FGEGT',
							'IGLJ6' : 'FGSGT',
							'IGLJ7' : 'FGGGT'},

			'macaque':	{	'IGHJ1' : 'WGQGA',
							'IGHJ2' : 'WGPGT',
							'IGHJ3' : 'WGQGL',
							'IGHJ4' : 'WGQGV',
							'IGHJ5' : 'WGPGV',
							'IGHJ6' : 'WGQGV',
							'IGKJ1' : 'FGQGT',
							'IGKJ2' : 'FGQGT',
							'IGKJ3' : 'FGPGT',
							'IGKJ4' : 'FGGGT',
							'IGKJ5' : 'FGQGT',
							'IGLJ1' : 'FGAGT',
							'IGLJ2' : 'FGGGT',
							'IGLJ3' : 'FGGGT',
							'IGLJ4' : 'FCGGT',
							'IGLJ5' : 'FGEGT',
							'IGLJ6' : 'FGSGT'},

			'mouse':	{	'IGHJ1' : 'WGAGT',
							'IGHJ2' : 'WGQGT',
							'IGHJ3' : 'WGQGT',
							'IGHJ4' : 'WGQGT',
							'IGKJ1' : 'FGGGT',
							'IGKJ2' : 'FGGGT',
							'IGKJ4' : 'FGSGT',
							'IGKJ5' : 'FGAGT',
							'IGLJ1' : 'FGGGT',
							'IGLJ2' : 'FGGGT',
							'IGLJ3' : 'FGSGT'},

			'rabbit':	{}}

		return germs[self.species][trunc_gene]




######################################################
#
#                AMBIGS AND INDELS
#
######################################################


#-----------------
#   VARIABLE
#-----------------


	def _fix_v_ambigs(self):
		for n in re.finditer('N|n', self.var_nt_seq):
			i = n.start()
			if self.var_nt_alignment[i] != '-':
				self.var_nt_seq = self.var_nt_seq[:i] + self.var_nt_alignment[i] + self.var_nt_seq[i+1:]
				self.var_nt_alignment = self.var_nt_alignment[:i] + '.' + self.var_nt_alignment[i+1:]


	def _v_indel_check(self):
		i = iter(self.chunks['alignment_summary'])
		row = next(i)
		while row.find('Total') == -1: 
			row = next(i)
		if int(row.split()[6]) > 0: return True
		return False


	def _fix_v_indels(self):
		self.var_ins = []
		self.var_del = []
		self._fix_combined_nfs_indels()
		self._v_retranslate()

	
	def _fix_combined_nfs_indels(self):
		'''Accounts for a quirk in IgBLAST (likely due to incorrect gap open/extend penalties) which
		   can result in a single non-frameshift indel being split into multiple samller, frameshift indels.
		'''
		full_nt_germline = self._full_v_nt_germline()
		nt_seq = self.var_nt_seq.replace('-', '')
		self._recalc_combined_indel_germline_alignment(nt_seq, full_nt_germline)
		ins = [i for i in re.finditer('-+', self.var_nt_alignment)]
		if len(ins) > 0:
			self._fix_v_ins(ins)
		dels = [d for d in re.finditer('-+', self.var_nt_seq)]
		if len(dels) > 0:
			self._fix_v_dels(dels)


	def _full_v_nt_germline(self):
		full_nt_germline = ''
		for i, nt in enumerate(self.var_nt_alignment):
			if nt == '.':
				full_nt_germline += self.var_nt_seq[i]
			else:
				full_nt_germline += self.var_nt_alignment[i]
		return full_nt_germline.replace('-', '')


	def _recalc_combined_indel_germline_alignment(self, nt_seq, full_nt_alignment):
		[alignment] = pairwise2.align.globalms(nt_seq, full_nt_alignment, 1, 0, -5, -1, one_alignment_only=True)
		aligned_seq = alignment[0]
		aligned_germ = alignment[1]
		dot_alignment = ''
		for i, g in enumerate(aligned_germ):
			if g == aligned_seq[i]:
				dot_alignment += '.'
			else:
				dot_alignment += g
		self.var_nt_alignment = dot_alignment

		self.var_nt_seq = aligned_seq

	
	def _fix_v_ins(self, ins):
		ins.sort(key=lambda x: x.start(), reverse=True)
		for i in ins:
			l = i.end() - i.start()
			if self.debug:
				print i.start()
				print i.end()
				print l
			if l % 3 == 0: 
				self._v_nfs_ins(i.start(), i.end())
			else: 
				self._v_fs_ins(i.start(), i.end())


	def _v_nfs_ins(self, s, e):
		self.var_ins.append('{0},{1},{2}'.format(s + self.var_germ_offset_nt, e-s, self.var_nt_seq[s:e]))
		self.var_nt_alignment = self.var_nt_alignment[:s] + '.'*(e-s) + self.var_nt_alignment[e:]
		aa_len = (e-s)/3
		self.var_germ_aa_seq_raw = self.var_germ_aa_seq_raw[:s] + '-'*aa_len + self.var_germ_aa_seq_raw[s+aa_len:]
		self.var_germ_aa_seq = self.var_germ_aa_seq_raw.replace(' ', '')


	def _v_fs_ins(self, s, e):
		self.var_nt_seq = self.var_nt_seq[:s] + self.var_nt_seq[e:]
		self.var_nt_alignment = self.var_nt_alignment[:s] + self.var_nt_alignment[e:]
		self.v_region_delims = self.v_region_delims[:s] + self.v_region_delims[e:]


	def _fix_v_dels(self, dels):
		dels.sort(key=lambda x: x.start(), reverse=True)
		for d in dels:
			l = d.end() - d.start()
			if self.debug:
				print d.start()
				print d.end()
				print l
			if l % 3 == 0: 
				self._v_nfs_del(d.start(), d.end())
			else: 
				self._v_fs_del(d.start(), d.end())


	def _v_nfs_del(self, s, e):
		self.var_del.append('{0},{1},{2}'.format(s + self.var_germ_offset_nt, e-s, self.var_nt_alignment[s:e]))
		self.var_nt_alignment = self.var_nt_alignment[:s] + self.var_nt_alignment[e:]
		self.var_nt_seq = self.var_nt_seq[:s] + self.var_nt_seq[e:]
		self.var_germ_aa_seq_raw = self.var_germ_aa_seq_raw[:s] + self.var_germ_aa_seq_raw[e:]
		self.var_germ_aa_seq = self.var_germ_aa_seq_raw.replace(' ', '')


	def _v_fs_del(self, s, e):
		self.var_nt_seq = self.var_nt_seq[:s] + self.var_nt_alignment[s:e] + self.var_nt_seq[e:]
		self.var_nt_alignment = self.var_nt_alignment[:s] + '.'*(e-s) + self.var_nt_alignment[e:]


	def _v_retranslate(self):
		self.var_aa_seq = Seq(self.var_nt_seq[self.var_readframe:len(self.var_nt_alignment)], IUPAC.ambiguous_dna).translate()
		if self.debug:
			print 'var_aa:', self.var_aa_seq



#-----------------
#   DIVERSITY
#-----------------


	def _fix_d_ambigs(self):

		for n in re.finditer('N|n', self.div_nt_seq):
			i = n.start()
			if self.div_nt_alignment[i] != '-':
				self.div_nt_seq = self.div_nt_seq[:i] + self.div_nt_alignment[i] + self.div_nt_seq[i+1:]
				self.div_nt_alignment = self.div_nt_alignment[:i] + '.' + self.div_nt_alignment[i+1:]


	def _d_indel_check(self):

		if '-' in self.div_nt_alignment: return True
		if '-' in self.div_nt_seq: return True
		return False


	def _fix_d_indels(self):

		ins = re.findall('-+', self.div_nt_alignment)
		dels = re.findall('-+', self.div_nt_seq)
		self.div_ins = []
		self.div_del = []

		if len(ins) > 0: self._fix_d_ins()
		if len(dels) > 0: self._d_fix_dels() 


	def _fix_d_ins(self):

		o = 0
		for i in re.finditer('-+', self.div_nt_alignment):
			s = i.start() - o
			e = i.end() - o
			l = e - s
			if l % 3 == 0: 
				self._d_nfs_ins(s, e)
			else: 
				self._d_fs_ins(s, e)
				o += l


	def _d_nfs_ins(self, s, e):

		self.div_ins.append('{0},{1},{2}'.format(s + self.div_germ_offset_nt, e-s, self.div_nt_seq[s:e]))
		self.div_nt_alignment = self.div_nt_alignment[:s] + '.'*(e-s) + self.div_nt_alignment[e:]
		if self.div_nt != self.div_nt_seq:
			self.join_start = self.var_germ_offset_nt + len(self.var_nt_seq + self.n1_nt + self.div_nt_seq + self.n2_nt)


	def _d_fs_ins(self, e, s):

		self.div_nt_seq = self.div_nt_seq[:s] + self.div_nt_seq[e:]
		self.div_nt_alignment = self.div_nt_alignment[:s] + self.div_nt_alignment[e:]


	def _fix_d_dels(self):

		o = 0
		for i in re.finditer('-+', self.div_nt_seq):
			s = i.start() - o
			e = i.end() - o
			l = e - s
			if l % 3 == 0: 
				self._d_nfs_del(s, e)
			else: 
				self._d_fs_del(s, e)
				o += l


	def _d_nfs_del(self, s, e):

		self.div_del.append('{0},{1},{2}'.format(s + self.div_start, e-s, self.div_nt_alignment[s:e]))
		self.div_nt_alignment = self.div_nt_alignment[:s] + self.div_nt_alignment[e:]
		self.div_nt_seq = self.div_nt_seq[:s] + self.div_nt_seq[e:]
		if self.div_nt != self.div_nt_seq:
			self.join_start = self.var_germ_offset_nt + len(self.var_nt_seq + self.n1_nt + self.div_nt_seq + self.n2_nt)


	def _d_fs_del(self, s, e):

		self.div_nt_seq = self.div_nt_seq[:s] + self.div_nt_alignment[s:e] + self.div_nt_seq[e:]
		self.div_nt_alignment = self.div_nt_alignment[:s] + '.'*(e-s) + self.div_nt_alignment[e:]



#-----------------
#   JOINING
#-----------------


	def _fix_j_ambigs(self):

		for n in re.finditer('N|n', self.join_nt_seq):
			i = n.start()
			if self.join_nt_alignment[i] != '-':
				self.join_nt_seq = self.join_nt_seq[:i] + self.join_nt_alignment[i] + self.join_nt_seq[i+1:]
				self.join_nt_alignment = self.join_nt_alignment[:i] + '.' + self.join_nt_alignment[i+1:]


	def _j_indel_check(self):

		if '-' in self.join_nt_alignment: return True
		if '-' in self.join_nt_seq: return True
		return False


	def _fix_j_indels(self):

		ins = len(re.findall('-+', self.join_nt_alignment))
		dels = len(re.findall('-+', self.join_nt_seq))
		self.join_ins = []
		self.join_del = []

		if ins > 0: self._fix_j_ins()
		if dels > 0: self._fix_j_dels()


	def _fix_j_ins(self):

		o = 0
		for i in re.finditer('-+', self.join_nt_alignment):
			s = i.start() - o
			e = i.end() - o
			l = e - s
			if l % 3 == 0: 
				self._j_nfs_ins(s, e)
			else: 
				self._j_fs_ins(s, e)
				o += l


	def _j_nfs_ins(self, s, e):

		self.join_ins.append('{0},{1},{2}'.format(s + self.join_start, e-s, self.join_nt_seq[s:e]))
		self.join_nt_alignment = self.join_nt_alignment[:s] + '.'*(e-s) + self.join_nt_alignment[e:]


	def _j_fs_ins(self, e, s):

		self.join_nt_seq = self.join_nt_seq[:s] + self.join_nt_seq[e:]
		self.join_nt_alignment = self.join_nt_alignment[:s] + self.join_nt_alignment[e:]


	def _fix_j_dels(self):

		o = 0
		for i in re.finditer('-+', self.join_nt_seq):
			s = i.start() - o
			e = i.end() - o
			l = e - s
			if l % 3 == 0: 
				self._j_nfs_del(s, e)
			else: 
				self._j_fs_del(s, e)
				o += l


	def _j_nfs_del(self, s, e):

		self.join_del_append('{0},{1},{2}'.format(s + self.join_start, e-s, self.join_nt_alignment[s:e]))
		self.join_nt_alignment = self.join_nt_alignment[:s] + self.join_nt_alignment[e:]
		self.join_nt_seq = self.join_nt_seq[:s] + self.join_nt_seq[e:]


	def _j_fs_del(self, s, e):

		self.join_nt_seq = self.join_nt_seq[:s] + self.join_nt_alignment[s:e] + self.join_nt_seq[e:]
		self.join_nt_alignment = self.join_nt_alignment[:s] + '.'*(e-s) + self.join_nt_alignment[e:]




######################################################
#
#                  VDJ ASSEMBLY
#
######################################################


	def _assemble_vdj(self):

		if self.chain == 'heavy':
			self.vdj_nt = self.var_nt_seq + self.n1_nt + self.div_nt + self.n2_nt + self.join_nt_seq
		else:
			self.join_nt_seq = self.join_nt_seq[self._j_overlap_offset:] # prevents double-counting of residues that IgBLAST identifies as 'overlapping' in the V-J junction
			self.vdj_nt = self.var_nt_seq + self.n_nt + self.join_nt_seq

		self.vdj_aa = str(Seq(self.vdj_nt[self.var_readframe:], IUPAC.ambiguous_dna).translate())




######################################################
#
#                     REGIONS
#
######################################################


	def _var_regions(self):

		regions = ['FWR1', 'CDR1', 'FWR2', 'CDR2', 'FWR3']
		r = {}
		for region in regions:
			r[region] = {}
			r[region]['start'] = 0
			r[region]['end'] = 0

		d = [len(i)+1 for i in self.v_region_delims.split('>')]
		first_region = 5-self.first_region

		if first_region == 5:
			r['FWR1']['end'] = d[first_region - 5]
		if first_region >= 4:
			r['CDR1']['start'] = r['FWR1']['end']
			r['CDR1']['end'] = r['CDR1']['start'] + d[first_region - 4]
		if first_region >= 3:
			r['FWR2']['start'] = r['CDR1']['end']
			r['FWR2']['end'] = r['FWR2']['start'] + d[first_region - 3]
		if first_region >= 2:
			r['CDR2']['start'] = r['FWR2']['end']
			r['CDR2']['end'] = r['CDR2']['start'] + d[first_region - 2]
		if first_region >= 1:
			r['FWR3']['start'] = r['CDR2']['end']
			r['FWR3']['end'] = r['FWR3']['start'] + d[first_region - 1]

		self._var_region_nt_seqs(r)
		self._var_region_aa_seqs()
		self._var_region_nt_muts(r)
		self._var_region_aa_muts(r)
		self._var_region_lengths()


	def _var_region_nt_seqs(self, r):

		self.FR1_nt = self.vdj_nt[r['FWR1']['start'] : r['FWR1']['end']]
		self.CDR1_nt = self.vdj_nt[r['CDR1']['start'] : r['CDR1']['end']]
		self.FR2_nt = self.vdj_nt[r['FWR2']['start'] : r['FWR2']['end']]
		self.CDR2_nt = self.vdj_nt[r['CDR2']['start'] : r['CDR2']['end']]
		self.FR3_nt = self.vdj_nt[r['FWR3']['start'] : r['FWR3']['end']]


	def _var_region_aa_seqs(self):

		self.FR1_aa = str(Seq(self.FR1_nt[self.var_readframe:], IUPAC.ambiguous_dna).translate())
		self.CDR1_aa = str(Seq(self.CDR1_nt, IUPAC.ambiguous_dna).translate())
		self.FR2_aa = str(Seq(self.FR2_nt, IUPAC.ambiguous_dna).translate())
		self.CDR2_aa = str(Seq(self.CDR2_nt, IUPAC.ambiguous_dna).translate())
		self.FR3_aa = str(Seq(self.FR3_nt, IUPAC.ambiguous_dna).translate())


	def _var_region_nt_muts(self, r):

		self.FR1_nt_muts = []
		self.FR2_nt_muts = []
		self.FR3_nt_muts = []
		self.CDR1_nt_muts = []
		self.CDR2_nt_muts = []

		o = self.var_germ_offset_nt
		for m in self.var_muts_nt:
			p = int(m.split(':')[0])
			if p < r['FWR1']['end'] + o: self.FR1_nt_muts.append(m)
			elif p < r['CDR1']['end'] + o: self.CDR1_nt_muts.append(m)
			elif p < r['FWR2']['end'] + o: self.FR2_nt_muts.append(m)
			elif p < r['CDR2']['end'] + o: self.CDR2_nt_muts.append(m)
			elif p < r['FWR3']['end'] + o: self.FR3_nt_muts.append(m)

		self.FR1_nt_muts_count = str(len(self.FR1_nt_muts))
		self.CDR1_nt_muts_count = str(len(self.CDR1_nt_muts))
		self.FR2_nt_muts_count = str(len(self.FR2_nt_muts))
		self.CDR2_nt_muts_count = str(len(self.CDR2_nt_muts))
		self.FR3_nt_muts_count = str(len(self.FR3_nt_muts))


	def _var_region_aa_muts(self, r):

		self.FR1_aa_muts = []
		self.FR2_aa_muts = []
		self.FR3_aa_muts = []
		self.CDR1_aa_muts = []
		self.CDR2_aa_muts = []

		o = self.var_germ_offset_nt
		for m in self.var_muts_aa:
			p = int(m.split(':')[0])
			if p < (r['FWR1']['end'] + o)/3: self.FR1_aa_muts.append(m)
			elif p < (r['CDR1']['end'] + o)/3: self.CDR1_aa_muts.append(m)
			elif p < (r['FWR2']['end'] + o)/3: self.FR2_aa_muts.append(m)
			elif p < (r['CDR2']['end'] + o)/3: self.CDR2_aa_muts.append(m)
			elif p < (r['FWR3']['end'] + o)/3: self.FR3_aa_muts.append(m)

		self.FR1_aa_muts_count = str(len(self.FR1_aa_muts))
		self.CDR1_aa_muts_count = str(len(self.CDR1_aa_muts))
		self.FR2_aa_muts_count = str(len(self.FR2_aa_muts))
		self.CDR2_aa_muts_count = str(len(self.CDR2_aa_muts))
		self.FR3_aa_muts_count = str(len(self.FR3_aa_muts))


	def _var_region_lengths(self):

		self.FR1_nt_length = str(len(self.FR1_nt))
		self.CDR1_nt_length = str(len(self.CDR1_nt))
		self.FR2_nt_length = str(len(self.FR2_nt))
		self.CDR2_nt_length = str(len(self.CDR2_nt))
		self.FR3_nt_length = str(len(self.FR3_nt))

		self.FR1_aa_length = str(len(self.FR1_aa))
		self.CDR1_aa_length = str(len(self.CDR1_aa))
		self.FR2_aa_length = str(len(self.FR2_aa))
		self.CDR2_aa_length = str(len(self.CDR2_aa))
		self.FR3_aa_length = str(len(self.FR3_aa))




######################################################
#
#                    IDENTITIES
#
######################################################


	def _get_identities(self):

		self._get_nt_identities()
		self._get_var_aa_identity()
		self._get_join_aa_identity()


	def _get_nt_identities(self):

		v_blocks = self._get_vdj_alignment_blocks(self.blocks, 'V')
		j_blocks = self._get_vdj_alignment_blocks(self.blocks, 'J')
		
		v = iter(v_blocks[0])
		row = next(v)
		while row == '': row = next(v)
		while row[0] != 'V': row = next(v)
		self.var_nt_identity = float(row.split()[1].replace('%', ''))
		
		if self.has_div:
			d_blocks = self._get_vdj_alignment_blocks(self.blocks, 'D')
			d = iter(d_blocks[0])
			row = next(d)
			while row == '': row = next(d)
			while row[0] != 'D': row = next(d)
			self.div_nt_identity = float(row.split()[1].replace('%', ''))
		
		j = iter(j_blocks[0])
		row = next(j)
		while row == '': next(j)
		while row[0] != 'J': row = next(j)
		self.join_nt_identity = float(row.split()[1].replace('%', ''))
		self.join_nt_length = int(row.split()[2].split('/')[1].replace(")",''))
		self.join_aa_length = int(math.floor(float(self.join_nt_length) / 3))


	def _get_var_aa_identity(self):
		matches = 0
		mismatches = 0
		[(a1, a2, score, begin, end)] = pairwise2.align.globalms(self.var_germ_aa_seq, self.var_aa_seq, 1, 0, -10, -1, one_alignment_only=True)
		if self.debug:
			print a1
			print a2
		for i, nt in enumerate(a1):
			if nt == '-' or a2[i] == '-':
				continue
			if nt == a2[i]:
				matches += 1
			else:
				mismatches += 1
		if self.var_ins:
			mismatches += len(self.var_ins)
		if self.var_del:
			mismatches += len(self.var_del)
		self.var_aa_identity = round(100. * matches / (matches + mismatches), 2)
		if self.debug:
			print self.var_aa_identity
		# self.var_aa_identity = round(100 * (float(len(self.var_aa_seq)) - self.var_muts_aa_count) / len(self.var_aa_seq), 2)


	def _get_join_aa_identity(self):

		germ = self._get_join_germ_seq(self.join_gene)
		self.join_aa = self.vdj_aa[-self.join_aa_length:]
		length = min(len(germ), len(self.join_aa))
		score = pairwise2.align.globalxx(self.join_aa, germ, one_alignment_only=1, score_only=1)
		self.join_aa_identity = round(100 * float(score) / length, 2)


	def _get_join_germ_seq(self, gene):

		trunc_gene = gene[:5]
		germs = {

			'human':	{	'IGHJ1' : 'AEYFQHWGQGTLVTVSS', 
							'IGHJ2' : 'YWYFDLWGRGTLVTVSS', 
							'IGHJ3' : 'DAFDVWGQGTMVTVSS', 
							'IGHJ4' : 'YFDYWGQGTLVTVSS',
							'IGHJ5' : 'NWFDSWGQGTLVTVSS',
							'IGHJ6' : 'YYYYYGMDVWGQGTTVTVSS',
							'IGKJ1' : 'WTFGQGTKVEIK',
							'IGKJ2' : 'YTFGQGTKLEIK',
							'IGKJ3' : 'FTFGPGTKVDIK',
							'IGKJ4' : 'LTFGGGTKVEIK',
							'IGKJ5' : 'ITFGQGTRLEIK',
							'IGLJ1' : 'YVFGTGTKVTVL',
							'IGLJ2' : 'VVFGGGTKLTVL',
							'IGLJ3' : 'VVFGGGTKLTVL',
							'IGLJ4' : 'FVFGGGTQLIIL',
							'IGLJ5' : 'WVFGEGTELTVL',
							'IGLJ6' : 'NVFGSGTKVTVL',
							'IGLJ7' : 'AVFGGGTQLTVL'},

			'macaque':	{	'IGHJ1' : 'AEYFEFWGQGALVTVSS',
							'IGHJ2' : 'YWYFDLWGPGTPITISS',
							'IGHJ3' : 'DAFDFWGQGLRVTVSS',
							'IGHJ4' : 'YFDYWGQGVLVTVSS',
							'IGHJ5' : 'NRFDVWGPGVLVTVSS',
							'IGHJ6' : 'YYGLDSWGQGVVVTVSS',
							'IGKJ1' : 'WTFGQGTKVEIK',
							'IGKJ2' : 'YSFGQGTKVEIK',
							'IGKJ3' : 'FTFGPGTKLDIK',
							'IGKJ4' : 'LTFGGGTKVEIK',
							'IGKJ5' : 'ITFGQGTRLEIK',
							'IGLJ1' : 'YIFGAGTRLTVL',
							'IGLJ2' : 'GLFGGGTRLTVL',
							'IGLJ3' : 'VLFGGGTRLTVL',
							'IGLJ4' : 'LVFCGGTQLTIV',
							'IGLJ5' : 'WVFGEGTKLTIL',
							'IGLJ6' : 'DVFGSGTKLTVL'},

			'mouse':	{	'IGHJ1' : 'YWYFDVWGAGTTVTVSS',
							'IGHJ2' : 'YFDYWGQGTTLTVSS',
							'IGHJ3' : 'WFAYWGQGTLVTVSA',
							'IGHJ4' : 'YYAMDYWGQGTSVTVSS',
							'IGKJ1' : 'WTFGGGTKLEIK',
							'IGKJ2' : 'YTFGGGTKLEIK',
							'IGKJ4' : 'FTFGSGTKLEIK',
							'IGKJ5' : 'LTFGAGTKLELK',
							'IGLJ1' : 'WVFGGGTKLTVL',
							'IGLJ2' : 'YVFGGGTKVTVL',
							'IGLJ3' : 'FIFGSGTKVTVL'},

			'rabbit':	{}}

		return germs[self.species][trunc_gene]




######################################################
#
#                     OUTPUT
#
######################################################


	def _get_log(self, log):

		if log == '':
			return sys.stdout
		else:
			open(log, 'w').write('')
			return open(log, 'a')


	def _insanity_output(self):

		return ['', '', '\n'.join(self.input)]


	def _exception_output(self, e):

		failed_seq_id = self.get_id()
		exc = '\n\n'
		exc += '-' * 60 
		exc += '\n'
		exc += failed_seq_id
		exc += '\n'
		exc += '-' * 60
		exc += '\n'
		exc += e
		exc += '\n'
		exc += '-' * 60 
		exc += '\n'
		exc += '\n'.join(self.input)
		exc += '\n' 
		exc += '\n\n'

		return ['', exc, '']


	def _build_output(self):

		# build the output for heavy chains
		if self.chain == 'heavy':
			if not self.tsv:
				output = collections.OrderedDict([
				('seq_id', self.seq_id),
				('uaid', self.uaid),
				('chain', self.chain),
				('v_gene', {'full' : self.var_gene, 
							'fam' : self.var_gene_fam, 
							'gene' : self.var_gene_gene, 
							'all' : self.var_gene_allele}),
 				('d_gene', {'full' : self.div_gene, 
 							'fam' : self.div_gene_fam, 
 							'gene' : self.div_gene_gene, 
 							'all' : self.div_gene_allele}),
 				('j_gene', {'full' : self.join_gene, 
 							'gene' : self.join_gene_gene, 
 							'all' : self.join_gene_allele}),
 				('bitscores', {'v' : self.var_bitscore, 
 							   'd' : self.div_bitscore, 
 							   'j' : self.join_bitscore}),
 				('e_values', {'v' : self.var_evalue, 
 							  'd' : self.div_evalue, 
 							  'j' : self.join_evalue}),
 				('nt_identity', {'v' : self.var_nt_identity, 
 								 'd' : self.div_nt_identity, 
 								 'j' : self.join_nt_identity}),
 				('aa_identity', {'v' : self.var_aa_identity, 
 								 'j' : self.join_aa_identity}),
 				('junc_aa', self.junction_aa),
 				('junc_len', len(self.junction_aa)),
 				('cdr3_aa', self.CDR3_aa),
 				('cdr3_len', len(self.CDR3_aa)),
 				('vdj_nt', self.vdj_nt),
 				('fr1_nt', self.FR1_nt),
 				('cdr1_nt', self.CDR1_nt),
 				('fr2_nt', self.FR2_nt),
 				('cdr2_nt', self.CDR2_nt),
 				('fr3_nt', self.FR3_nt),
 				('junc_nt', self.junction_nt),
 				('var_muts_nt', {'num': len(self.var_muts_nt), 
 								 'muts': [{'loc' : v.split(":")[0], 
 								 		   'mut' : v.split(":")[1]} for v in self.var_muts_nt]}),
 				('div_muts_nt', {'num': len(self.div_muts_nt), 
 								 'muts': [{'loc' : d.split(":")[0], 
 								 		   'mut' : d.split(":")[1]} for d in self.div_muts_nt]}),
 				('join_muts_nt', {'num': len(self.join_muts_nt), 
 								  'muts': [{'loc' : j.split(":")[0], 
 								  			'mut' : j.split(":")[1]} for j in self.join_muts_nt]}),
 				('mut_count_nt', len(self.var_muts_nt) + len(self.div_muts_nt) + len(self.join_muts_nt)),
 				('vdj_aa', self.vdj_aa),
 				('fr1_aa', self.FR1_aa),
 				('cdr1_aa', self.CDR1_aa),
 				('fr2_aa', self.FR2_aa),
 				('cdr2_aa', self.CDR2_aa),
 				('fr3_aa', self.FR3_aa),
 				('v_region_len_aa', {	'fr1' : int(self.FR1_aa_length), 
 										'cdr1' : int(self.CDR1_aa_length), 
 										'fr2' : int(self.FR2_aa_length), 
 										'cdr2' : int(self.CDR2_aa_length), 
 										'fr3' : int(self.FR3_aa_length)}),
 				('v_region_len_nt', {	'fr1' : int(self.FR1_nt_length), 
 										'cdr1' : int(self.CDR1_nt_length), 
 										'fr2' : int(self.FR2_nt_length), 
 										'cdr2' : int(self.CDR2_nt_length), 
 										'fr3' : int(self.FR3_nt_length)}),
 				('var_muts_aa', {'num': len(self.var_muts_aa), 'muts': [{'loc' : int(v.split(":")[0]), 'mut' : v.split(":")[1]} for v in self.var_muts_aa]}),
 				('v_ins', [{'loc' : int(i.split(",")[0]), 'len' : int(i.split(",")[1]), 'seq' : i.split(",")[2]} for i in self.var_ins]),
 				('v_del', [{'loc' : int(d.split(",")[0]), 'len' : int(d.split(",")[1]), 'seq' : d.split(",")[2]} for d in self.var_del]),
 				('d_ins', [{'loc' : int(i.split(",")[0]), 'len' : int(i.split(",")[1]), 'seq' : i.split(",")[2]} for i in self.div_ins]),
 				('d_del', [{'loc' : int(d.split(",")[0]), 'len' : int(d.split(",")[1]), 'seq' : d.split(",")[2]} for d in self.div_del]),
 				('j_ins', [{'loc' : int(i.split(",")[0]), 'len' : int(i.split(",")[1]), 'seq' : i.split(",")[2]} for i in self.join_ins]),
 				('j_del', [{'loc' : int(d.split(",")[0]), 'len' : int(d.split(",")[1]), 'seq' : d.split(",")[2]} for d in self.join_del]),
				('v_region_muts_nt', {	'fr1' : {'num' : int(self.FR1_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR1_nt_muts]}, 
										'cdr1' : {'num' : int(self.CDR1_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR1_nt_muts]}, 
										'fr2' : {'num' : int(self.FR2_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR2_nt_muts]}, 
										'cdr2' : {'num' : int(self.CDR2_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR2_nt_muts]}, 
										'fr3' : {'num' : int(self.FR3_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR3_nt_muts]}}),
				('v_region_muts_aa', {	'fr1' : {'num' : int(self.FR1_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR1_aa_muts]}, 
										'cdr1' : {'num' : int(self.CDR1_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR1_aa_muts]}, 
										'fr2' : {'num' : int(self.FR2_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR2_aa_muts]}, 
										'cdr2' : {'num' : int(self.CDR2_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR2_aa_muts]}, 
										'fr3' : {'num' : int(self.FR3_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR3_aa_muts]}}),
				('prod', self.productivity.lower()),
				('raw_input', self.raw_input),
				('alignment', self.pretty_alignment)
				])
		
		else:			
			if not self.tsv:
				output = collections.OrderedDict([
				('seq_id', self.seq_id),
				('uaid', self.uaid),
				('chain', self.chain),
				('v_gene', {'full' : self.var_gene, 
							'fam' : self.var_gene_fam, 
							'gene' : self.var_gene_gene, 
							'all' : self.var_gene_allele}),
 				('j_gene', {'full' : self.join_gene, 
 							'gene' : self.join_gene_gene, 
 							'all' : self.join_gene_allele}), 				
 				('bitscores', {'v' : self.var_bitscore, 
 							   'j' : self.join_bitscore}),
 				('e_values', {'v' : self.var_evalue, 
 							  'j' : self.join_evalue}),
				('nt_identity', {'v' : self.var_nt_identity, 
								 'j' : self.join_nt_identity}),
				('aa_identity', {'v' : self.var_aa_identity, 
								 'j' : self.join_aa_identity}),
 				('junc_aa', self.junction_aa),
 				('junc_len', len(self.junction_aa)),
 				('cdr3_aa', self.CDR3_aa),
 				('cdr3_len', len(self.CDR3_aa)),
 				('vdj_nt', self.vdj_nt),
 				('fr1_nt', self.FR1_nt),
 				('cdr1_nt', self.CDR1_nt),
 				('fr2_nt', self.FR2_nt),
 				('cdr2_nt', self.CDR2_nt),
 				('fr3_nt', self.FR3_nt),
 				('junc_nt', self.junction_nt),
 				('var_muts_nt', {'num': len(self.var_muts_nt), 
 								 'muts': [{'loc' : v.split(":")[0], 
 								 		   'mut' : v.split(":")[1]} for v in self.var_muts_nt]}),
 				('join_muts_nt', {'num': len(self.join_muts_nt), 
 								  'muts': [{'loc' : j.split(":")[0], 
 								  			'mut' : j.split(":")[1]} for j in self.join_muts_nt]}),
 				('mut_count_nt', len(self.var_muts_nt) + len(self.join_muts_nt)),
 				('vdj_aa', self.vdj_aa),
 				('fr1_aa', self.FR1_aa),
 				('cdr1_aa', self.CDR1_aa),
 				('fr2_aa', self.FR2_aa),
 				('cdr2_aa', self.CDR2_aa),
 				('fr3_aa', self.FR3_aa),
 				('v_region_len_aa', {	'fr1' : int(self.FR1_aa_length), 
 										'cdr1' : int(self.CDR1_aa_length), 
 										'fr2' : int(self.FR2_aa_length), 
 										'cdr2' : int(self.CDR2_aa_length), 
 										'fr3' : int(self.FR3_aa_length)}),
 				('v_region_len_nt', {	'fr1' : int(self.FR1_nt_length), 
 										'cdr1' : int(self.CDR1_nt_length), 
 										'fr2' : int(self.FR2_nt_length), 
 										'cdr2' : int(self.CDR2_nt_length), 
 										'fr3' : int(self.FR3_nt_length)}),
 				('var_muts_aa', {'num': len(self.var_muts_aa), 'muts': [{'loc' : int(v.split(":")[0]), 'mut' : v.split(":")[1]} for v in self.var_muts_aa]}),
 				('v_ins', [{'loc' : int(i.split(",")[0]), 'len' : int(i.split(",")[1]), 'seq' : i.split(",")[2]} for i in self.var_ins]),
 				('v_del', [{'loc' : int(d.split(",")[0]), 'len' : int(d.split(",")[1]), 'seq' : d.split(",")[2]} for d in self.var_del]),
 				('j_ins', [{'loc' : int(i.split(",")[0]), 'len' : int(i.split(",")[1]), 'seq' : i.split(",")[2]} for i in self.join_ins]),
 				('j_del', [{'loc' : int(d.split(",")[0]), 'len' : int(d.split(",")[1]), 'seq' : d.split(",")[2]} for d in self.join_del]),
				('v_region_muts_nt', {	'fr1' : {'num' : int(self.FR1_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR1_nt_muts]}, 
										'cdr1' : {'num' : int(self.CDR1_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR1_nt_muts]}, 
										'fr2' : {'num' : int(self.FR2_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR2_nt_muts]}, 
										'cdr2' : {'num' : int(self.CDR2_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR2_nt_muts]}, 
										'fr3' : {'num' : int(self.FR3_nt_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR3_nt_muts]}}),
				('v_region_muts_aa', {	'fr1' : {'num' : int(self.FR1_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR1_aa_muts]}, 
										'cdr1' : {'num' : int(self.CDR1_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR1_aa_muts]}, 
										'fr2' : {'num' : int(self.FR2_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR2_aa_muts]}, 
										'cdr2' : {'num' : int(self.CDR2_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.CDR2_aa_muts]}, 
										'fr3' : {'num' : int(self.FR3_aa_muts_count), 'muts' : [{'loc' : int(m.split(":")[0]), 'mut' : m.split(":")[1]} for m in self.FR3_aa_muts]}}),
				('prod', self.productivity.lower()),
				('raw_input', self.raw_input),
				('alignment', self.pretty_alignment)
				])

		for i in output.keys():
			if output[i] == "":
				del output[i]
			elif output[i] == []:
				del output[i]
			elif output[i] == {}:
				del output[i]

		json_output = json.dumps(output) + '\n'
		
		return [json_output, '', '']




