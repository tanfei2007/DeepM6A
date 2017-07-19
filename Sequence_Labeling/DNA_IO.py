##copied dna_io.py

from collections import OrderedDict

import numpy as np
import numpy.random as npr
from sklearn import preprocessing


################################################################################
# dna_one_hot
#
# Input
#  seq:
#
# Output
#  seq_vec: Flattened column vector
################################################################################

def dna_one_hot(seq, seq_len=None, flatten=False):
	if seq_len == None:
		seq_len = len(seq)
		seq_start = 0
	else:
		if seq_len <= len(seq):
			# trim the sequence
			seq_trim = (len(seq)-seq_len)/2
			seq = seq[seq_trim:seq_trim+seq_len]
			seq_start = 0
		else:
			seq_start = (seq_len-len(seq))/2
	
	seq = seq.upper()
	seq = seq.replace('A','0')
	seq = seq.replace('C','1')
	seq = seq.replace('G','2')
	seq = seq.replace('T','3')
	# map nt's to a matrix 4 x len(seq) of 0's and 1's.
	#  dtype='int8' fails for N's
	seq_code = np.zeros((4,seq_len), dtype='float16')
	for i in range(seq_len):
		if i < seq_start:
			seq_code[:,i] = 0.25
		else:
			try:
				seq_code[int(seq[i-seq_start]),i] = 1
			except:
				seq_code[:,i] = 0.25
	
	# flatten and make a column vector 1 x len(seq)
	if flatten:
		seq_vec = seq_code.flatten()[None,:]
		return seq_vec
	else:
		return seq_code
#

	
################################################################################
# hash_sequences_1hot
#
# Input
#  fasta_file:  Input FASTA file.
#  extend_len:  Extend the sequences to this length.
#
# Output
#  seq_vecs:    Dict mapping FASTA headers to sequence representation vectors.
################################################################################
def hash_sequences_1hot(fasta_file, extend_len=None):
	# determine longest sequence
	if extend_len is not None:
		seq_len = extend_len
	else:
		seq_len = 0
		seq = ''
		for line in open(fasta_file):
			if line[0] == '>':
				if seq:
					seq_len = max(seq_len, len(seq))
				#
				header = line[1:].rstrip()
				seq = ''
			else:
				seq += line.rstrip()
		#
		if seq:
			seq_len = max(seq_len, len(seq))
	#
	# load and code sequences
	seq_vecs = OrderedDict()
	seq = ''
	for line in open(fasta_file):
		if line[0] == '>':
			if seq:
				seq_vecs[header] = dna_one_hot(seq, seq_len)
			#
			header = line[1:].rstrip()
			seq = ''
		else:
			seq += line.rstrip()
	#
	if seq:
		seq_vecs[header] = dna_one_hot(seq, seq_len)
	#
	return seq_vecs

#

