#*****************
# import modules
#****************
from __future__ import division
import pandas as pd
import numpy as np
import h5py
import json
import itertools
#******************
# load data
#******************
file_name = '../seq/m6A_30'
trainmat = h5py.File(file_name + '.train_1.hdf5', 'r')
validmat = h5py.File(file_name + '.valid_1.hdf5', 'r')
testmat = h5py.File(file_name + '.test_1.hdf5', 'r')

X_train = np.transpose(np.array(trainmat['x_train']),axes=(0,2, 1))
y_train = np.array(trainmat['y_train'])
X_valid = np.transpose(np.array(validmat['x_train']),axes=(0,2, 1))
y_valid = np.array(validmat['y_train'])
X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])


X = np.concatenate((X_train, X_valid), axis = 0)
X = np.concatenate((X, X_test), axis = 0)
y = np.concatenate((y_train, y_valid), axis = 0)
y = np.concatenate((y, y_test), axis = 0)

idx_pos = np.where(y == 1)
X_pos = X[idx_pos[0],:,:]
idx_neg	= np.where(y == 0)
X_neg =	X[idx_neg[0],:,:]

del trainmat, validmat, testmat, X_train, X_valid, X_test, y_train, y_valid, y_test

#**************
#analysis
#**************

def letter_to_vec(letter):
	if letter == 'A':
		vec = [1,0,0,0]
	elif letter == 'C':
		vec = [0,1,0,0]
	elif letter == 'G':
		vec = [0,0,1,0]
	elif letter == 'T':
		vec = [0,0,0,1]
	
	return vec

def mletters_to_array(letters):
	n = len(letters)
	arr = np.zeros((n, 4))
	
	i = 0
	for letter in letters:
		arr[i,:] = letter_to_vec(letter)
		i = i + 1
	
	return arr


def gene_motifs(length = 4):
        
	letters	= ['A', 'C', 'G', 'T']
	N = 4**length	
	
	motifs = [p for p in itertools.product(letters, repeat=length)]	
	motifs = np.array(motifs).reshape(N, length)
	
	return motifs


def gene_positions(start=28, end=39, gap=4):
	
	N = end-start-gap+2
	positions = np.zeros((N, gap), dtype = int)
	
	for i in range(N):
		positions[i,:] = range(start+i,start+i+gap)
	
	return positions


def odds_ratio(X_pos, X_neg, positions, letters):
	"""
	positions = [29, 30, 31, 32]
	letters = ['G', 'A', 'G', 'G']
	"""
	n_pos = X_pos.shape[0]
	n_neg = X_neg.shape[0]
	

	#number of methylated sites with motif
	pos_with_motif = (X_pos[:, positions, :] ==  mletters_to_array(letters)).all(axis=(1,2)).sum()
	
	#number of methylated sites without motif
        pos_without_motif = n_pos - pos_with_motif
	
	#number of non-methylated sites with motif
        neg_with_motif = (X_neg[:, positions, :] ==  mletters_to_array(letters)).all(axis=(1,2)).sum()

	#number of non-methylated sites without motif
	neg_without_motif = n_neg - neg_with_motif

	#print str(pos_with_motif), str(pos_without_motif), str(neg_with_motif), str(neg_without_motif)

	if (pos_with_motif == 0) or (neg_with_motif == 0):
		oddr = 0
	else:
		oddr = (pos_with_motif/pos_without_motif) / (neg_with_motif/neg_without_motif)
	
	return oddr, pos_with_motif, neg_with_motif


gaps = [4,3]
for gap in gaps:
	positions = gene_positions(start=28, end=39, gap=gap)
	motifs = gene_motifs(gap)

	rslt = {}
	for i in range(positions.shape[0]):
		for j in range(motifs.shape[0]):
			position = positions[i,:]
			motif = motifs[j,:]
			key = ''.join(motif) + '_' + ''.join(map(str, position))
			rslt[key] = odds_ratio(X_pos, X_neg, position, motif)
			print i, j
	print(len(rslt))

	with open('motif_' + str(gap) + '.json', 'w') as fp:
    		json.dump(rslt, fp)






