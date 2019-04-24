#*********************
#   Fei	Tan
#   ft54@njit.edu
#********************

## import modules 
from __future__ import print_function
import numpy as np
from keras.models import load_model
from keras import backend as K
from theano import tensor as T
from motif_viz import compile_saliency_function
import theano
import h5py
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import math
_EPSILON = K.epsilon()

#***********************
# configuration 
# data and model loading
#***********************
## configuration
model_file = '../../model/30/bestmodel1.hdf5'
seq_file_train = '../data/frac_m6A_30_train.pos_1.fa.hdf5'
seq_file_valid = '../data/frac_m6A_30_valid.pos_1.fa.hdf5'
seq_file_test = '../data/frac_m6A_30_test.pos_1.fa.hdf5'


##load data
trainmat = h5py.File(seq_file_train, 'r')
X_train = np.transpose(np.array(trainmat['x_train']),axes=(0,2, 1))
y_train = np.array(trainmat['y_train'])
frac_train = np.array(trainmat['frac'])

validmat = h5py.File(seq_file_valid, 'r')
X_valid = np.transpose(np.array(validmat['x_train']),axes=(0,2, 1))
y_valid = np.array(validmat['y_train'])
frac_valid = np.array(validmat['frac'])

testmat = h5py.File(seq_file_test, 'r')
X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])
frac_test = np.array(testmat['frac'])

X_test = np.concatenate((X_train, X_valid, X_test), axis = 0)
y_test = np.concatenate((y_train, y_valid, y_test), axis = 0)
frac_test = np.concatenate((frac_train, frac_valid, frac_test), axis = 0)

print('x shape:',X_test.shape)
print('y shape:',y_test.shape)
print('frac shape:',frac_test.shape)

del trainmat, validmat, testmat
del X_train, y_train, X_valid, y_valid, frac_train, frac_valid

idx = np.where((frac_test > 0.8) & (frac_test <= 1))
X_test = X_test[idx]
y_test = y_test[idx]
frac_test = frac_test[idx]
level = 'high'
print(y_test.shape)

## load keras model 
model = load_model(model_file)
model.summary()

#*******************
# compute salient 
# motif in sequence
#*******************
n = 2
idxs = np.array_split(range(X_test.shape[0]), n)

for i in range(len(idxs)):
	print(i)
	idx = idxs[i]
	X_input = X_test[idx,:,:]
	rslt_saliency_i, rslt_class_i = compile_saliency_function(model)([X_input, 0])
	if i == 0:
		rslt_saliency = rslt_saliency_i
		rslt_class = rslt_class_i
	else:
		rslt_saliency = np.concatenate((rslt_saliency, rslt_saliency_i), axis = 0)
		rslt_class = np.concatenate((rslt_class, rslt_class_i), axis = 0)


#***********************
#compute
#sequence frequency 
#**********************
seq_freq = np.zeros(X_test.shape[1:])
seq_freq_num = np.zeros(X_test.shape[1:])
for i in range(rslt_saliency.shape[0]):
	if rslt_class[i] == 1 and y_test[i] ==1 :
	#if rslt_class[i] == 1:
		saliency = rslt_saliency[i,:,:]
		input = X_test[i,:,:]
		saliency_input = saliency * input
		saliency_input_max = saliency_input.max(1)
		input_max = input.max(1)
                #idx1 = np.where(saliency_input_max > 0)
                #idx1 = idx1[0].tolist()
                idx2 = np.where(input_max == 1)
                idx2 = idx2[0].tolist()
                #idx = list(set(idx1).intersection(idx2))
		idx = idx2
		seq_freq[idx,:] += X_test[i,idx,:] * rslt_saliency[i,idx,:]
		seq_freq_num[idx,:] += X_test[i,idx,:]
		#print(i)


seq_freq = seq_freq/seq_freq_num
seq_freq = seq_freq.clip(min=0)
seq_freq[np.isnan(seq_freq)] = 0 

print(seq_freq)
print(seq_freq.shape)

seq_freq = np.transpose(seq_freq)
seq_freq_laplace = seq_freq + 0.0000000001
seq_freq_relative = seq_freq_laplace / seq_freq_laplace.sum(0) 


#*************************************
# plot position importance
# weighted by saliency map with motifs
#*************************************
seq_freq_cum = seq_freq.sum(axis = 0)
seq_freq_cum_norm = seq_freq_cum/seq_freq_cum.max()
position_score_motif = seq_freq_relative * seq_freq_cum_norm

prefix = level
np.savetxt(prefix + '.txt', position_score_motif, delimiter='\t')
Rcmd = 'Rscript plot_position_score_motif.r %s' % prefix
subprocess.call(Rcmd, shell = True)
