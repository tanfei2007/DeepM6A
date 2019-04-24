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
model_file = '../model/30/bestmodel1.hdf5'
seq_file_train = '../seq/hdf5/30/m6A_30.train_1.hdf5'
seq_file_valid = '../seq/hdf5/30/m6A_30.valid_1.hdf5'
seq_file_test = '../seq/hdf5/30/m6A_30.test_1.hdf5'


##load data
trainmat = h5py.File(seq_file_train, 'r')
X_train = np.transpose(np.array(trainmat['x_train']),axes=(0,2, 1))
y_train = np.array(trainmat['y_train'])

validmat = h5py.File(seq_file_valid, 'r')
X_valid = np.transpose(np.array(validmat['x_train']),axes=(0,2, 1))
y_valid = np.array(validmat['y_train'])


testmat = h5py.File(seq_file_test, 'r')
X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])


X_test = np.concatenate((X_train, X_valid, X_test), axis = 0)
y_test = np.concatenate((y_train, y_valid, y_test), axis = 0)


print('x shape:',X_test.shape)
print('y shape:',y_test.shape)

del trainmat, validmat, testmat
del X_train, y_train, X_valid, y_valid


## load keras model 
model = load_model(model_file)
model.summary()

#*******************
# compute salient 
# motif in sequence
#*******************
n = 20
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

#******************************
# plot culumlative saliency map
# via heat map
#*****************************
prefix = 'saliency_seq'
np.savetxt(prefix + '.txt', seq_freq_relative, delimiter='\t')
Rcmd = 'Rscript seq_logo_plot.r %s' % prefix
subprocess.call(Rcmd, shell = True)

#*******************************
# plot position importance
# weighted by saliency map
#*******************************
seq_freq_cum = seq_freq.sum(axis = 0)
seq_freq_cum_norm = seq_freq_cum/seq_freq_cum.max()

prefix = 'seq_freq_cum'
np.savetxt(prefix + '.txt', seq_freq_cum_norm, delimiter='\t')
Rcmd = 'Rscript seq_freq_plot.r %s' % prefix
subprocess.call(Rcmd, shell = True)


#*************************************
# plot position importance
# weighted by saliency map with motifs
#*************************************
seq_freq_cum = seq_freq.sum(axis = 0)
seq_freq_cum_norm = seq_freq_cum/seq_freq_cum.max()
position_score_motif = seq_freq_relative * seq_freq_cum_norm

prefix = 'position_score_motif'
np.savetxt(prefix + '.txt', position_score_motif, delimiter='\t')
Rcmd = 'Rscript plot_position_score_motif.r %s' % prefix
subprocess.call(Rcmd, shell = True)



""""""
#*********************
# functions 
# 
#********************

def m1hot_to_dec(seq_1hot):
	""" convert multiple
	one-hot encoding to decimal number
	"""
        rslt = ''
        for i in range(seq_1hot.shape[0]):
                rslt += format(int(np.where(seq_1hot[i,:] == 1)[0]),'02b')
        dec = int(rslt,2)
        return dec


def motif_scan(input_seq, idx_positive, size_motif, saliency_input):
        motif = np.zeros((X_test.shape[1]-size_motif+1, 4 ** size_motif))
        for idx in idx_positive:
                #check if the consecutive motifs lie in idx_positive
                motif_idx = range(idx, idx+size_motif)
                if set(motif_idx) < set(idx_positive):
                        combind_idx = m1hot_to_dec(input_seq[motif_idx,:])
                        motif[idx,combind_idx] = saliency_input[motif_idx,:].sum()
        return motif

def dec_to_string(idx, key):
        letters = {'00':'A', '01':'C', '10':'G', '11':'T'}
        spec = '0'+str(key*2)+'b'
        idxs = format(idx, spec)
        rslt = ''
        for i in range(key):
                rslt += letters[idxs[(i*2):(i*2+2)]]
        return rslt

def decs_to_string(idx_list, key):
        rslt = []
        for idx in idx_list:
                rslt.append(dec_to_string(idx, key))
        return rslt


def cal_IC(motif_freq_key, key):
	"""calculate quantaties associated with information content
	"""
        motif_freq_key = np.transpose(motif_freq_key)
	nums = np.sum(motif_freq_key, axis = 0) #number of samples located within each starting site
	s = motif_freq_key.shape[0] #the number of candidates in a specific position
	ens = (1/(math.log(2))) * ((s-1)/(2*nums))  #The approximation for the small-sample correction
	

        motif_freq_key +=  _EPSILON  #additive/laplace smoothing
        position_weight = np.sum(motif_freq_key, axis = 0)
        position_motif_weight = motif_freq_key/position_weight
	
	H = -position_motif_weight * np.log2(position_motif_weight)
	H = np.sum(H, axis = 0)
	R = math.log(s,2) - (H + ens) #total bits for each site
	R = R * position_weight/position_weight.max()
	R_motif = R * position_motif_weight	#bits/length for each motif at each site
	R_letter = R_motif/key 	#the average bits for letters in each motif at each site
	
	
        return(position_weight, position_motif_weight, motif_freq_key, R, R_motif, R_letter)


def plot_position_weight_heat(position_weight, file_name):
        weight_range = position_weight.max()
        #
	position_weight = position_weight.reshape(position_weight.shape[0],1)
        sns.set(font_scale=1)
        plt.figure(figsize=(5, position_weight.shape[0]))
        sns.heatmap(position_weight, cmap='PRGn', linewidths=0.1, vmin=0, vmax=weight_range)
        ax = plt.gca()
        ax.set_xticklabels([1])
        ylabels = range(1, position_weight.shape[0]+1)
        ylabels.reverse()
        ax.set_yticklabels(ylabels, rotation='horizontal')
        plt.savefig(file_name)
        plt.close()


def plot_position_motif_weight(position_motif_weight, file_name, key):
        motif_weight_range = position_motif_weight.max()
        motif_number = position_motif_weight.shape[0]
        #
        sns.set(font_scale=5)
        plt.figure(figsize=(position_motif_weight.shape[1]+10, motif_number))
        sns.heatmap(position_motif_weight, cmap='PRGn', linewidths=0.1, vmin=0, vmax=motif_weight_range)
        ax = plt.gca()
        #
	ylabels = decs_to_string(range(0, motif_number), key)
        ylabels.reverse()
        ax.set_xticklabels(range(1, position_motif_weight.shape[1]+1))
        ax.set_yticklabels(ylabels, rotation='horizontal')
        ax.tick_params(labeltop=True, labelright=True)
        #plt.colorbar(ax)
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()


def plot_all(motif_freq_key, key):
	position_weight, position_motif_weight, motif_freq_key, R, R_motif, R_letter = cal_IC(motif_freq_key, key = 4)
        app_str =  '_' + str(key) + '.pdf'
        plot_position_weight_heat(position_weight, 'position_weight' + app_str)
        #plot_position_motif_weight(position_motif_weight, 'position_motif_weight' + app_str, key)
        #plot_position_motif_weight(motif_freq_key, 'motif_freq_key'+app_str, key)

        plot_position_weight_heat(R, 'R' + app_str)
        #plot_position_motif_weight(R_motif, 'R_motif' + app_str, key)
        #plot_position_motif_weight(R_letter, 'R_letter' + app_str, key)

	#**************************
	# print out the most
	# salient motifs
	#**************************
        motif_num = R_motif.shape[0]
        site_len = R_motif.shape[1]

        ylabels = decs_to_string(range(0, motif_num), key)
        idx = np.dstack(np.unravel_index(np.flipud(np.argsort(R_motif.ravel())), (motif_num, site_len)))

        #print top 20
        top = 20
        file_name = open('top20_' + str(key) + '.txt', 'w')
        print(('Site', 'Motif', 'IC', 'Weight', 'Freq'), file = file_name)
        for i in range(top):

        	motif = ylabels[idx[0,i,0]]
        	site = idx[0,i,1] + 1
        	ic = R_motif[idx[0,i,0], idx[0,i,1]]
        	weight = position_motif_weight[idx[0,i,0], idx[0,i,1]]
       		freq = motif_freq_key[idx[0,i,0], idx[0,i,1]]
        	print((site, motif, ic, weight, freq), file = file_name)


#*********************
#compute motif
# frequency
#*********************

motif_freq={}
keys = [3]
for key in keys:
	motif_freq[key] = np.zeros((X_test.shape[1]-key+1, 4**key))

for i in range(rslt_saliency.shape[0]):
        if rslt_class[i] == 1 and y_test[i] == 1:
                saliency = rslt_saliency[i,:,:]
                input = X_test[i,:,:]
                saliency_input = saliency * input
                saliency_input_max = saliency_input.max(1)
		input_max = input.max(1)
                
		idx1 = np.where(saliency_input_max > 0)
		idx1 = idx1[0].tolist()
		idx2 = np.where(input_max == 1)
		idx2 = idx2[0].tolist()
		idx = list(set(idx1).intersection(idx2))


		for key in motif_freq.keys():
			temp = motif_scan(input, idx, key, saliency_input)
			motif_freq[key] += temp
                print(i)




#plot_all(motif_freq[1], key = 1)
#plot_all(motif_freq[3], key = 3)
#plot_all(motif_freq[5], key = 5)
""""""



