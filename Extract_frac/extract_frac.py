#***************************************
# extract frac information from .fa file
# Xiurui Hou
# xh256@njit.edu
#***************************************
from DNA_IO import *;
import h5py
import sys
from joblib import Parallel, delayed
import multiprocessing

def getMatrix(file_name, chr_pos_frac ):
    seq_vecs = hash_sequences_1hot(file_name)
    seq_headers = seq_vecs.keys()
    train_seqs = []
    train_frac = []
    train_scores = []
    for header in seq_headers:
        train_seqs.append(seq_vecs[header])
        train_scores.append(1)
    count1 = 0
    count2 = 0
    index = 0
    for header in seq_headers:
        seq_range = header.split(':')[1]
        number = seq_range.split('-')
        chr_name = header.split(':')[0].replace('>','')
        #find the middle position
        if len(number) == 3:
            middle_pos = 30 - int(number[1])
        else:
            middle_pos = 30 + int(number[0])

        if chr_pos_frac[chr_name][middle_pos] > 1:
            train_frac.append('na')
            count2 += 1
        else:
            train_frac.append(chr_pos_frac[chr_name][middle_pos])
            count1 += 1
    print file_name + ' : ' + str(count1) + ' has frac'
    print file_name + ' : ' + str(count2) + ' has no frac\n'
    train_seqs = np.array(train_seqs)
    train_scores = np.array(train_scores)
    train_frac = np.array(train_frac)

    h5f = h5py.File('./frac/frac_' + file_name+'.hdf5', 'w')
    h5f.create_dataset('x_train', data=train_seqs)
    h5f.create_dataset('y_train', data=train_scores)
    h5f.create_dataset('frac', data=train_frac)
    h5f.close()
chr_pos_frac = {}
file = open('ara_m6A.gff', 'r')
line = file.readline()
count = 0
count1 = 0
count2 = 0
while line:
    count += 1
    position = line.split()[3]
    chr_name = line.split()[0]
    frac = line.split(';')[3].split('=')[1]
    if line.split(';')[3].split('=')[0] == 'frac':
        count1 += 1
        if chr_name not in chr_pos_frac.keys():
            chr_pos_frac[chr_name] = {}
        chr_pos_frac[chr_name][int(position)] = float(frac)
    else:
        count2 += 1
        if chr_name not in chr_pos_frac.keys():
            chr_pos_frac[chr_name] = {}
        chr_pos_frac[chr_name][int(position)] = 1.1
    line = file.readline()

print 'total : ' + str(count)
print 'has frac : ' + str(count1)
print 'has no frac : ' + str(count2) + '\n'

getMatrix('m6A_30_test.pos_1.fa', chr_pos_frac)
getMatrix('m6A_30_valid.pos_1.fa', chr_pos_frac)
getMatrix('m6A_30_train.pos_1.fa', chr_pos_frac)



