#***************************************************************
# labeling sequence according to site information (Example: Ara)
# Xiurui Hou
# xh256@njit.edu
#***************************************************************

from DNA_IO import *;
import h5py
import sys
from joblib import Parallel, delayed
import multiprocessing
import csv

def getMatrix(num, nature, pos_nature):
    seq_vecs = hash_sequences_1hot('/cstor/xsede/users/xs-ttgump/Ara/genome/'+num+'.'+nature+'.fa')
    seq_headers = seq_vecs.keys()
    n = len(seq_headers)
    print num+'.'+nature+'.fa has ' + str(n) + ' seqs'

    train_seqs = []
    train_scores = []
    for header in seq_headers:
        train_seqs.append(seq_vecs[header])
    n2 = 0
    for header in seq_headers:
        seq_range = header.split(':')[1]
        number = seq_range.split('-')
        
        #find the middle position
        if len(number) == 3:
            middle_pos = 30 - int(number[1])
        else:
            middle_pos = 30 + int(number[0])

        #check the middle position and nature in csv file
        if str(middle_pos) in pos_nature:
            if (nature == 'neg' and pos_nature[str(middle_pos)]=='-') or (nature == 'pos' and pos_nature[str(middle_pos)]=='+'):
                n2 += 1
                train_scores.append(1)
            else:
                train_scores.append(0)
        else:
            train_scores.append(0)
    print num+'_'+nature+'.fa has ' + str(n2) + ' seqs in the csv'

    train_seqs = np.array(train_seqs)
    train_scores = np.array(train_scores)
    
    h5f = h5py.File('./h5_output/' + num + '.' + nature+'.hdf5', 'w')
    h5f.create_dataset('x_train', data=train_seqs)
    h5f.create_dataset('y_train', data=train_scores)
    h5f.close()

chr_num = ['NC_003070.9','NC_003071.7','NC_003074.8','NC_003075.7','NC_003076.8']

seq_chr = {}
seq_chr['NC_003070.9'] = 'chr1'
seq_chr['NC_003071.7'] = 'chr2'
seq_chr['NC_003074.8'] = 'chr3'
seq_chr['NC_003075.7'] = 'chr4'
seq_chr['NC_003076.8'] = 'chr5' 
for num in chr_num:
    csv_file = file('ara.csv', 'r')
    reader = csv.reader(csv_file)
    pos_nature = {}

    for line in reader:
        name, pos, nature, seq = line
        if name == seq_chr[num]:
            pos_nature[pos] = nature
    csv_file.close()
    print 'pos_nature length is : ' + str(len(pos_nature))
    getMatrix(num, 'pos', pos_nature)
    getMatrix(num, 'neg', pos_nature)

