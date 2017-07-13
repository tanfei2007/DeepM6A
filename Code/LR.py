#**************************
# Logistic Regression
# for DNA N6-Adenine Methylation
# Tian Tian
# tt72@njit.edu
#**************************

import sys
import numpy as np
import itertools
import multiprocessing

#**************************
# import modules
#**************************
from joblib import Parallel, delayed
from collections import OrderedDict
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import *
from sklearn.preprocessing import StandardScaler
from sklearn.externals import joblib
from numpy import genfromtxt

#**************************
# load data
#**************************
name = sys.argv[1]
index = sys.argv[2]
print 'loading kmer data'
train = h5py.File('kmer/LR_m6A_'+name+'.train_'+index+'.hdf5', 'r')
test = h5py.File('kmer/LR_m6A_'+name+'.test_'+index+'.hdf5', 'r')
valid = h5py.File('kmer/LR_m6A_'+name+'.valid_'+index+'.hdf5', 'r')

Y_train = np.array(train['Y'])
Y_test = np.array(test['Y'])
Y_valid = np.array(valid['Y'])

X_train = np.transpose(np.array(train['X']))
X_valid = np.transpose(np.array(valid['X']))

Y_train_array = np.concatenate((Y_train,Y_valid), axis=0)
Y_train_array = (Y_train_array.ravel()).tolist()

X_train_array = np.concatenate((X_train,X_valid), axis=0)
X_test_array = np.transpose(np.array(test['X']))
Y_test_array = (Y_test.ravel()).tolist()

#**************************
# Data standardization
# Same scaling was applied for both test and train data
#**************************
scaler = StandardScaler()
scaler.fit(X_train_array)
X_train_array = scaler.transform(X_train_array)
X_test_array = scaler.transform(X_test_array)

#**************************
# Logistic Regression
#**************************
print 'training LR model'
clf = LogisticRegressionCV(Cs=10.0**-np.arange(-4,4), cv=5, penalty="l2", solver="sag", n_jobs=-1, max_iter=2000, verbose=1)
LR_fit = clf.fit(X_train_array, Y_train_array)
predict_proba = LR_fit.predict_proba(X_test_array)
predict_class = LR_fit.predict(X_test_array)
df = DataFrame({'groundtruth':Y_test_array, 'predict_proba':predict_proba[:,1].flatten().tolist()}, index = range(len(predict_proba[:,1].flatten().tolist())))
df.to_csv('rslt_len'+str(length)+'_index'+str(indx)+'.csv', index = False)
false_positive_rate, true_positive_rate, thresholds = roc_curve(Y_test_array, predict_proba[:,1])
auc_test = auc(false_positive_rate, true_positive_rate)
mcc_test = matthews_corrcoef(Y_test_array, predict_class)
prf = precision_recall_fscore_support(Y_test_array, predict_class)

#**************************
# Output results
#**************************
print 'auc:', auc_test
print 'mcc:', mcc_test
print 'precision, recall, f1score:', prf

joblib.dump(clf, 'LR_len_'+name+'_'+index+'.pkl')