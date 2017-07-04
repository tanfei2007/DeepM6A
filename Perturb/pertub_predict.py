#************************
# Deep Residual Learning 
# for DNA Methylation
# Fei Tan
# ft54@njit.edu
#************************


#*************************
# import modules 
# 
#*************************
import os
os.environ['THEANO_FLAGS'] = "device=gpu0"
import sys
sys.setrecursionlimit(15000)
import numpy as np
import h5py
from sklearn.metrics import roc_auc_score, matthews_corrcoef, precision_recall_fscore_support
from pandas import DataFrame
from random import randint

from keras.preprocessing import sequence
from keras.optimizers import RMSprop
from keras.optimizers import SGD
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, MaxoutDense
from keras.layers.convolutional import Convolution1D, MaxPooling1D, AveragePooling1D
from keras.layers.local import LocallyConnected1D
from keras.layers.pooling import AveragePooling1D
from keras.regularizers import l2, activity_l1
from keras.constraints import maxnorm
from keras.layers.recurrent import LSTM, GRU
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers.wrappers import Bidirectional
from keras.regularizers import l1,l2
from keras.layers.normalization import BatchNormalization
from keras.models import load_model
from keras.layers.advanced_activations import LeakyReLU, PReLU
#from keras.utils.visualize_util import plot


#***********************************
# prediction on test data
#***********************************
model = Sequential()

#*************************
# 1 convolutional layer
#*************************
NUM_FILTER1 = 80
model.add(Convolution1D(input_dim=4,
                        input_length=61,
                        nb_filter=NUM_FILTER1,
                        filter_length=4,
                        border_mode="valid",
			activation="linear",
                        subsample_length=1,
			#W_regularizer = l2(0.01),
			init='he_normal',
			name = "cov1"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
#model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.2))

#*******************************
# 2 convolutional layer
#*******************************
model.add(Convolution1D(nb_filter=80,
                        filter_length=2,
                        border_mode="valid",
			activation="linear",
                        subsample_length=1, init='he_normal',
			name = "cov2"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
#model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.2))

#*******************************
# 3 convolutional layer
#*******************************

model.add(Convolution1D(nb_filter=80,
                        filter_length=4,
                        border_mode="valid",
                        activation="linear",
                        subsample_length=1, init='he_normal',
			name="cov3"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
#model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.2))

#*******************************
# 4 convolutional layer
#*******************************
model.add(Convolution1D(nb_filter=80,
                        filter_length=4,
                        border_mode="valid",
                        activation="linear",
                        subsample_length=1, init='he_normal',
			name = "cov4"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
#model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.2))


#*******************************
# 5 convolutional layer
#*******************************
model.add(Convolution1D(nb_filter=80,
                        filter_length=4,
                        border_mode="valid",
                        activation="linear",
                        subsample_length=1, init='he_normal',
			name = "cov5"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
#model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.5))


#**********************************
# post layers
#**********************************

model.add(Flatten())

#*****
# FC1
#*****
model.add(Dense(output_dim=100, init='he_normal'))
#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
model.add(Dropout(0.5))

model.add(Dense(output_dim=1))
model.add(Activation('sigmoid'))

model = load_model('./bestmodel_configure_30.hdf5')

model.summary()
#************************************
# load data
#************************************
np.random.seed(1337) # for reproducibility
print 'loading data'
file_name = '/cstor/xsede/users/xs-ttgump/C.elegan/seq/seq/m6A_30.test_'


idxs = range(1, 11)
#idxs = [1]
n_sample = len(idxs)
startCoors = [0, 1, 29, 29, 34, 41]
endCoors = [0, 28, 40, 33, 40, 61]
#startCoors = [1]
#endCoors = [28]
n_region = len(startCoors)


auc = np.zeros((n_sample, n_region))
mcc = np.zeros((n_sample, n_region))

for idx in idxs:
  testmat = h5py.File(file_name + str(idx) + '.hdf5', 'r')
  X = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
  y = np.array(testmat['y_train'])
  print('test_label: count', np.unique(y,  return_counts=True))

  #************************************
  # perturb data
  #************************************
  # letters = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
  
  for count in range(0, n_region):

    X_test = np.copy(X)
    y_test = np.copy(y)
    startCoor = startCoors[count]-1
    endCoor = endCoors[count]-1
    
    for i in range(0, X_test.shape[0]):
      for j in range(startCoor, endCoor+1):
        if j != (X_test.shape[1]-1)/2:
          if (X_test[i, j] == [1, 0, 0, 0]).all():
                letters = [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                X_test[i, j] = letters[randint(0,2)]
          elif (X_test[i, j] == [0, 1, 0, 0]).all():
                letters = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                X_test[i, j] = letters[randint(0,2)]
          elif (X_test[i, j] == [0, 0, 1, 0]).all():
                letters = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
                X_test[i, j] = letters[randint(0,2)]
          elif (X_test[i, j] == [0, 0, 0, 1]).all():
                letters = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
                X_test[i, j] = letters[randint(0,2)]


    print '**************prediction results on test dataset************'
    pred_prob_test = model.predict(X_test, verbose=1, batch_size = 100)
    pred_class_test = model.predict_classes(X_test, verbose=1, batch_size = 100)
    auc_test = roc_auc_score(y_test, pred_prob_test)
    mcc_test = matthews_corrcoef(y_test, pred_class_test)
    #prfs_test = precision_recall_fscore_support(y_test, pred_class_test)

    auc[idx-1, count] = auc_test
    mcc[idx-1, count] = mcc_test

startCoors = map(str, startCoors)
endCoors = map(str, endCoors)
col_names = ['_'.join(x) for x in zip(startCoors, endCoors)]
auc_df = DataFrame(data=auc, index=range(n_sample), columns=col_names)
mcc_df = DataFrame(data=mcc, index=range(n_sample), columns=col_names)


auc_df.to_csv('auc.csv', index = False)
mcc_df.to_csv('mcc.csv', index = False)


##print '************************'
##print 'auc:', auc_test
##print 'mcc:', mcc_test
##print 'negative ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_test[0][0], prfs_test[1][0], prfs_test[2][0], prfs_test[3][0])
##print 'positive ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_test[0][1], prfs_test[1][1], prfs_test[2][1], prfs_test[3][1])
##print '************************'





