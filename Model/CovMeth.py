#************************
# Deep Learning 
# for DNA Methylation on
# N6-Adenine
# Fei Tan
# ft54@njit.edu
#************************


#*************************
# import modules 
# 
#*************************
import os
#os.environ['THEANO_FLAGS'] = "device=gpu0"
import sys
sys.setrecursionlimit(15000)
import numpy as np
import h5py
from sklearn.metrics import roc_auc_score, matthews_corrcoef, precision_recall_fscore_support
from pandas import DataFrame

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
from residual_blocks import building_residual_block
from keras.models import load_model
from keras.layers.advanced_activations import LeakyReLU, PReLU
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#************************************
# load data
#************************************
np.random.seed(1337) # for reproducibility
print 'loading data'
file_name = '../seq/m6A_30'
trainmat = h5py.File(file_name + '.train_1.hdf5', 'r')
validmat = h5py.File(file_name + '.valid_1.hdf5', 'r')
testmat = h5py.File(file_name + '.test_1.hdf5', 'r')
#testmat = h5py.File('/cstor/xsede/users/xs-ttgump/C.elegan/seq_neighbour/seq/m6A_30.test_1.hdf5', 'r')


X_train = np.transpose(np.array(trainmat['x_train']),axes=(0,2, 1))
y_train = np.array(trainmat['y_train'])

X_valid = np.transpose(np.array(validmat['x_train']),axes=(0,2, 1))
y_valid = np.array(validmat['y_train'])

X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])


print(X_train.shape)
print('train_label: count', np.unique(y_train,  return_counts=True))
print('valid_label: count', np.unique(y_valid,  return_counts=True))
print('test_label: count', np.unique(y_valid,  return_counts=True))


#**************************************
# build model
# convoluation
# residual 
#**************************************

print 'building model...............'
model = Sequential()

#*************************
# 1 convolutional layer
#*************************
NUM_FILTER1 = 80
model.add(Convolution1D(input_dim=4,
                        input_length=X_train.shape[1],
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

"""
#*******************************
# 6 convolutional layer
#*******************************
model.add(Convolution1D(nb_filter=40,
                        filter_length=4,
                        border_mode="valid",
                        activation="linear",
                        subsample_length=1, init='he_normal',
			name = "cov6"))

#model.add(Activation('relu'))
model.add(LeakyReLU(alpha=.001))
model.add(MaxPooling1D(pool_length=2, stride=2))
model.add(Dropout(0.2))
"""
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


#***********************************
# model training  
#
#
#***********************************
model.summary()
print 'compiling and fitting model...........'


#************************
# use pretrained weights
#************************
"""model_old = load_model('bestmodel_configure.hdf5')
layer_dict = dict([(layer.name, layer) for layer in model_old.layers])
for i in layer_dict.keys():
	 try:
	 	weight_old = layer_dict[i].get_weights()
	 	model.get_layer(i).set_weights(weight_old)
	 except:
		pass
         print i
del model_old
"""
sgd = SGD(lr=0.01, momentum=0.9, decay=1e-6, nesterov=True)
checkpointer = ModelCheckpoint(filepath="./bestmodel.hdf5", verbose=1, save_best_only=True)
earlystopper = EarlyStopping(monitor='val_loss', patience=50, verbose=1)

model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])

Hist = model.fit(X_train, y_train, batch_size=256, nb_epoch=500, shuffle=True, 
	  validation_data=(X_valid, y_valid), callbacks=[checkpointer,earlystopper])

print 'training done!'

#*******************
# plot both training
# and validation loss
#*******************
"""
k = 50
loss = Hist.history['loss']
val_loss = Hist.history['val_loss']
epoch = range(1,len(loss)+1)
plt.plot(epoch[k:], loss[k:])
plt.plot(epoch[k:], val_loss[k:])
plt.legend(['train_loss', 'valid_loss'], loc = 'upper right')
plt.xlabel('epoch')
plt.ylabel('loss')
plt.savefig('monitor.png')
"""



#***********************************
# prediction on validation and 
# test data
#***********************************
model = load_model('./bestmodel.hdf5')

print '**************vadiation results on validation dataset************'
keras_valid_valid = model.evaluate(X_valid, y_valid)
print keras_valid_valid

pred_prob_valid = model.predict(X_valid, verbose=1)
pred_class_valid = model.predict_classes(X_valid, verbose=1)


auc_valid = roc_auc_score(y_valid, pred_prob_valid)
mcc_valid = matthews_corrcoef(y_valid, pred_class_valid)
prfs_valid = precision_recall_fscore_support(y_valid, pred_class_valid)

print '************************'
print 'auc:', auc_valid
print 'mcc:', mcc_valid
print 'negative ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_valid[0][0], prfs_valid[1][0], prfs_valid[2][0], prfs_valid[3][0])
print 'positive ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_valid[0][1], prfs_valid[1][1], prfs_valid[2][1], prfs_valid[3][1])
print '************************'




print '**************prediction results on test dataset************'
keras_eval_test = model.evaluate(X_test, y_test)
print keras_eval_test

pred_prob_test = model.predict(X_test, verbose=1)
pred_class_test = model.predict_classes(X_test, verbose=1)


auc_test = roc_auc_score(y_test, pred_prob_test)
mcc_test = matthews_corrcoef(y_test, pred_class_test)
prfs_test = precision_recall_fscore_support(y_test, pred_class_test)

print '************************'
print 'auc:', auc_test
print 'mcc:', mcc_test
print 'negative ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_test[0][0], prfs_test[1][0], prfs_test[2][0], prfs_test[3][0])
print 'positive ---> precision:%s, recall:%s, f1score:%s, support:%s' %(prfs_test[0][1], prfs_test[1][1], prfs_test[2][1], prfs_test[3][1])
print '************************'



df = DataFrame({'true':y_test.flatten().tolist(),  'pred' :pred_prob_test.flatten().tolist()}, 
		index = range(len(y_test.flatten().tolist())))
df.to_csv('rslt.csv', index = False)
