import numpy as np
import h5py
import pandas as pd
from sklearn.manifold import TSNE

from keras import backend as K

from sklearn.metrics import roc_auc_score
from keras.preprocessing import sequence
from keras.optimizers import RMSprop
from keras.optimizers import SGD
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.regularizers import l2, activity_l1
from keras.constraints import maxnorm
from keras.layers.recurrent import LSTM, GRU
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.models import load_model

testmat1 = h5py.File('data/frac_m6A_30_test.pos_1.fa.hdf5', 'r')
X_test1 = np.transpose(np.array(testmat1['x_train']),axes=(0,2, 1))
y_test1 = np.array(testmat1['y_train'])
frac_test1 = np.array(testmat1['frac'])

testmat2 = h5py.File('data/m6A_30.test_1.hdf5', 'r')
X_test2 = np.transpose(np.array(testmat2['x_train']),axes=(0,2, 1))
y_test2 = np.array(testmat2['y_train']).ravel()

X_test = np.concatenate((X_test1, X_test2[X_test2.shape[0]/2+1:,:,:]), axis=0)
y_test = np.concatenate((y_test1, y_test2[y_test2.shape[0]/2+1:]), axis=0)
frac_test = np.concatenate((frac_test1, np.zeros(X_test2[X_test2.shape[0]/2+1:,:,:].shape[0])), axis=0)

model = load_model('model/bestmodel_configure_30.hdf5')
layer_name = 'dropout_6'
get_layer_output = K.function([model.layers[0].input, K.learning_phase()], [model.get_layer(name=layer_name).output])
layer_output = get_layer_output([X_test, 0])[0]

# layer_output_array = layer_output.reshape((X_test.shape[0], -1), order='F')

layer_output_array = layer_output.copy()
feat_cols = [ 'seq'+str(i) for i in range(layer_output_array.shape[1]) ]

df = pd.DataFrame(layer_output_array,columns=feat_cols)
df['label'] = y_test
df['label'] = df['label'].apply(lambda i: str(i))
df['frac'] = frac_test
df['frac'] = df['frac'].apply(lambda i: str(i))

# Sample 10000 sequences to plot
np.random.seed(1337) # for reproducibility
n_sne = 10000
rndperm = np.random.permutation(df.shape[0])


tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(df.loc[rndperm[:n_sne],feat_cols].values)

df_tsne = df.loc[rndperm[:n_sne],:].copy()
df_tsne['x-tsne'] = tsne_results[:,0]
df_tsne['y-tsne'] = tsne_results[:,1]

df_tsne.to_csv('layer_dropout_6_seq_tsne1.csv', index = False)

df_tsne2 = df_tsne.loc[:,['label','frac','x-tsne','y-tsne']]
df_tsne2.to_csv('layer_dropout_6_seq_tsne2.csv', index = False)