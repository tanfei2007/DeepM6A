import numpy as np
import h5py
import pandas as pd
from sklearn.manifold import TSNE
from keras import backend as K
from keras.models import load_model

# Read testing data and model
testmat = h5py.File('data/m6A_30.test_1.hdf5', 'r')
X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])

model = load_model('model/bestmodel_configure_30.hdf5')

# Select layer to plot
layer_name = 'dropout_6'
get_layer_output = K.function([model.layers[0].input, K.learning_phase()], [model.get_layer(name=layer_name).output])
# If the output of a layer is a 3-D array, it should be transfer to 2-D array first
layer_output_array = get_layer_output([X_test, 0])[0]

feat_cols = [ 'seq'+str(i) for i in range(layer_output_array.shape[1]) ]

df = pd.DataFrame(layer_output_array,columns=feat_cols)
df['label'] = y_test
df['label'] = df['label'].apply(lambda i: str(i))

# Sample 10000 sequences to plot
np.random.seed(1337) # for reproducibility
n_sne = 10000
rndperm = np.random.permutation(df.shape[0])


tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(df.loc[rndperm[:n_sne],feat_cols].values)

df_tsne = df.loc[rndperm[:n_sne],'label'].copy()
df_tsne['x-tsne'] = tsne_results[:,0]
df_tsne['y-tsne'] = tsne_results[:,1]

df_tsne.to_csv('Dropout_6_layer_tsne.csv', index = False)