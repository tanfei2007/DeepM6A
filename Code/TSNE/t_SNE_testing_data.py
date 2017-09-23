import numpy as np
import h5py
import pandas as pd
from sklearn.manifold import TSNE

# Read testing data
testmat = h5py.File('data/m6A_30.test_1.hdf5', 'r')
X_test = np.transpose(np.array(testmat['x_train']),axes=(0,2, 1))
y_test = np.array(testmat['y_train'])
# Transfer 3-D array to 2-D array
X_test_array = X_test.reshape((X_test.shape[0], -1), order='F')

feat_cols = [ 'seq'+str(i) for i in range(X_test_array.shape[1]) ]

df = pd.DataFrame(X_test_array,columns=feat_cols)
df['label'] = y_test
df['label'] = df['label'].apply(lambda i: str(i))

# Sample 10000 sequences to plot
np.random.seed(1337) # for reproducibility
n_sne = 10000
rndperm = np.random.permutation(df.shape[0])


tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(df.loc[rndperm[:n_sne],feat_cols].values)

df_tsne = df.loc[rndperm[:n_sne],:].copy()
df_tsne['x-tsne'] = tsne_results[:,0]
df_tsne['y-tsne'] = tsne_results[:,1]

df_tsne.to_csv('Raw_test_seq_tsne.csv', index = False)
