import numpy as np


def make_noise(X, y, fdr=0.5):
	print(y.shape)
	idx_pos = np.where(y==1)
	idx_neg = np.where(y==0)

	pos_len = int(np.ceil(idx_pos[0].shape[0]/2.0))
	false_pos_length = int(pos_len*fdr)
	true_pos_length = pos_len - false_pos_length

	idx_true_pos = idx_pos[0][:true_pos_length]
	idx_true_neg = idx_neg[0][:pos_len]
	idx_false_pos= idx_neg[0][pos_len:(pos_len+false_pos_length)]

	idx = np.concatenate([idx_true_pos, idx_true_neg, idx_false_pos])
	X = X[idx]
	y[idx_false_pos] = 1
	y = y[idx] 
	
	return X, y
