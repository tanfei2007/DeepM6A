import h5py
import numpy as np

def gene_negative(size = 10000, length = 30):
	letters = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
	X = np.zeros((size,2*length+1,4))
	np.random.seed(526)
	for i in range(0,size):
		 letter = np.random.choice(range(0,4), 2*length+1, replace=True)
		 X[i,:,:] = letters[letter,:]
		 X[i,length,:] = letters[0,:]
	y = np.zeros(size)
	return(X, y)



def gene_positive(size = 10000, length = 30):
        letters = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
	#signal_letters = np.array([letters, letters[[1,2,3,0],:], letters[[2,3,0,1],:],letters[[3,0,1,2],:] ])
	signal_letters = np.zeros((1, 4, 4))
	for letter0 in range(0,4):
		for letter1 in range(0,4):
			if letter1 != letter0:
				for letter2 in range(0,4):
					if letter2 not in [letter0, letter1]:
						for letter3 in range(0,4):
							if letter3 not in [letter0, letter1, letter2]:
								print(letter0, letter1, letter2, letter3)
								one_hot_letter0 = [0,0,0,0]
								one_hot_letter1 = [0,0,0,0]
								one_hot_letter2 = [0,0,0,0]
								one_hot_letter3 = [0,0,0,0]
								print(one_hot_letter0)
								one_hot_letter0[letter0] = 1
								one_hot_letter1[letter1] = 1	
								one_hot_letter2[letter2] = 1
								one_hot_letter3[letter3] = 1
								print(one_hot_letter0)
								signal_letter = np.array([[one_hot_letter0, one_hot_letter1, one_hot_letter2, one_hot_letter3]])
								signal_letters = np.concatenate((signal_letters, signal_letter))
	signal_letters = signal_letters[1:,:,:]
	print(signal_letters)
	print(signal_letters.shape)		
	
	X = np.zeros((size,2*length+1,4))
	np.random.seed(526)
        for i in range(0,size):
                 letter = np.random.choice(range(0,4), 2*length+1, replace=True)
                 X[i,:,:] = letters[letter,:]
		 X[i,length,:] = letters[0,:]
		 
		 singal_letter = np.random.choice(range(0,24), 1, replace=True)
		 #if i%2 != 0:
		 #X[i,[length+1, length+3, length+5, length+7],:] = signal_letters[singal_letter,:,:]
		 #else:
		 X[i,[length-1, length-3, length-5, length-7],:] = letters

	y = np.ones(size)
	return (X, y)


#*****************
#generate data
#****************

#train
X_neg_train, y_neg_train = gene_negative(20000,30)
X_pos_train, y_pos_train = gene_positive(20000,30)
X_train = np.concatenate((X_neg_train, X_pos_train), axis = 0)
y_train = np.concatenate((y_neg_train, y_pos_train), axis = 0)

#valid
X_neg_valid, y_neg_valid = gene_negative(4000,30)
X_pos_valid, y_pos_valid = gene_positive(4000,30)
X_valid	= np.concatenate((X_neg_valid, X_pos_valid), axis = 0)
y_valid	= np.concatenate((y_neg_valid, y_pos_valid), axis = 0)

#test
X_neg_test, y_neg_test = gene_negative(4000,30)
X_pos_test, y_pos_test = gene_positive(4000,30)
X_test	= np.concatenate((X_neg_test, X_pos_test), axis = 0)
y_test	= np.concatenate((y_neg_test, y_pos_test), axis = 0)



h5f = h5py.File('./tgca_hdf5/train.hdf5','w')
h5f.create_dataset('x_train', data = X_train, compression="gzip", compression_opts=9)
h5f.create_dataset('y_train', data = y_train, compression="gzip", compression_opts=9)
h5f.close()

h5f = h5py.File('./tgca_hdf5/valid.hdf5','w')
h5f.create_dataset('x_train', data = X_valid, compression="gzip", compression_opts=9)
h5f.create_dataset('y_train', data = y_valid, compression="gzip", compression_opts=9)
h5f.close()

h5f = h5py.File('./tgca_hdf5/test.hdf5','w')
h5f.create_dataset('x_train', data = X_test, compression="gzip", compression_opts=9)
h5f.create_dataset('y_train', data = y_test, compression="gzip", compression_opts=9)
h5f.close()
