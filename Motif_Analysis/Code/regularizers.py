import numpy as np
from keras import backend as K

from losses import Loss




def normalize(seq, value):
	
	return value / np.prod(seq._keras_shape[1:])



class LPNorm(Loss):
	def __init__(self, seq_input, p=2.):
		super(LPNorm, self).__init__()
		if p < 1:
			raise ValueError('p value should range between [1, inf)')
		self.name = "L-{} Norm Penalty".format(p)
		self.p = p
		self.seq = seq_input

	def build_loss(self):
		if np.isinf(self.p):
			value = K.max(self.seq)
		else:
            		value  = K.pow(K.sum(K.pow(K.abs(self.seq), self.p)), 1./self.p)
			
	    	return normalize(self.seq, value)
