#*********************
#   Fei Tan
#   ft54@njit.edu
#   March 28, 2017
#********************



from keras import backend as K


class Loss(object):
	
	def __init__(self):
		self.name = 'Unamed Loss'

	def __str__(self):
		return self.name
	
	def build_loss(self):
		raise NotImplementedError()




class ActivationMaximization(Loss):
	
	def __init__(self, layer, filter_indices):
		super(ActivationMaximization, self).__init__()
		self.name = 'ActivationMaximization Loss'
		self.layer = layer
		self.filter_indices = filter_indices

	def build_loss(self):
		layer_output = self.layer.output

		#for all other layers it is 3
		is_dense = K.ndim(layer_output) == 2
		
		loss = 0.

		for idx in self.filter_indices:
			if is_dense:
				loss += - K.mean(layer_output[:,idx])
			else:
				loss += -K.mean(layer_output[:,:,idx])
		
		return(loss)
