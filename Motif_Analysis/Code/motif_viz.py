#*********************
#   Fei Tan
#   ft54@njit.edu
#********************



## import modules 
import numpy as np
import subprocess
from keras.models import load_model
from keras import backend as K
from theano import tensor as T
import pprint
from losses import ActivationMaximization
from optimizer import Optimizer
from regularizers import LPNorm


def get_num_filters(layer):
	isDense = K.ndim(layer.outpu) == 2

	if isDense:
		return layer.output.shape[1]
	else:
		return layer.output._keras_shape[2]	


def viz_activation(model, layer_idx, filter_indices=None,
                         seed_seq=None, act_max_weight=1, lp_norm_weight=1, **optimizer_params):

	print("Working on filters: {}".format(pprint.pformat(filter_indices)))
	
	optimizer_params_default = {
		'seed_seq': seed_seq,
		'max_iter': 1000,
		'lr': 0.1,
		'verbose': True
	}
	
	optimizer_params_default.update(optimizer_params)
	optimizer_params = optimizer_params_default
	
	losses = [(ActivationMaximization(model.layers[layer_idx], filter_indices), act_max_weight),
		(LPNorm(model.input, 2), lp_norm_weight)]

	
	opt = Optimizer(model.input, losses)
	seq = opt.minimize(**optimizer_params)
	
	return seq
	
	
def compile_saliency_function(model):
    """
    Compiles a function to compute the saliency maps and predicted classes
    for a given minibatch of input genomic sequences.
    """

    inp = model.layers[0].input
    outp = model.layers[-2].output
    max_outp = T.max(outp, axis=1)
    print(inp._keras_shape, outp._keras_shape)
    saliency = T.grad(max_outp.sum(), wrt=inp)
    max_class = T.cast(max_outp > 0, 'int16')
    return K.function([inp, K.learning_phase()], [saliency, max_class])
