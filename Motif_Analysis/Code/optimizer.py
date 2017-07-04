#*********************
#   Fei Tan
#   ft54@njit.edu
#   March 28, 2017
#********************


import numpy as np
import pprint

from keras import backend as K
from collections import OrderedDict

class Optimizer(object):
	
	def __init__(self, seq_input, losses, wrt=None):

		self.seq = seq_input
		self.loss_functions = []
		self.wrt = self.seq if wrt is None else wrt

		overall_loss = K.variable(0.)

		for loss, weight in losses:
			if weight != 0:
				loss_fn = weight * loss.build_loss()
				overall_loss += loss_fn
				# learning phase: 0 (test), 1(train)
				self.loss_functions.append( ( loss.name, K.function( [self.seq, K.learning_phase()], [loss_fn]) ) )

		grads = K.gradients(overall_loss, self.wrt)[0]
		grads = grads / (K.sqrt(K.mean(K.square(grads))) + K.epsilon())

		self.overall_loss_grad_wrt_fn = K.function([self.seq, K.learning_phase()], [overall_loss, grads, self.wrt])
		
	
	def eval_losses(self, seq):

		losses = OrderedDict()
		for name, fn in self.loss_functions:
			losses[name] = fn([seq, 0])
		return losses

	def rmsprop(self, grads, cache=None, decay_rate=0.95):
		if cache is None:
			cache = np.zeros_like(grads)
		cache = decay_rate * cache + (1 - decay_rate) * (grads ** 2)
		step = -grads / np.sqrt(cache + K.epsilon())
		return step, cache

	def get_seed_seq(self, seed_seq):
		sample_size, filter_length, nucleotide_size = self.seq._keras_shape
		print(filter_length, nucleotide_size)
		if seed_seq is None:
			seed_seq = np.ones((filter_length, nucleotide_size))/4.

		seed_seq = np.array([seed_seq], dtype=np.float32)

		return seed_seq


	def minimize(self, seed_seq=None, max_iter=200, lr = 0.1, verbose=True):
		"""Performs gradient descent on the input sequenc  with respect to defined losses"""
		
		seed_seq = self.get_seed_seq(seed_seq)
		cache = None
		best_loss = float('inf')
		best_seq = None

		grads = None

		
		for i in range(max_iter):
			overall_loss, grads, wrt_value = self.overall_loss_grad_wrt_fn([seed_seq, 0])

			if verbose:
				losses = self.eval_losses(seed_seq)
				print('Interation: {}, losses: {}, overall loss: {}'.format(i+1, pprint.pformat(losses), overall_loss))
	
			
			#gradient descent update
			if self.wrt is self.seq:
				step, cache = self.rmsprop(grads, cache)
				
				#to be revised later
				seed_seq += lr * step
			
			if overall_loss < best_loss:
				best_loss = overall_loss
				best_seq = seed_seq.copy()


		return best_seq[0], grads, wrt_value
