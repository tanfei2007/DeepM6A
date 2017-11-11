import numpy as np


x = np.array([3,5,7,1,9,8,6,6])
y = np.array([2,1,5,10,100,6])


def find_index(wg, ps):
  """wg: whole genome
     ps: positive sites
  """
  index = np.argsort(wg)
  sorted_wg = wg[index]
  sorted_index = np.searchsorted(sorted_wg, ps)

  ps_index = np.take(index, sorted_index, mode="clip")
  
  return ps_index
