import numpy as np

def find_index(wg, ps):
  """wg: whole genome
     ps: positive sites
  """
  index = np.argsort(wg)
  sorted_wg = wg[index]
  sorted_index = np.searchsorted(sorted_wg, ps)
  ps_index = np.take(index, sorted_index, mode="clip")
  return ps_index
