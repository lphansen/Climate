#%%
import numpy as np
import pickle
import scipy
#%%
csr_mat = pickle.load(open("csr_mat", "rb"))
csr_mat
#%%
b = np.load("b.npy")
b.shape

#%%
v_cpp = np.load("v_cpp.npy")
v_cpp.shape

#%%
res = np.max(np.abs(csr_mat.multiply(v_cpp) - b))
res