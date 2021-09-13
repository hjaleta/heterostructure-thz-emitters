from math import ldexp
import numpy as np
# import matplotlib.pyplot as plt

# arr1 = np.loadtxt("data/FePt_bilayer-open/output/ntotal_up.out")
# arr2 = np.loadtxt("data/FePt_bilayer-open/output/ntotal_down.out")

# arr_tot = arr1 - arr2

# fig = plt.plot(arr_tot)

# print(len(fig))
# plt.savefig("arr_tot.jpg")
# plt.close()

# arr_flux = np.loadtxt("data/FePt_bilayer-open/flux.out")
# one_layer = arr_flux[:,4]
# # print(arr_flux[0:5,:])
# # print(np.gradient(arr_flux)[0][0:5,:])
# plt.plot(one_layer)
# plt.savefig("one_layer_flux.jpg")

# grad = np.gradient(arr_flux, axis = 0)
# one_grad = grad[:,4]
# plt.close()
# plt.plot(one_grad)
# plt.savefig("One_grad.jpg")

from math import sqrt, floor

print(type(floor(-0.1)))