from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

# loads data file into a 2D array (arr[i][j] is the element at line i, column j in the tsv file)
arrData = np.loadtxt("../Results/caca_4.tsv")
x = arrData[:,0]
y = arrData[:,1]
z = arrData[:,2]

# loads info file containing dimensions of the finite-volume cells, timestep and current time
# arrInfo = np.loadtxt("info.tsv",usecols=(1,2,3))

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ax.scatter(x,y,z)

plt.show()
