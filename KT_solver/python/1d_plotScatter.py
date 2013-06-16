from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

# name specified as filename in the C++ code
name = 'burgers'

# loads info file containing dimensions of the finite-volume cells, timestep and current time
# arrInfo = np.loadtxt(name + '_info.tsv',usecols=(1,2,3))
filenumb = 400

# loads data file into a 2D array (arr[i][j] is the element at line i, column j in the tsv file)
for i in xrange(0, filenumb):
	filename = '../Results/burgers_' + str(i) + '.tsv'
	arrData = np.loadtxt(filename)
	x = arrData[:,0]
	y = arrData[:,1]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(x, y, c="tomato", s=20)

plt.show()
