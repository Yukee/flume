#!/usr/bin/env python
"""
usage: python plotContour.py filename filenumber
"""
from os import system
import pylab
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from sys import argv

# get the path of the file which will be loaded

filename = argv[1]
number = int(argv[2])

loadpath = "../Results/" + filename + "_" + str(number) + ".tsv"
infopath = "../Results/" + filename + "_infos.tsv"

# unpack = True option transposes the returned array. dtype option must always be used, otherwise python will search for datatypes across the whole file, which can be time consuming.

x, y, z = np.genfromtxt(loadpath, unpack = True, dtype = 'f8')

# interpolate z(x,y) at points xi, yi using Delaunay triangulation aka natural neighbour interpolation ('nn' option)
#z = mlab.griddata(x, y, z, x, y, interp = 'nn')

z = z.reshape(20,400)
x = np.linspace(-20,-0.05,400)
y = np.linspace(0,0.95,20)

# array containing value of levels
levels = np.linspace(0.2, 1, 100)

fig = plt.figure()
CS1 = plt.contourf(x, y, z, levels, cmap = plt.cm.coolwarm, vmax = np.max(levels), vmin = np.min(levels))
fig.colorbar(CS1)

plt.show()
