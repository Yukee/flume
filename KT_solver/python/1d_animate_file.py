"""
  Module for testing population sim over time for a universal scale.

  Use matplotlib to display population changes.
"""
#!/usr/bin/env python

import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import numpy as np
import time

class ScatterAnim(object):
    def __init__(self, name, filenumb, xinterval, lowerleftcorner, sleep=False):
	self.filename = '../Results/' + name + '_'
	self.framesnumber = filenumb
	self.xinterval = xinterval
	self.llc = lowerleftcorner

        self.tol = 0.1  #0.05
        self.x = 1.0
        self.y = 1.0
        self.x1 = np.loadtxt(self.filename + str(0) + '.tsv')[:,0]

        self.xwork = np.copy(self.x1)
        self.ywork = np.copy(self.y1(0))

	self.ymax = max(self.ywork)
	self.ymin = min(self.ywork)

        self.fig = pyplot.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()

        self.ani = animation.FuncAnimation(self.fig, self.update, interval=20,
            init_func=self.setup_plot, blit=True,
            repeat=False, frames=self.framesnumber) # limit frames

        self.first = sleep

    def y1(self, i):
	return np.loadtxt(self.filename + str(i) + '.tsv')[:,1]

    def setup_plot(self):
        self.ax.axis([-0.1 + self.llc, self.llc + self.xinterval + 0.1, 1.1*self.ymin, 1.1*self.ymax])
        self.scat = self.ax.scatter(self.xwork, self.ywork, c="tomato", s=20, animated=True)

        return self.scat,

    def update(self, i):
        if self.first:
            time.sleep(10)
            self.first = False

        self.ywork = self.analog_filter(self.ywork, self.y1(i), 1) # setting factor to 1 deactivates the filter
        self.scat = self.ax.scatter(self.xwork, self.ywork, c="tomato", s=20, animated=True)

        return self.scat,

    def show(self):
        try:
            pyplot.show()
        except AttributeError as e:
            # eat error from lib
            pass
            print "Exception {0}".format(e)

    def analog_filter(self, old, new, factor=0.1):
        old = old * (1.0 - factor) + new * factor
        return old

def main():

    # not plotting anything in movie
    # Don't try this, it is busted and does not work.  I used CamStudio to make the video on youtube.
    #a.ani.save('scatter_transform.mp4', fps=30, codec='libx264')

    # name specified as filename in the C++ code
    name = 'burgers'

    # loads info file containing dimensions of the finite-volume cells, timestep and current time
    # arrInfo = np.loadtxt(name + '_infos.tsv',usecols=(1,2,3))
    filenumb = 1000
    xinterval = 2 * np.pi
    lowerleftcorner = 0

    # use this to delay for making a video (gives you time to manually start the recording)
    #a = ScatterAnim(name, filenumb, sleep=True)
    # use this for normal viewing
    a = ScatterAnim(name, filenumb, xinterval, lowerleftcorner)

    a.show()

if __name__ == '__main__':
    main()

