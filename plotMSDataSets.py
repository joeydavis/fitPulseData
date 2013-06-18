# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:05:35 2013

@author: jhdavis
"""

import os
import pprint
import random
import wx
import csv
import string
import qMS
import sets
import numpy
import matplotlib as mpl
import matplotlib.pyplot as pylab
from mpl_toolkits.mplot3d import Axes3D

def randrange(n, vmin, vmax):
    return (vmax-vmin)*numpy.random.rand(n) + vmin

def readFile(datapath):
    data = list(csv.reader( open(datapath, 'rU'), delimiter=' '))
    xs = []
    ys = []
    for line in data:
        xs.append(float(line[0]))
        ys.append(float(line[1]))
    return [xs, ys, datapath]

if __name__=="__main__":
    path = "/home/jhdavis/scripts/python/figures/fitPulseData/rawMSData/"

fileList = ["a.txt", "b.txt", "c.txt", "d.txt", "e.txt", "f.txt"]

pylab.close('all')
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')
low = 39
high = -82
[xs, ys, name] = readFile(path+fileList[2])
xtest = numpy.array(xs[low:high])
ytest = numpy.array(ys[low:high])
zs = [3]*len(xtest)
ztest = numpy.array(zs)
ax.plot(xtest,ztest,ytest, color='#005824', lw=1.5, label="0 mins")

top = float(max(ytest))

[xs, ys, name] = readFile(path+fileList[1])
xtest = numpy.array(xs[low:high])
ytest = numpy.array(ys[low:high])
zs = [2]*len(xtest)
ztest = numpy.array(zs)
ax.plot(xtest,ztest,ytest, color = '#41AE76', lw=1.5, label="x mins")

top = max([top, float(max(ytest))])

[xs, ys, name] = readFile(path+fileList[0])
xtest = numpy.array(xs[low:high])
ytest = numpy.array(ys[low:high])
zs = [1]*len(xtest)
ztest = numpy.array(zs)
ax.plot(xtest,ztest,ytest, color = '#99D8C9', lw=1.5, label="y mins")

top = max([top, float(max(ytest))])

ax.w_xaxis.pane.set_visible(False)
ax.w_yaxis.pane.set_visible(False)
ax.w_zaxis.pane.set_visible(False)

#ax.w_xaxis.gridlines.set_visible(False)
#ax.w_yaxis.gridlines.set_visible(False)
#ax.w_zaxis.gridlines.set_visible(False)

ax.w_xaxis.gridlines.set_linewidth(2)
ax.w_yaxis.gridlines.set_linewidth(2)
ax.w_zaxis.gridlines.set_linewidth(2)

[i.set_linewidth(2) for i in ax.w_xaxis.get_ticklines()]
[i.set_linewidth(2) for i in ax.w_yaxis.get_ticklines()]
[i.set_linewidth(2) for i in ax.w_zaxis.get_ticklines()]

ax.w_xaxis.line.set_linewidth(2.5)
ax.w_yaxis.line.set_linewidth(2.5)
ax.w_zaxis.line.set_linewidth(2.5)
           
ax.set_yticks([1,2,3])
ax.set_yticklabels(['26','52','78'])


ax.set_zticks([0, top/3, 2*top/3, top])
ax.set_zticklabels(['','',''])
ax.set_zlim3d([0, top])
#ax.set_zlabels(['', '', ''])

ax.set_xticks([544, 546, 548, 550, 552])
ax.set_xlim3d([544, 552])

ax.set_xlabel("mass")
ax.set_zlabel("intensity")

ax.view_init(10, -75)

pylab.ylim(0.5, 3.5)
pylab.show()
