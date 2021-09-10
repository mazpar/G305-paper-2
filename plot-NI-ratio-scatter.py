import ast
import numpy as np
from scipy import stats
import math
import matplotlib.pyplot as plt
dirName = '/home/pmazumdar/Documents/LASMA/Reduction/scripts/'

data=np.genfromtxt(dirName+'nco_intensity.dat',dtype=float, usecols=(0,1))

plt.rcParams.update({'font.size': 15})
#plt.scatter(data[:,0],data[:,1])
#bins = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10,20]

median,bin_edge,binnumber = stats.binned_statistic(data[:,0],data[:,1],statistic='median',bins=100)
mean,mean_bin_edge,mean_binnumber = stats.binned_statistic(data[:,0],data[:,1],statistic='mean',bins=100)
plt.rcParams.update({"font.size":9})
plt.xlim(0.5,15)
plt.ylim(4e15,2e17)
plt.yscale('log')
plt.xscale('log')
plt.hlines(median,bin_edge[:-1],bin_edge[1:],colors='g', lw=5, label='median')
plt.hlines(mean,bin_edge[:-1],bin_edge[1:],colors='r', lw=5, label='mean')
plt.ylabel(r'$N/I\,(^{13}CO)\,[mol.cm^{-2}.(K.km.s^{-1})^{-1}]$', fontsize=9)
plt.xlabel(r'$I(^{13}CO) [K.km.s^{-1}]$', fontsize=9)
3plt.scatter(data[:,0],data[:,1],s=0.3,c='cyan')
plt.legend(loc='upper right')
plt.savefig("NI-ratio-scatter.eps",dvi=300,overwrite=True,format='eps')