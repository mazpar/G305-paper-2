#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 12:09:31 2021

Draw histogram of noise distribution in fits map

@author: pmazumdar
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import seaborn as sns
import scipy.stats as sst

#filepath = "/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/"
filepath = "/home/pmazumdar/Documents/LASMA/Reduction/datacubes/G346/lmv_maps/"

#Load the FITS files

#noise_hdu_12 = fits.open(filepath+"G305_12CO_noise.fits")[0]
#noise_hdu_13 = fits.open(filepath+"G305_13CO_noise.fits")[0]
noise_hdu_12 = fits.open(filepath+"G346_12CO_noise.fits")[0]
noise_hdu_13 = fits.open(filepath+"G346_13CO_noise.fits")[0]
noise_hdu_SED = fits.open(filepath+"G346_SEDIGISM_13CO21_noise.fits")[0]


#Flatten the dataset

noise_12 = noise_hdu_12.data.flatten()
noise_13 = noise_hdu_13.data.flatten()
noise_SED = noise_hdu_SED.data.flatten()

#only for G346
noise_12[noise_12>0.6] = np.nan
noise_13[noise_13>0.6] = np.nan

#Plot the histogram
plt.rcParams.update({'font.size':14})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))

sns.histplot(noise_12,\
             ax=ax,\
             fill=0,\
             alpha=0.6,\
             log_scale=10,\
             #binwidth=0.05,\
             element = 'step',\
             color = "tab:blue",\
             label = r'$^{12}$CO (3-2)',\
             legend=False)
sns.histplot(noise_13,\
             ax=ax,\
             fill=0,\
             alpha=0.6,\
             #binwidth=0.05,\
             log_scale=10,\
             element = 'step',\
             color = "firebrick",\
             label = r'$^{13}$CO (3-2)',\
             legend=False)
sns.histplot(noise_SED,\
             ax=ax,\
             fill=0,\
             alpha=0.6,\
             #binwidth=0.05,\
             log_scale=10,\
             element = 'step',\
             color = "tab:green",\
             label = r'$^{13}$CO (2-1)',\
             legend=False)

#G305
#ax.axvline(0.13,linewidth=1.2,linestyle='dotted',color='tab:blue') 
#ax.axvline(0.29,linewidth=1.2,linestyle='dotted',color='firebrick') 

#G346
ax.axvline(0.19,linewidth=1.2,linestyle='dotted',color='tab:blue') 
ax.text(0.19,100,'0.19 K',color='tab:blue',rotation=90)
ax.axvline(0.22,linewidth=1.2,linestyle='dotted',color='firebrick') 
ax.text(0.22,100,'0.22 K',color='firebrick',rotation=90)
ax.axvline(0.58,linewidth=1.2,linestyle='dotted',color='tab:green') 
ax.text(0.58,100,'0.58 K',color='tab:green',rotation=90)
#ax[0].set_xlim(0.01,3)
#ax[1].set_xlim(0.01,3)
#ax.set_xlim(0.02,3)

#ax[0].set_yscale('log')
#ax[1].set_yscale('log')

#ax[1].set_ylabel(r'Number of pixels')
#ax[0].set_ylabel('')
#ax[0].set_xlabel(r'RMS Noise $^{12}$CO [K]')
#ax[1].set_xlabel(r'RMS Noise $^{13}$CO [K]')

ax.set_xlabel(r'RMS Noise [K]')
ax.set_ylabel(r'Number of pixels')

plt.legend(loc='upper right',frameon=False, labelspacing=0.2, fontsize=12)
plt.tight_layout()

#plt.savefig(filepath+"plots/noise_histogram.pdf")
plt.savefig(filepath+"noise_histogram_G346.pdf")
