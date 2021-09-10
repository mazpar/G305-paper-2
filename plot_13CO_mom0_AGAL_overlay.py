#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 11:57:12 2021

@author: pmazumdar
"""

import matplotlib.pyplot as plt
import aplpy

path = "/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/"

plt.rcParams.update({'font.size': 10})
fig= plt.figure(figsize=(6,5))

F = aplpy.FITSFigure("G305_13CO_mom_0.fits",figure=fig)
F.show_colorscale(cmap='Blues',stretch='log')
F.set_nan_color('white')
F.add_colorbar()
F.colorbar.set_ticks([3,10,30,70])
F.colorbar.set_width(0.1)
F.colorbar.set_axis_label_pad(0.1)
F.colorbar.set_axis_label_text("K. km. s$^{-1}$")
F.show_contour("G305-agal-planck.fits",\
               levels=[0.67,1,3],\
               smooth=1,\
               linewidths=0.5)
F.add_beam()
F.beam.set_facecolor('k')
F.beam.set_corner('bottom right')

F.save(path+"plots/G305_13CO_mom_0_AGAL_overlay.eps")

