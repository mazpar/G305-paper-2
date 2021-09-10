#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 14:16:04 2021

Make RGB Images

@author: Rohit
@modified by: Parichay
"""

import aplpy
import numpy as np
import matplotlib.pyplot as plt

rgb_path = "/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/"
plot_path = "/home/pmazumdar/Documents/LASMA/G305 Papers/Paper 2/"

aplpy.make_rgb_cube([rgb_path+'G305-MSX-21-reproject.fits', \
                     rgb_path+'G305-GLM-8-reproject.fits', \
                     rgb_path+'G305-HPACS-70-reproject.fits'],\
                    rgb_path+'G305_rgb.fits', north=False)

aplpy.make_rgb_image(rgb_path+'G305_rgb.fits', rgb_path+'example_cube.eps',
                     vmin_r=-1e-6, vmax_r=4.5e-5, stretch_r='linear',
                     vmin_g=70, vmax_g=550, stretch_g='sqrt',
                     vmin_b=0.3, vmax_b=2, stretch_b='log')

###############################################################################################

aplpy.make_rgb_cube([rgb_path+'G305-MSX-21-reproject.fits', \
                     rgb_path+'G305-GLM-8-reproject.fits', \
                     rgb_path+'G305_13CO_mom_0.fits'],\
                    rgb_path+'G305_rgb_13CO.fits', north=False)

aplpy.make_rgb_image(rgb_path+'G305_rgb_13CO.fits', rgb_path+'example_cube_13CO.eps',
                     vmin_r=-1e-6, vmax_r=4.5e-5, stretch_r='linear',
                     vmin_g=70, vmax_g=550, stretch_g='sqrt',
                     vmin_b=0.5, vmax_b=45, stretch_b='sqrt')

plt.rcParams.update({'font.size':12})

Frgb = aplpy.FITSFigure(rgb_path+'G305_rgb_13CO_2d.fits')
Frgb.show_rgb(rgb_path+'example_cube_13CO.eps')

x_Danks = np.array([305.3384,305.3934])
y_Danks = np.array([+00.0719,+00.0874])
radii_Danks = np.array([0.018,0.026])
Frgb.show_circles(x_Danks,y_Danks,radius=radii_Danks,edgecolor='white',facecolor='none',linewidths=1)

x_wr48 = 305.361
y_wr48 = +00.056

Frgb.show_markers(x_wr48,y_wr48,edgecolor='white',marker='*', facecolor='none',s=100,linewidth=1)

Frgb.savefig(plot_path+'G305_RGB.eps')
#Frgb.savefig(plot_path+'G305_RGB.png',dpi=400)

###   Test RGB Range: ##########################################################################
  
Fr = aplpy.FITSFigure(rgb_path+'G305-MSX-21-reproject.fits')
Fr.show_colorscale(vmin=0.00001,vmax=0.000175, stretch='log')
Fr.add_colorbar()

Fg = aplpy.FITSFigure(rgb_path+'G305-GLM-8-reproject.fits')
Fg.show_colorscale(vmin=70,vmax=600,stretch='log')
Fg.add_colorbar()

Fb = aplpy.FITSFigure(rgb_path+'G305-HPACS-70-reproject.fits')
Fb.show_colorscale(vmin=0.5,vmax=7,stretch='log')
Fb.add_colorbar()


