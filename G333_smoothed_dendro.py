#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:26:10 2021

G333 Smoothing and Dendrogram

@author: pmazumdar
"""

import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.convolution import Gaussian1DKernel
import matplotlib.pyplot as plt
from astropy.wcs import WCS


fits_path = "/home/pmazumdar/Documents/LASMA/Ancillary_Data/"

## get the cube
cube = SpectralCube.read(fits_path+"SEDIGISM_G333_13CO21.fits")

## change velocity units to km/s
cube2 = cube.with_spectral_unit(u.km/u.s)

## extract slice from -70 to -20 km/s
subcube = cube2.spectral_slab(-70*u.km/u.s,-20*u.km/u.s)

## Extract Spatial Subcube
subcube_final = subcube.subcube(xlo = 332.9*u.deg, xhi = 333.9*u.deg,
                                ylo = -0.5*u.deg, yhi = 0.3*u.deg,
                                zlo = 'min', zhi = 'max')


## Spectral smoothing

fwhm_factor = np.sqrt(8*np.log(2))
current_resolution = 0.25 * u.km/u.s
target_resolution = 0.5 * u.km/u.s
pixel_scale = 0.25 * u.km/u.s
gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                  pixel_scale / fwhm_factor)
kernel = Gaussian1DKernel(gaussian_width.value)
new_cube = subcube_final.spectral_smooth(kernel)
new_cube.write('smoothed_G333.fits',overwrite=True)

## Run Dendrogram

#load_newcube
new_cube = SpectralCube.read(fits_path+"smoothed_G333.fits")


from astropy.io import fits
from astrodendro.pruning import all_true, min_vchan, min_delta, min_area
from astropy import constants as const
import aplpy
import seaborn as sns
import scipy.stats as sst
import radio_beam
from astrodendro import Dendrogram, ppv_catalog, structure
from astrodendro.pruning import all_true, min_vchan, min_delta, min_area
from astropy import constants as const
import aplpy
import seaborn as sns
import scipy.stats as sst


data = new_cube.hdu.data
hd = new_cube.hdu.header
wc = WCS(hd)

##    Custom Definitions for the Dendrogram    ##
rms = 0.45 # rms noise
bmaj = hd['bmaj'] # beam_major
bmin = hd['bmin'] # beam_minor
cdelt1 = hd['cdelt1'] # delta_x
cdelt2 = hd['cdelt2'] # delta_y
deltav_kms = abs(hd['CDELT3']) # vel res in kmps
ppb = abs((bmaj*bmin)/(cdelt1*cdelt2)*2*np.pi/(8*np.log(2))) # pixel_per_beam


is_independent = all_true((min_delta(5*rms), min_area(1*ppb), min_vchan(12)))

dG333 = Dendrogram.compute(data, min_value=5*rms, wcs=wc, is_independent = is_independent, verbose=1)

dG333.save_to('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_13CO_smoothed_dendro.fits')

#make catalog

metadata = {}
metadata['data_unit'] = u.Jy
metadata['beam_major'] = (bmaj * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['beam_minor'] = (bmin * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['velocity_scale'] = 0.25 * u.km/u.s # v_res

cat = ppv_catalog(dG333,metadata)



dists = 3600*u.parsec # Distance of the source.
x2 = 5.4 # x2 = X_CO / 2E20 [cm^-2 / K km s^-1] , where X_CO is a CO-H2 conversion factor.
#x2 = 1 # for 12CO
sigma_v = np.array(cat['v_rms'])
sigma_x = np.array(cat['radius'])*(((1*u.arcsecond).to(u.rad)).value)*dists.value
eta = 1.91 # conversion factor. R = eta * sigma_r
G = 4.302e-3 # units of pc (M_sun)^-1 (km/s)^2
deltax_pc = abs(np.pi/180.*hd['CDELT1']*dists.value) # delta x in pc
deltay_pc = abs(np.pi/180.*hd['CDELT2']*dists.value) # delta y in pc
sigma_majs = cat['major_sigma']
sigma_mins = cat['minor_sigma']
mu = 2.72 # Average H2 mass including He fraction
mp = 8.4089382e-58*u.solMass # Proton Mass
nu_12CO = 345.79598990*u.GHz
nu_13CO = 330.58796530*u.GHz
delta_nu_12 = 0.0011534512649414282*u.GHz
delta_nu_13 = 0.0011027227552869259*u.GHz

##
##                                Additions to Integrated Intensity Catalog
##


##    adding a radius column to the catalog
cat['radius_pc'] = np.zeros(len(cat),dtype=float)
cat['radius_pc'] = eta*np.sqrt((sigma_majs*deltax_pc)*(sigma_mins*deltay_pc))
cat['radius_pc'].unit = u.parsec


##    adding a luminosity column to the catalog
cat['luminosity']=np.zeros(len(cat),dtype=float)
cat['luminosity'] = cat['flux']*deltav_kms*deltax_pc*deltay_pc
cat['luminosity'].unit = u.K * u.km / u.s * u.pc * u.pc


##    adding a mass column to the catalog
cat['Mass']=np.zeros(len(cat),dtype=float)
cat['Mass'] = cat['luminosity']*4.4*x2
cat['Mass'].unit = u.solMass


##    adding a surface density column to the catalog
cat['Sigma_exact']=np.zeros(len(cat),dtype=float)
cat['Sigma_exact'] = cat['Mass']/(cat['area_exact']*deltax_pc*deltay_pc)
cat['Sigma_exact'].unit = u.solMass/(u.pc*u.pc)

cat['Sigma_ellipse'] = np.zeros(len(cat),dtype=float)
cat['Sigma_ellipse'] = cat['Mass']/(np.pi*cat['radius_pc']**2)
cat['Sigma_ellipse'].unit = u.solMass/(u.pc*u.pc)

##    calculating virial parameter alpha
cat['virial_parameter'] = np.zeros(len(cat),dtype=float)
cat['virial_parameter'] = (5*((sigma_v)**2)*cat['radius_pc'])/(4.4*x2*cat['luminosity']*G)
cat['virial_parameter'].unit = ''



##    Read the 8um Map

hdu_g333_glm8 = fits.open('/home/pmazumdar/Documents/LASMA/Ancillary_Data/GLIMPSE/GLM_33300+0000_mosaic_I4.fits')[0]
data_g333_glm8 = hdu_g333_glm8.data
w_g333_glm8 = WCS(hdu_g333_glm8.header)


# Create 2D cutout
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

position = SkyCoord('333.4deg -0.1deg', frame='galactic')
size = u.Quantity((0.8, 1.0), u.deg) # ny,nx order
cutout = Cutout2D(data_g333_glm8, position, size, wcs=w_g333_glm8)

# Put the cutout image in the FITS HDU
hdu_g333_glm8.data = cutout.data

# Update the FITS header with the cutout WCS
hdu_g333_glm8.header.update(cutout.wcs.to_header())

hdu_g333_glm8.writeto('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_8um_cutout.fits', overwrite=True)

##    Reproject the cube on 8um Map

from reproject import reproject_interp

hdu_T = new_cube[0,:,:].hdu
array, footprint = reproject_interp(hdu_g333_glm8, hdu_T.header)
fits.writeto('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_8um_reproject.fits', array, hdu_T.header, clobber=True)

g333_hdu_dust = fits.open('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_8um_reproject.fits')[0]
g333_w_dust = WCS(g333_hdu_dust.header)




#---- Define Threshold values and make mask

threshold = 200 #lower threshold
g333_glm8_mask = g333_hdu_dust.data > threshold # 8um mask




##   Adding mean 8um Flux
cat['8um Flux Mean']=np.zeros(len(cat),dtype=float)

for l in dG333.all_structures:    ## adding average 8um flux of leaves
  leaf_mask_3d = l.get_mask()
  leaf_mask = np.max(leaf_mask_3d,axis=0)
  masked_8um_map = g333_hdu_dust.data[leaf_mask]
  mean_8um_flux = np.nanmean(masked_8um_map)
  cat[l.idx]['8um Flux Mean']=mean_8um_flux


### Mask Tags to Catalog

cat['Mask_Tag']=np.zeros(len(cat),dtype='str')

cat['Structure Tag'] = np.zeros(len(cat),dtype='str')


for l in dG333.all_structures:

  sr_mask_3d = l.get_mask()   # Make a structure mask
  sr_mask = np.max(sr_mask_3d,axis=0)  # Project it to 2D
  and_mask = sr_mask & g333_glm8_mask  # An overlap mask of the two masks

  # 8um Tag
  if 100.*(np.sum(and_mask)/np.sum(sr_mask))>67:  # > 2/3rd area I the mask
    cat[l.idx]['Mask_Tag']= "I"
  elif 100.*(np.sum(and_mask)/np.sum(sr_mask))<=10:  # > 1/10 area I the mask
    cat[l.idx]['Mask_Tag']= "O"
  else:
    cat[l.idx]['Mask_Tag']= "P"

  # Leaf or Branch or Trunk
  if l.is_leaf:
    cat[l.idx]['Structure Tag']="l"
  elif l.is_branch:
    cat[l.idx]['Structure Tag']="b"
  else:
    cat[l.idx]['Structure Tag']="t"

cat['Trunk Tag'] = np.zeros(len(cat),dtype=int)
cat['Trunk Tag'] = 0
for t in dG333.trunk:
  cat[t.idx]['Trunk Tag']=1  



# save catalog
cat.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_13CO_smoothed_cat.fits',overwrite=True)



#
#  PLOTTING
# 

################################
# Plotting Scaling Relations
################################

L_CO = cat['luminosity'].astype(np.float64)
radius_pc = cat['radius_pc'].astype(np.float64)
vrms = sigma_v.astype(np.float64)
mass = cat['Mass'].astype(np.float64)
Sigma_exact = cat['Sigma_exact'].astype(np.float64)
Sigma_ellipse = cat['Sigma_ellipse'].astype(np.float64)
alpha = cat['virial_parameter'].astype(np.float64)
mask_tag = cat['Mask_Tag']
sr_type = cat['Structure Tag']
trunk_tag = cat['Trunk Tag']
y_var = vrms/np.sqrt(radius_pc)
glm8_flux =np.array(cat['8um Flux Mean'])



#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&%    Heyer et al. Plot
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%



#         Separate variables for leaves, trunks and trunks
# --------------------------------------------------------------------------------------------------
mask_tag_leaves = mask_tag[(sr_type=="l")]
mass_leaves = mass[(sr_type=="l")]
sigma_leaves = Sigma_ellipse[(sr_type=="l")]
y_var_leaves = y_var[(sr_type=="l")]

mask_tag_branches = mask_tag[(sr_type=="b")&(trunk_tag==0)]
mass_branches = mass[(sr_type=="b")&(trunk_tag==0)]
sigma_branches = Sigma_ellipse[(sr_type=="b")&(trunk_tag==0)]
y_var_branches = y_var[(sr_type=="b")&(trunk_tag==0)]

mask_tag_trunks = mask_tag[(sr_type=="b")&(trunk_tag==1)]
mass_trunks = mass[(sr_type=="b")&(trunk_tag==1)]
sigma_trunks = Sigma_ellipse[(sr_type=="b")&(trunk_tag==1)]
y_var_trunks = y_var[(sr_type=="b")&(trunk_tag==1)]
  

plt.rcParams.update({'font.size': 14})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,7))

#         kde plots for the leaves
#

sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="I")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="I")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='firebrick',\
            fill=0,\
            alpha=1,\
            label='Mostly Inside')
            
sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="P")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="P")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='tab:green',\
            fill=0,\
            alpha=1,\
            label='Partly Inside')
            
sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="O")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="O")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='tab:blue',\
            fill=0,\
            alpha=1,\
            label='Outside')
            

#         Scatter plot of branches based on I and O the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(sigma_branches[(mask_tag_branches=="I")], y_var_branches[(mask_tag_branches=="I")],\
                         color='firebrick',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_branches[(mask_tag_branches=="P")], y_var_branches[(mask_tag_branches=="P")],\
                         color='tab:green',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_branches[(mask_tag_branches=="O")], y_var_branches[(mask_tag_branches=="O")],\
                         color='tab:blue',marker='o',s=10,alpha=0.7)        # for column density dendrogram


#         Scatter plot of trunks based on I and O the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(sigma_trunks[(mask_tag_trunks=="I")], y_var_trunks[(mask_tag_trunks=="I")],\
                         color='firebrick',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_trunks[(mask_tag_trunks=="P")], y_var_trunks[(mask_tag_trunks=="P")],\
                         color='tab:green',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_trunks[(mask_tag_trunks=="O")], y_var_trunks[(mask_tag_trunks=="O")],\
                         color='tab:blue',marker='*',s=50,alpha=0.7)        # for column density dendrogram


x2data = np.logspace(0,4.5,1000)
y2data = np.logspace(-1,1,1000)

X2, Y2 = np.meshgrid(x2data, y2data)

def balance (x,y,press=0,vir=1):   # pressure in units of K/cm^3
    conv_factor = 0.0020399266017097576 # convert to consistent unit
    press_conv = press*conv_factor
    return y**2-((4*1.9/3)*press_conv/x)-(1.9*np.pi*G*x*vir/5)

vir1=plt.contour(X2, Y2, balance(X2,Y2,press=0,vir=1), levels=[0], colors=['#014182'])
fmt = {}
strs = [r'$\alpha_{vir}=1$']
for l, s in zip(vir1.levels, strs):
    fmt[l] = s
manual_location=[(10,0.19)]
plt.clabel(vir1, fmt=fmt, inline=True, fontsize=11, manual=manual_location)

vir2=plt.contour(X2, Y2, balance(X2,Y2,press=0,vir=2), fmt='-.', levels=[0],colors=['#014182'])
fmt = {}
strs = [r'$\alpha_{vir}=2$']
for l, s in zip(vir2.levels, strs):
    fmt[l] = s
manual_location=[(6,0.3)]
plt.clabel(vir2, fmt=fmt, inline=True, fontsize=11, manual=manual_location)

con1=plt.contour(X2, Y2, balance(X2,Y2,press=1e3), levels=[0], colors=['black'])
fmt = {}
strs = [r"$P=10^3 \, \rm{K}\cdot\rm{cm}^{-3}$"]
for l, s in zip(con1.levels, strs):
    fmt[l] = s
manual_location=[(3,1.5)]
plt.clabel(con1, fmt=fmt, inline=True, fontsize=11, manual=manual_location)

con2=plt.contour(X2, Y2, balance(X2,Y2,press=1e4), levels=[0], colors=['black'])
fmt = {}
strs = [r"$P=10^4 \, \rm{K}\cdot\rm{cm}^{-3}$"]
for l, s in zip(con2.levels, strs):
    fmt[l] = s
manual_location=[(3.5,3)]
plt.clabel(con2,fmt=fmt,inline=1,fontsize=11, manual=manual_location)

con3=plt.contour(X2, Y2, balance(X2,Y2,press=1e5), levels=[0], colors=['black'])
fmt = {}
strs = [r"$P=10^5 \, \rm{K}\cdot\rm{cm}^{-3}$"]
for l, s in zip(con3.levels, strs):
    fmt[l] = s
manual_location=[(20,7)]
plt.clabel(con3, fmt=fmt, inline=True, fontsize=11, manual=manual_location)

plt.ylim(bottom=0.1,top=10)
plt.xlim(left=1,right=10**4.4)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\sigma_v / R^{0.5} \,\,\, [km.s^{-1}.pc^{-0.5}]$')
plt.xlabel(r'$\Sigma \,\,\, [M_{\odot}.pc^{-2}]$')
#plt.legend(loc='lower right',ncol=2,framealpha=0.5)

labels = ["Branches","Trunks"]
mtype = ["o","*"]
for i, stype in enumerate(labels):
    ax.scatter([], [], marker=mtype[i], c='k',label=str(labels[i]))
ax.legend(scatterpoints=1, frameon=False, labelspacing=0.8, loc='lower right',fontsize=13)

plt.show()
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G333_Dynamical_State_vs_8um.png", format="png",dpi=800)

