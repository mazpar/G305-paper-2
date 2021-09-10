#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 13:28:03 2021

@author: pmazumdar
"""

##############################################
##   Smooth 13CO map to CHIMPS Resolution   ##
##############################################
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.table import Table
from astrodendro import Dendrogram, ppv_catalog, structure
from astropy import wcs
from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astrodendro.pruning import all_true, min_vchan, min_delta, min_area
from astropy import constants as const
import aplpy
import seaborn as sns
import scipy.stats as sst
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve



cube = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO_resample.fits')
cube_cd = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal.fits')
beam = radio_beam.Beam(major=27.4*u.arcsec, minor=27.4*u.arcsec, pa=0*u.deg)
new_cube = cube.convolve_to(beam)
new_cube_cd = cube_cd.convolve_to(beam)

hdu_13CO = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO-moment0.fits')[0]
mask_nan = ~np.isnan(hdu_13CO.data) # only include non nan pixels

masked_cube = new_cube.with_mask(mask_nan)   # apply mask to spectral cube
masked_cube_cd = new_cube_cd.with_mask(mask_nan)   # apply mask to spectral cube

data = masked_cube.hdu.data
hd = masked_cube.hdu.header
wc = wcs.WCS(hd)
data_cd = masked_cube_cd.hdu.data
hd_cd = masked_cube_cd.hdu.header
wc_cd = wcs.WCS(hd_cd)


##    Custom Definitions for the Dendrogram    ##
rms = 0.15 # rms noise
rms_cd = 1.6e7
cd_min = 3.37e11
bmaj = hd['bmaj'] # beam_major
bmin = hd['bmin'] # beam_minor
cdelt1 = hd['cdelt1'] # delta_x
cdelt2 = hd['cdelt2'] # delta_y
deltav_kms = abs(hd['CDELT3']/1000.) # vel res in kmps
ppb = abs((bmaj*bmin)/(cdelt1*cdelt2)*2*np.pi/(8*np.log(2))) # pixel_per_beam

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Creating the Dendrogram
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


# Prune leaves below a given height:

#def custom_independent(structure, index=None, value=None):
#    peak_index, peak_value = structure.get_peak()
#    return peak_value > 5

is_independent = all_true((min_delta(5*rms), min_area(1*ppb), min_vchan(6)))
#is_independent_cd = all_true((min_delta(3*cd_min), min_area(1*ppb), min_vchan(2)))

d = Dendrogram.compute(data, min_value=5*rms, wcs=wc, is_independent = is_independent, verbose=1)
d.save_to('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_smoothed_dendro.fits')

#d_cd = Dendrogram.compute(data_cd, min_value=5*cd_min, wcs=wc_cd, is_independent = is_independent_cd, verbose=1)
#d_cd.save_to('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_column_densty_smoothed_dendro.fits')


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% START HERE IF DENDROGRAM ALREADY RUN
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#d = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_dendro.fits')
#d_cd = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_column_densty_dendro.fits')

leaf_id=np.array([])
for s in d.all_structures:
    if s.is_leaf:
        leaf_id=np.append(leaf_id,s.idx)
print(leaf_id)
leaf_id.sort()

leaf_cd_id=np.array([])
#for s in d_cd.all_structures:
#    if s.is_leaf:
#        leaf_cd_id=np.append(leaf_cd_id,s.idx)
#print(leaf_cd_id)
#leaf_cd_id.sort()

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Viewing and Plotting the Dendrogram
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

p = d.plotter()

#p_cd = d_cd.plotter()

#fig,ax = plt.subplots(nrows = 2, ncols=1)
#plt.rcParams.update({"font.size":6})

#ax[0] = fig.add_subplot(1, 1, 1)
#ax[0].set_ylabel('$^{13}$CO Peak Intensity')
#p.plot_tree(ax[0],color='seagreen',lw=0.5)
#p_cd.plot_tree(ax[1],color='orange',lw=0.5)
#ax[1].set_yscale('log')
#ax[1].set_xlabel('Index of Structure')
#ax[1].set_ylabel('$^{13}$CO Column Density')
#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Dendrogram_G305.eps",dpi=300)
#plt.close()


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Creating the ppv_catalog
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

metadata = {}
metadata['data_unit'] = u.Jy
metadata['beam_major'] = (bmaj * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['beam_minor'] = (bmin * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['velocity_scale'] = 0.5 * u.km/u.s # v_res

cat = ppv_catalog(d,metadata)
#cat_cd = ppv_catalog(d_cd,metadata)

# Note : Catalog of Column Density : Flux per pixel represents Col. Density per delX.delY.delV ##
#


dists = 3800*u.parsec # Distance of the source.
dist_gc = 6598.5452942296305*u.parsec
x2 = 18 # x2 = X_CO / 2E20 [cm^-2 / K km s^-1] , where X_CO is a CO-H2 conversion factor.
#x2 = 1 # for 12CO
sigma_v = np.array(cat['v_rms'])
sigma_x = np.array(cat['radius'])*(((1*u.arcsecond).to(u.rad)).value)*dists.value
#sigma_v_cd = np.array(cat_cd['v_rms'])
#sigma_x_cd = np.array(cat_cd['radius'])*(((1*u.arcsecond).to(u.rad)).value)*dists.value
eta = 1.91 # conversion factor. R = eta * sigma_r
G = 4.302e-3 # units of pc (M_sun)^-1 (km/s)^2
deltax_pc = abs(np.pi/180.*hd['CDELT1']*dists.value) # delta x in pc
deltay_pc = abs(np.pi/180.*hd['CDELT2']*dists.value) # delta y in pc
sigma_majs = cat['major_sigma']
sigma_mins = cat['minor_sigma']
#sigma_majs_cd = cat_cd['major_sigma']
#sigma_mins_cd = cat_cd['minor_sigma']
R12_13 = (6.21*(dist_gc.value/1000.0))+18.71 # C12/C13 abundance ratio (Milam et al. 2005)
R12 = 8.5e-5 # C12/H2 abundance (Frerking et al.)
R13_inv =  R12_13/R12
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
#cat['M_lum']=np.zeros(len(cat),dtype=float)
#cat['M_lum'] = cat['luminosity']*4.4*x2
#cat['M_lum'].unit = u.solMass

cat['Mass']=np.zeros(len(cat),dtype=float)
cat['Mass'].unit = u.solMass
for s in d.all_structures:
  str_mask = s.get_mask()
  str_cube = data_cd[str_mask]
  total_cd = np.nansum(str_cube)
  mass = mu*mp.value*R13_inv*total_cd*deltax_pc*deltay_pc*(u.parsec.to(u.centimeter))**2
  cat[s.idx]['Mass'] = mass


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


##
##                                Additions to Integrated Col. Den. Catalog
##

##    adding a radius column to the catalog
#cat_cd['radius_pc'] = np.zeros(len(cat_cd),dtype=float)
#cat_cd['radius_pc'] = eta*np.sqrt((sigma_majs_cd*deltax_pc)*(sigma_mins_cd*deltay_pc))
#cat_cd['radius_pc'].unit = u.parsec


##    adding a mass column to the catalog
#cat_cd['Mass'] = np.zeros(len(cat_cd),dtype=float)
#cat_cd['Mass'] = mu*mp*R13_inv*cat_cd['flux']*deltav_kms*deltax_pc*deltay_pc*(u.parsec.to(u.centimeter))**2
#cat_cd['Mass'].unit = u.solMass


##    adding a surface density column to the catalog
#cat_cd['Sigma_exact']=np.zeros(len(cat_cd),dtype=float)
#cat_cd['Sigma_exact'] = cat_cd['Mass']/(cat_cd['area_exact']*deltax_pc*deltay_pc)
#cat_cd['Sigma_exact'].unit = u.solMass/(u.pc*u.pc)

#cat_cd['Sigma_ellipse'] = np.zeros(len(cat_cd),dtype=float)
#cat_cd['Sigma_ellipse'] = cat_cd['Mass']/(np.pi*cat_cd['radius_pc']**2)
#cat_cd['Sigma_ellipse'].unit = u.solMass/(u.pc*u.pc)


##    calculating virial parameter alpha
#cat_cd['virial_parameter'] = np.zeros(len(cat_cd),dtype=float)
#cat_cd['virial_parameter'] = (3*((sigma_v_cd)**2)*cat_cd['radius_pc'])/(cat_cd['Mass']*G)
#cat_cd['virial_parameter'].unit = ''

##    Save Catalogs   ##
cat.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_smoothed_cat.fits',overwrite=True)
#cat_cd.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_column_density_smoothed_cat.fits',overwrite=True)

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Loading the ppv_catalog
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#cat = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_cat.fits')
#cat_cd = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_column_density_cat.fits')

#
# Variables for plotting leaves
#


L_CO = np.array([])
radius_pc = np.array([])
vrms = np.array([])
mass = np.array([])
Sigma_exact = np.array([])
Sigma_ellipse = np.array([])
alpha = np.array([])


#radius_pc_cd = np.array([])
#vrms_cd = np.array([])
#mass_cd = np.array([])
#Sigma_exact_cd = np.array([])
#Sigma_ellipse_cd = np.array([])
#alpha_cd = np.array([])

for i in leaf_id:
    L_CO = np.append(L_CO, cat[int(i)]['luminosity'])
    radius_pc = np.append(radius_pc,cat[int(i)]['radius_pc'])
    vrms = np.append(vrms,sigma_v[int(i)])
    mass = np.append(mass, cat[int(i)]['M_lum'])
    #mass = np.append(mass, cat[i]['Mass'])
    Sigma_exact = np.append(Sigma_exact, cat[int(i)]['Sigma_exact'])
    Sigma_ellipse = np.append(Sigma_ellipse,cat[int(i)]['Sigma_ellipse'])
    alpha = np.append(alpha,cat[int(i)]['virial_parameter'])
#for i in leaf_cd_id:
#    radius_pc_cd = np.append(radius_pc_cd,cat_cd[int(i)]['radius_pc'])
#    vrms_cd = np.append(vrms_cd,sigma_v_cd[int(i)])
#    mass_cd = np.append(mass_cd, cat_cd[int(i)]['Mass'])
#    Sigma_exact_cd = np.append(Sigma_exact_cd, cat_cd[int(i)]['Sigma_exact'])
#    Sigma_ellipse_cd = np.append(Sigma_ellipse_cd,cat_cd[int(i)]['Sigma_ellipse'])
#    alpha_cd = np.append(alpha_cd,cat_cd[int(i)]['virial_parameter'])


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% CHIMPS and ATLASGAL Data Load
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
cmp_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/table4.fits")
cmp_mass = cmp_data['Mass']
cmp_radius = cmp_data['R_eq']
cmp_alpha = cmp_data['alpha'].astype(np.float64)
cmp_flag = cmp_data['Flag_R']
cmp_vrms = cmp_data['sigma_v']

AGAL_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/asu.fits")
AGAL_mass = AGAL_data['logMclump'].astype(np.float64)
AGAL_alpha = AGAL_data['alpha'].astype(np.float64)
AGAL_radius = AGAL_data['Rad'].astype(np.float64)
AGAL_vrms = AGAL_data['NH3LW'].astype(np.float64)

plt.rcParams.update({'font.size':10})
fig,ax = plt.subplots(ncols=4,nrows=1,figsize=(15,4))

sns.kdeplot(data=np.log10(mass),ax=ax[0],c='goldenrod',fill=1,alpha=0.5,label='G305')
#sns.kdeplot(data=np.log10(mass_cd),ax=ax[0],c='gray',fill=0,alpha=0.5,label='N$_{cd}$')
sns.kdeplot(data=np.log10(cmp_mass[(cmp_flag==3)]),ax=ax[0],c='blue',fill=1,alpha=0.5,label='CHIMPS')
sns.kdeplot(data=AGAL_mass[AGAL_mass>0],ax=ax[0],c='maroon',fill=0,alpha=0.5,linestyle='dashed',label='AGAL')
ax[0].set_xlabel(r'log$_{10}$(M/M$_{\odot}$)')
ax[0].legend(loc='upper left')

sns.kdeplot(data=np.log10(radius_pc),ax=ax[1],c='goldenrod',fill=1,alpha=0.5)
#sns.kdeplot(data=np.log10(radius_pc_cd),ax=ax[1],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_radius[(cmp_flag==3)]),ax=ax[1],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_radius[AGAL_radius>0]),ax=ax[1],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[1].set_xlabel(r'log$_{10}$(R$_{eq}$/pc)')
ax[1].set_ylabel('')

sns.kdeplot(data=np.log10(alpha),ax=ax[2],c='goldenrod',fill=1,alpha=0.5)
#sns.kdeplot(data=np.log10(alpha_cd),ax=ax[2],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_alpha[(cmp_flag==3)]),ax=ax[2],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_alpha[AGAL_alpha>0]),ax=ax[2],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[2].set_xlabel(r'log$_{10}(\alpha_{vir})$')
ax[2].set_ylabel('')

sns.kdeplot(data=np.log10(vrms),ax=ax[3],c='goldenrod',fill=1,alpha=0.5)
#sns.kdeplot(data=np.log10(vrms_cd),ax=ax[3],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_vrms[(cmp_flag==3)]),ax=ax[3],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_vrms[AGAL_vrms>0]),ax=ax[3],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[3].set_xlabel(r'log$_{10}\,v_{rms}$')
ax[3].set_ylabel('')


plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/clump_properties_comparison.png",format="png", overwrite=True)
