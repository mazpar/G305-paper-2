#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:42:52 2021

Compare G305 Clumps Properties with Galactic Clumps

@author: pmazumdar
"""
from astropy.table import Table
from astrodendro import Dendrogram, ppv_catalog, structure
from astropy import units as u
from astropy import wcs
from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astrodendro.pruning import all_true, min_vchan, min_delta, min_area
from astropy import constants as const
from spectral_cube import SpectralCube
import radio_beam
import aplpy
import seaborn as sns
import scipy.stats as sst


figpath = "/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/"





#&%&%&%&%&%&%&%&%&%&%&%&%&%&
#&%
#&% Load the fits file
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&
cube = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO_resample.fits')
cube_cd = SpectralCube.read('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal.fits')

hdu_13CO = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO-moment0.fits')[0]
mask_nan = ~np.isnan(hdu_13CO.data) # only include non nan pixels

masked_cube = cube.with_mask(mask_nan)   # apply mask to spectral cube
#masked_cube_cd = cube_cd.with_mask(mask_nan)   # apply mask to spectral cube

data = masked_cube.hdu.data
hd = masked_cube.hdu.header
wc = wcs.WCS(hd)
#data_cd = masked_cube_cd.hdu.data
#hd_cd = masked_cube_cd.hdu.header
#wc_cd = wcs.WCS(hd_cd)





############################
##    Read the 8um Map    ##
############################
hdu_glm8 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
data_glm8 = hdu_glm8.data
w_glm8 = wcs.WCS(hdu_glm8.header)





#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% DENDROGRAM LOAD
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
d = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_dendro.fits')
d_sm = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_smoothed_dendro.fits')






#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Loading the ppv_catalog
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
cat = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_cat.fits')
cat_sm = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_smoothed_cat.fits')


dists = 3800*u.parsec # Distance of the source.
dist_gc = 6598.5452942296305*u.parsec
x2 = 18 # x2 = X_CO / 2E20 [cm^-2 / K km s^-1] , where X_CO is a CO-H2 conversion factor.
#x2 = 1 # for 12CO
sigma_v = np.array(cat['v_rms'])
sigma_x = np.array(cat['radius'])*(((1*u.arcsecond).to(u.rad)).value)*dists.value
eta = 1.91 # conversion factor. R = eta * sigma_r
G = 4.302e-3 # units of pc (M_sun)^-1 (km/s)^2
deltax_pc = abs(np.pi/180.*hd['CDELT1']*dists.value) # delta x in pc
deltay_pc = abs(np.pi/180.*hd['CDELT2']*dists.value) # delta y in pc
deltav_kms = abs(hd['CDELT3']/1000.) # vel res in kmps
nu_12CO = 345.79598990*u.GHz
nu_13CO = 330.58796530*u.GHz
delta_nu_12 = 0.0011534512649414282*u.GHz
delta_nu_13 = 0.0011027227552869259*u.GHz
sigma_majs = cat['major_sigma']
sigma_mins = cat['minor_sigma']
R12_13 = (6.21*(dist_gc.value/1000.0))+18.71 # C12/C13 abundance ratio (Milam et al. 2005)
R12 = 8.5e-5 # C12/H2 abundance (Frerking et al.)
R13_inv =  R12_13/R12
mu = 2.72 # Average H2 mass including He fraction
mp = 8.4089382e-58*u.solMass # Proton Mass
K_Jy = 37
kmps_to_GHz_13 = 1000*(nu_13CO.si.value)/(const.c.value)
conversion_sollum = K_Jy*kmps_to_GHz_13



# Structure tags 


##    adding a surface density column to the catalog
cat['n_H2'] = np.zeros(len(cat),dtype=float)
#cat['n_H2'] = (cat['Mass']/(mu*mp*(4/3)*np.pi*cat['radius_pc']*cat['radius_pc']*cat['radius_pc'])).to(1/u.cm**3)
cat['n_H2'] = 15.1*cat['Mass']/((4/3)*np.pi*cat['radius_pc']*cat['radius_pc']*cat['radius_pc'])
cat['n_H2'].unit = 1/(u.cm**3)


# Normal Catalog
cat['Structure Tag'] = np.zeros(len(cat),dtype=str)

for l in d.all_structures:
  # Add a tag to each structure of leaf or branch or trunk
  if l.is_leaf:
    cat[l.idx]['Structure Tag']="l"
  else:
    cat[l.idx]['Structure Tag']="b"


cat['Trunk Tag'] = np.zeros(len(cat),dtype=int)
cat['Trunk Tag'] = 0
for t in d.trunk:
  cat[t.idx]['Trunk Tag']=1  

cat.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_cat.fits',overwrite=True)


# Smoothed Catalog

##    adding a surface density column to the catalog
cat_sm['n_H2'] = np.zeros(len(cat_sm),dtype=float)
cat_sm['n_H2'] = (cat_sm['Mass']/(mu*mp*(4/3)*np.pi*cat_sm['radius_pc']*cat_sm['radius_pc']*cat_sm['radius_pc'])).to(1/u.cm**3)
cat_sm['n_H2'].unit = 1/(u.cm**3)


cat_sm['Structure Tag'] = np.zeros(len(cat_sm),dtype=str)

for l in d_sm.all_structures:
  # Add a tag to each structure of leaf or branch or trunk
  if l.is_leaf:
    cat_sm[l.idx]['Structure Tag']="l"
  else:
    cat_sm[l.idx]['Structure Tag']="b"


cat_sm['Trunk Tag'] = np.zeros(len(cat_sm),dtype=int)
cat_sm['Trunk Tag'] = 0
for t in d_sm.trunk:
  cat_sm[t.idx]['Trunk Tag']=1  

cat_sm.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_smoothed_cat.fits',overwrite=True)




#
#Load Catalog properties
#

sr_type = cat['Structure Tag']
trunk_tag = cat['Trunk Tag']

LAS_mass = cat['Mass'].astype(np.float64)
LAS_mass_leaves = LAS_mass[(sr_type=="l")]
LAS_rad = cat['radius_pc'].astype(np.float64)
LAS_rad_leaves = LAS_rad[(sr_type=="l")][LAS_mass_leaves>0]
LAS_vrms = cat['v_rms'].astype(np.float64)
LAS_vrms_leaves = LAS_vrms[(sr_type=="l")][LAS_mass_leaves>0]
LAS_nH2 = cat['n_H2'].astype(np.float64)
LAS_nH2_leaves = LAS_nH2[(sr_type=="l")][LAS_mass_leaves>0]
LAS_alpha = cat['virial_parameter'].astype(np.float64)
LAS_alpha_leaves = LAS_alpha[(sr_type=="l")][LAS_mass_leaves>0]
LAS_lum_leaves = cat['luminosity'][(sr_type=="l")][LAS_mass_leaves>0]
LAS_mass_leaves = LAS_mass_leaves[LAS_mass_leaves>0]
LAS_mass_leaves_sol = LAS_mass_leaves*conversion_sollum


sr_type_sm = cat_sm['Structure Tag']
trunk_tag_sm = cat_sm['Trunk Tag']

LAS_mass_sm = cat_sm['Mass'].astype(np.float64)
LAS_mass_sm_leaves = LAS_mass_sm[(sr_type_sm=="l")]
LAS_rad_sm = cat_sm['radius_pc'].astype(np.float64)
LAS_rad_leaves_sm = LAS_rad_sm[(sr_type_sm=="l")][LAS_mass_sm_leaves>0]
LAS_vrms_sm = cat_sm['v_rms'].astype(np.float64)
LAS_vrms_leaves_sm = LAS_vrms_sm[(sr_type_sm=="l")][LAS_mass_sm_leaves>0]
LAS_nH2_sm = cat_sm['n_H2'].astype(np.float64)
LAS_nH2_leaves_sm = LAS_nH2_sm[(sr_type_sm=="l")][LAS_mass_sm_leaves>0]
LAS_alpha_sm = cat_sm['virial_parameter'].astype(np.float64)
LAS_alpha_sm_leaves = LAS_alpha_sm[(sr_type_sm=="l")][LAS_mass_sm_leaves>0]
LAS_mass_sm_leaves = LAS_mass_sm_leaves[LAS_mass_sm_leaves>0]
LAS_mass_sm_leaves_sol = LAS_mass_sm_leaves*conversion_sollum



cat_g333 = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G333_13CO_smoothed_cat.fits')
sr_type_g333 = cat_g333['Structure Tag']
G333_mass = cat_g333['Mass'].astype(np.float64)
G333_mass_leaves = G333_mass[(sr_type_g333=="l")]


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% CHIMPS and ATLASGAL Data Load and Compare
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
cmp_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/table4.fits")
cmp_mass = cmp_data['Mass'].astype(np.float64)
cmp_radius = cmp_data['R_eq'].astype(np.float64)
cmp_alpha = cmp_data['alpha'].astype(np.float64)
cmp_flag = cmp_data['Flag_R'].astype(np.float64)
cmp_vrms = cmp_data['sigma_v'].astype(np.float64)
cmp_nH2 = cmp_data['nH2'].astype(np.float64)


AGAL_2_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/asu.fits")
#AGAL_mass = AGAL_data['logMclump'].astype(np.float64)
#AGAL_alpha = AGAL_2_data['alpha'][AGAL_2_data['Clump']==AGAL_names]
#AGAL_radius = AGAL_data['Rad'].astype(np.float64)
#AGAL_vrms = AGAL_data['NH3LW'].astype(np.float64)
#AGAL_lum = AGAL_data['logLbol'].astype(np.float64)

AGAL_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/csc_properties_lasma.csv")
AGAL_names = AGAL_data['csc_name'][AGAL_data['clump_mass']!="NULL"]
AGAL_mass = AGAL_data['clump_mass'][AGAL_data['clump_mass']!="NULL"].astype(np.float64)
AGAL_radius = AGAL_data['radius_pc'][AGAL_data['clump_mass']!="NULL"].astype(np.float64)
AGAL_lmratio = AGAL_data['l/m'][AGAL_data['clump_mass']!="NULL"].astype(np.float64)
AGAL_dist = AGAL_data['distance'][AGAL_data['clump_mass']!="NULL"].astype(np.float64)
AGAL_lum = AGAL_data['lum_bol'][AGAL_data['clump_mass']!="NULL"].astype(np.float64)

AGAL_maj = AGAL_data['semi_maj'][AGAL_data['clump_mass']!="NULL"]
AGAL_min = AGAL_data['semi_min'][AGAL_data['clump_mass']!="NULL"]
deltax_pc = abs((np.pi/180.)*(0.00028*AGAL_dist*1000))
deltay_pc = abs((np.pi/180.)*(0.00028*AGAL_dist*1000))
AGAL_rad_pc = 1.91*np.sqrt((AGAL_maj*deltax_pc)*(AGAL_min*deltay_pc))


AGAL_alpha = np.array([])
AGAL_vrms = np.array([])
for name in AGAL_names:
  alpha = AGAL_2_data['alpha'][AGAL_2_data['Clump']==name]
  vrms = AGAL_2_data['NH3LW'][AGAL_2_data['Clump']==name]
  AGAL_alpha = np.append(AGAL_alpha,alpha)
  AGAL_vrms = np.append(AGAL_vrms,vrms)
  

SED_data = Table.read("/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/SEDIGISM_complete_cloud_catalogue_05May20_final.fits")
SED_mass = SED_data['Mass']
SED_den = SED_data['radius_eq_as']
SED_dist = SED_data['dist_kpc']

  
##  Only G333 Values
AGAL_mass_G333 = np.array([])
#AGAL_lum_G333 = np.array([])


for i in range(len(AGAL_data)):
  if (AGAL_dist[i]<3.62) & (AGAL_dist[i]>3.55):
    if "AGAL333" in str(AGAL_names[i]):
      AGAL_mass_G333 = np.append(AGAL_mass_G333,AGAL_mass[i])
      #AGAL_lum_G333 = np.append(AGAL_lum_G333,AGAL_lum[i])
      #print(AGAL_data['VLSR'][i])
     
      
## Dist Limited Sample
AGAL_names_dlim = AGAL_names[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
#AGAL_data_dlim = AGAL_data[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_dist_dlim = AGAL_dist[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_radius_dlim = AGAL_radius[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_rad_pc_dlim = AGAL_rad_pc[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_mass_dlim = AGAL_mass[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_alpha_dlim = AGAL_alpha[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_vrms_dlim = AGAL_vrms[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
AGAL_lmratio_dlim = AGAL_lmratio[(AGAL_dist>3.6)&(AGAL_dist<4.5)]
#AGAL_lum_dlim = AGAL_lum[(AGAL_data['Dist']>3.6)&(AGAL_data['Dist']<4.5)]
#AGAL_nH2_dlim = 15.1 * 10**(AGAL_mass_dlim) / ((4/3)*(np.pi*AGAL_radius_dlim**3))
AGAL_nH2_dlim = 15.1 * (AGAL_mass_dlim) / ((4/3)*(np.pi*AGAL_rad_pc_dlim**3))


cmp_radius_dlim = cmp_radius[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]
cmp_vrms_dlim = cmp_vrms[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]
cmp_mass_dlim = cmp_mass[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]
cmp_alpha_dlim = cmp_alpha[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]
cmp_flag_dlim = cmp_flag[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]
cmp_nH2_dlim = cmp_nH2[(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]


##  Only G305 Values
AGAL_mass_G305 = np.array([])
AGAL_mass_non_G305 = np.array([])

AGAL_lmratio_G305 = np.array([])
AGAL_lmratio_non_G305 = np.array([])

AGAL_rad_G305 = np.array([])
AGAL_rad_non_G305 = np.array([])

AGAL_vrms_non_G305 = np.array([])
AGAL_nH2_non_G305 = np.array([])
AGAL_alpha_non_G305 = np.array([])

for i in range(len(AGAL_dist_dlim)):
  if AGAL_dist_dlim[i]==3.8:
    if "AGAL305" in str(AGAL_names_dlim[i]):
      AGAL_mass_G305 = np.append(AGAL_mass_G305,AGAL_mass_dlim[i])
      AGAL_lmratio_G305 = np.append(AGAL_lmratio_G305,AGAL_lmratio_dlim[i])
      AGAL_rad_G305 = np.append(AGAL_rad_G305,AGAL_rad_pc_dlim[i])
    else:
      AGAL_rad_non_G305 = np.append(AGAL_rad_non_G305,AGAL_rad_pc_dlim[i])
      AGAL_lmratio_non_G305 = np.append(AGAL_lmratio_non_G305,AGAL_lmratio_dlim[i])
      AGAL_mass_non_G305 = np.append(AGAL_mass_non_G305, AGAL_mass_dlim[i])
      AGAL_vrms_non_G305 = np.append(AGAL_vrms_non_G305, AGAL_vrms_dlim[i])
      AGAL_nH2_non_G305 = np.append(AGAL_nH2_non_G305, AGAL_nH2_dlim[i])
      AGAL_alpha_non_G305 = np.append(AGAL_alpha_non_G305, AGAL_alpha_dlim[i])
  else:
    AGAL_rad_non_G305 = np.append(AGAL_rad_non_G305,AGAL_rad_pc_dlim[i])
    AGAL_lmratio_non_G305 = np.append(AGAL_lmratio_non_G305,AGAL_lmratio_dlim[i])
    AGAL_mass_non_G305 = np.append(AGAL_mass_non_G305, AGAL_mass_dlim[i])
    AGAL_vrms_non_G305 = np.append(AGAL_vrms_non_G305, AGAL_vrms_dlim[i])
    AGAL_nH2_non_G305 = np.append(AGAL_nH2_non_G305, AGAL_nH2_dlim[i])
    AGAL_alpha_non_G305 = np.append(AGAL_alpha_non_G305, AGAL_alpha_dlim[i])


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% PLOTTING
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


#
#  Compare Distribution of Clump Properties
#

def plot_clump_props_AGAL(save=False):
  plt.rcParams.update({'font.size':16})
  fig,ax = plt.subplots(ncols=4,nrows=1,figsize=(20,5))

  #rad_corr = 1.91/2.4  
  # Radius
  sns.histplot(np.log10(LAS_rad_leaves),\
              ax=ax[0],\
              #cut=0,\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(AGAL_rad_non_G305),\
              ax=ax[0],\
              #cut=0,\
              kde=True,\
              stat='density',\
              color = 'tab:blue')
  #sns.histplot(np.log10(rad_corr*AGAL_radius_dlim[AGAL_mass_dlim>0.0]),\
  #            ax=ax[0],\
  #            #cut=0,\
  #            kde=True,\
  #            stat='density',\
  #            color = 'tab:blue')
  
  # V+rms
  sns.histplot(np.log10(LAS_vrms_leaves),\
              ax=ax[1],\
              #cut=0,\
              binwidth=0.05,\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(AGAL_vrms_non_G305[AGAL_vrms_non_G305>0]/2.35),\
              ax=ax[1],\
              #cut=0,\
              binwidth=0.05,\
              kde=True,\
              stat='density',\
              color = 'tab:blue')
  
  
  #Mass
  sns.histplot(np.log10(LAS_mass_leaves),\
              ax=ax[2],\
              #cut=0,\
              binwidth=0.25,\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(AGAL_mass_non_G305),\
              ax=ax[2],\
              binwidth=0.25,\
              #cut=0,\
              kde=True,\
              stat='density',\
              color = 'tab:blue')
  

  #nH2
  #sns.histplot(np.log10(LAS_nH2_leaves),\
  #            ax=ax[3],\
  #            #cut=0,\
  #            kde=True,\
  #            stat='density',\
  #            color = 'tab:red')
  #sns.histplot(np.log10(AGAL_nH2_non_G305),\
  #            ax=ax[3],\
  #            #cut=0,\
  #            kde=True,\
  #            stat='density',\
  #            color = 'tab:blue')
  
  #Virial Parameter
  sns.histplot(np.log10(LAS_alpha_leaves),\
              ax=ax[3],\
              #cut=0,\
              binwidth=0.25,\
              kde=True,\
              stat='density',\
              color = 'tab:red',\
              label='LAsMA')
  sns.histplot(np.log10(AGAL_alpha_non_G305[AGAL_vrms_non_G305>0.0]),\
              ax=ax[3],\
              #cut=0,\
              binwidth=0.25,\
              kde=True,\
              stat='density',\
              color = 'tab:blue',\
              label = 'ATLASGAL')
  
  #labels and legends
  ax[0].set_xlabel('Log(Radius) ($pc$)')
  ax[1].set_xlabel(r'Log($\sigma_{\rm{rms}}$) ($km\,s^{-1}$)')
  ax[2].set_xlabel(r'Log(Mass) ($M_{\odot}$)')
  #ax[3].set_xlabel(r'Log($n_{H_2}$) ($cm^{-3}$)')
  ax[3].set_xlabel(r'Log($\alpha$)')
  
  ax[0].set_ylabel('Density')
  ax[1].set_ylabel('')
  ax[2].set_ylabel('')
  ax[3].set_ylabel('')
  #ax[4].set_ylabel('')
  
  ax[3].legend(loc='upper right',frameon=False, labelspacing=0.8, fontsize=16)
  
  plt.tight_layout()

  if save:
    plt.savefig(figpath+"Clump_Properties_Comparison_AGAL.pdf",format='pdf',dpi=800)


def plot_clump_props_CHIMPS(save=False):
  plt.rcParams.update({'font.size':16})
  fig,ax = plt.subplots(ncols=4,nrows=1,figsize=(20,5))

  # Radius
  sns.histplot(np.log10(LAS_rad_leaves_sm),\
              ax=ax[0],\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(cmp_radius_dlim[cmp_flag_dlim==3]),\
              ax=ax[0],\
              kde=True,\
              stat='density',\
              color = 'goldenrod')
  
  # V+rms
  sns.histplot(np.log10(LAS_vrms_leaves_sm),\
              ax=ax[1],\
              binwidth=0.05,\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(cmp_vrms_dlim[cmp_flag_dlim==3]),\
              ax=ax[1],\
              binwidth=0.05,\
              kde=True,\
              stat='density',\
              color = 'goldenrod')
  
  
  #Mass
  sns.histplot(np.log10(LAS_mass_sm_leaves),\
              ax=ax[2],\
              binwidth=0.25,\
              kde=True,\
              stat='density',\
              color = 'tab:red')
  sns.histplot(np.log10(cmp_mass_dlim[cmp_flag_dlim==3]),\
              ax=ax[2],\
              binwidth=0.25,\
              kde=True,\
              stat='density',\
              color = 'goldenrod')

  #nH2
  #sns.histplot(np.log10(LAS_nH2_leaves_sm),\
  #            ax=ax[3],\
  #            kde=True,\
  #            stat='density',\
  #            color = 'tab:red')
  #sns.histplot(np.log10(cmp_nH2_dlim[cmp_flag_dlim==3]),\
  #            ax=ax[3],\
  #            kde=True,\
  #            stat='density',\
  #            color = 'goldenrod')


  #Virial Parameter
  sns.histplot(np.log10(LAS_alpha_sm_leaves),\
              ax=ax[3],\
              binwidth=0.125,\
              kde=True,\
              stat='density',\
              color = 'tab:red',\
              label='LAsMA')
  sns.histplot(np.log10(cmp_alpha_dlim[cmp_flag_dlim==3]),\
              ax=ax[3],\
              binwidth=0.125,\
              kde=True,\
              stat='density',\
              color = 'goldenrod',\
              label='CHIMPS')
  
  
  #labels and legends
  ax[0].set_xlabel('Log(Radius) ($pc$)')
  ax[1].set_xlabel(r'Log($\sigma_{\rm{rms}}$) ($km\,s^{-1}$)')
  ax[2].set_xlabel(r'Log(Mass) ($M_{\odot}$)')
  #ax[3].set_xlabel(r'Log($n_{H_2}$) ($cm^{-3}$)')
  ax[3].set_xlabel(r'Log($\alpha$)')
  
  ax[0].set_ylabel('Density')
  ax[1].set_ylabel('')
  ax[2].set_ylabel('')
  ax[3].set_ylabel('')
  #ax[4].set_ylabel('')
  
  ax[3].legend(loc='upper right',frameon=False, labelspacing=0.8, fontsize=16)
  
  plt.tight_layout()

  if save:
    plt.savefig(figpath+"Clump_Properties_Comparison_CHIMPS.pdf",dpi=800)



# KS_Statistics of the properties

ks_rad_agal = sst.ks_2samp(LAS_rad_leaves, AGAL_rad_non_G305)
ks_rad_cmp = sst.ks_2samp(LAS_rad_leaves_sm,cmp_radius_dlim[cmp_flag_dlim==3])
ks_vrms_agal = sst.ks_2samp(LAS_vrms_leaves, AGAL_vrms_non_G305[AGAL_vrms_non_G305>0])
ks_vrms_cmp = sst.ks_2samp(LAS_vrms_leaves_sm,cmp_vrms_dlim[cmp_flag_dlim==3])
ks_mass_agal = sst.ks_2samp(LAS_mass_leaves,AGAL_mass_non_G305)
ks_mass_cmp = sst.ks_2samp(LAS_mass_sm_leaves,cmp_mass_dlim[cmp_flag_dlim==3])
ks_nH2_agal = sst.ks_2samp(LAS_nH2_leaves,AGAL_nH2_non_G305)
ks_nH2_cmp = sst.ks_2samp(LAS_nH2_leaves_sm,cmp_nH2_dlim[cmp_flag_dlim==3])
ks_vir_agal = sst.ks_2samp(LAS_alpha_leaves,AGAL_alpha_non_G305[AGAL_alpha_non_G305>0])
ks_vir_cmp = sst.ks_2samp(LAS_alpha_sm_leaves,cmp_alpha_dlim[cmp_flag_dlim==3])


#
#  CDF of Clump Mass
#
plt.rcParams.update({'font.size':18})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(9,9))

sns.ecdfplot(data=np.log10(LAS_mass_leaves[LAS_mass_leaves!=0]),\
             complementary=True,\
             ax=ax,\
             c='red',\
             label='G305 Clumps')
sns.ecdfplot(data=np.log10(AGAL_mass_non_G305),\
             complementary=True,\
             ax=ax,\
             c='blue',\
             label='AGAL Clumps')

sns.ecdfplot(data=np.log10(LAS_mass_sm_leaves[LAS_mass_sm_leaves>1]),\
             complementary=True,\
             ax=ax,\
             c='red',\
             linestyle='dotted',\
             label='G305 Smoothed Clumps')
sns.ecdfplot(data=np.log10(cmp_mass[(cmp_flag==3)&(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]),\
             complementary=True,\
             ax=ax,\
             c='blue',\
             linestyle='dotted',\
             label='CHIMPS Clumps')
#sns.ecdfplot(data=np.log10(G333_mass_leaves[G333_mass_leaves!=0]),\
#             complementary=True,\
#             ax=ax,\
#             c='magenta',\
#             linestyle='dotted',\
#             label='G333 Clumps')

ax.set_xlabel(r'log$_{10}$(M/M$_{\odot}$)')
ax.set_ylabel(r'CDF [N(>M)]')
ax.set_yscale('log')
ax.legend(loc='lower left',frameon=False, labelspacing=0.8)

plt.tight_layout()
plt.savefig(figpath+"Clump_Mass_CDF_comparison.pdf", format="pdf", overwrite=True, dpi=800)
#plt.savefig(figpath+"Clump_Mass_CDF_comparison_with_G333.pdf", format="pdf", overwrite=True, dpi=200)

ks_stat_AGAL = sst.ks_2samp(np.log10(AGAL_mass_dlim),np.log10(LAS_mass_leaves[LAS_mass_leaves!=0]))
print("KS_AGAL_pvalue ="+str(ks_stat_AGAL[1]))
ks_stat_CMP = sst.ks_2samp(np.log10(cmp_mass[(cmp_flag==3)&(cmp_data['Dist']>3.5)&(cmp_data['Dist']<4.5)]),\
                                             np.log10(LAS_mass_sm_leaves[LAS_mass_sm_leaves>1]))
print("KS_CMP_pvalue ="+str(ks_stat_CMP[1]))




#
#  CDF of Luminosity to Mass Ratio
#

fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(9,9))

sns.ecdfplot(data=AGAL_lmratio_G305,complementary=True,ax=ax,color='goldenrod',label='ATLASGAL G305 Clumps')
#sns.ecdfplot(data=AGAL_lum_G333-AGAL_mass_G333,complementary=True,ax=ax,color='magenta',linestyle='dotted',label='ATLASGAL G333 Clumps')
#sns.ecdfplot(data=AGAL_lum_dlim-AGAL_mass_dlim,complementary=True,ax=ax,color='blue', label='ATLASGAL Galactic Clumps')
#sns.ecdfplot(data=np.log((LAS_lum_leaves/LAS_mass_leaves)),complementary=True,ax=ax,color='red', label='LAsMA G305 Clumps')
sns.ecdfplot(data=AGAL_lmratio_non_G305,complementary=True,ax=ax,color='blue', label='ATLASGAL Galactic Clumps')

ax.set_xlabel(r'log(L$_{\odot}$/M$_{\odot}$)')
ax.legend(frameon=False, labelspacing=0.8)
ax.set_yscale('log')
ax.set_xscale('log')

plt.tight_layout()
plt.savefig(figpath+"lum_to_mass_cdf.pdf",format="pdf",dpi=800)
#plt.savefig(figpath+"lum_to_mass_cdf_with_G333.pdf",format="pdf",dpi=200)


#ks_stat = sst.ks_2samp(AGAL_lum_G305-AGAL_mass_G305,AGAL_lum_dlim-AGAL_mass_dlim)
#ks_stat_2 = sst.ks_2samp(AGAL_lum_G333-AGAL_mass_G333,AGAL_lum_dlim-AGAL_mass_dlim)
ks_stat_lm = sst.ks_2samp(AGAL_lmratio_G305,AGAL_lmratio_non_G305)

print("KS_pvalue ="+str(ks_stat[1]))
print("KS_pvalue ="+str(ks_stat_2[1]))












