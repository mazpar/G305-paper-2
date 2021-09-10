#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 13:26:02 2021

Scaling Relations Plotting

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

hdu_13CO = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO-moment0.fits')[0]
mask_nan = ~np.isnan(hdu_13CO.data) # only include non nan pixels

masked_cube = cube.with_mask(mask_nan)   # apply mask to spectral cube


data = masked_cube.hdu.data
hd = masked_cube.hdu.header
wc = wcs.WCS(hd)






############################
##    Read the 8um Map    ##
############################

hdu_glm8 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
data_glm8 = hdu_glm8.data
w_glm8 = wcs.WCS(hdu_glm8.header)







########################################
##    Make Masks based on 8um flux    ##
########################################

#---- Define Threshold values and make mask

threshold = 120 #lower threshold
glm8_mask = data_glm8 > threshold # 8um mask








##################################################
##    Make Masks for Radio Continuum Sources    ##
##################################################

from astropy.coordinates import SkyCoord
from regions import EllipseSkyRegion, CircleSkyRegion, PixCoord, write_fits_region

workfile = "/home/pmazumdar/Documents/LASMA/Dendrogram/tables_etc/Radio_Continuum_Sources.dat"

RC_sources = np.genfromtxt(workfile,usecols=(1,2,3)) # major dimension as rotation angle not given

RC_regions = [CircleSkyRegion\
              (center=SkyCoord(x,y,unit='deg', frame='galactic'),radius=r*u.arcsecond)\
              for x,y,r in RC_sources]

RC_regions_pix = [RC_region.to_pixel(w_glm8) for RC_region in RC_regions]

RC_masks = [RC_region_pix.to_mask().to_image(hdu_glm8.data.shape)\
            for RC_region_pix in RC_regions_pix]

RC_mask_sum = RC_masks[0]

for i in range(len(RC_masks)-1):
  RC_mask_sum = np.logical_or(RC_mask_sum,RC_masks[i+1])

RC_mask_sum = RC_mask_sum.astype(float)









#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% START HERE IF DENDROGRAM ALREADY RUN
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

##    Custom Definitions for the Dendrogram    ##
rms = 0.12 # rms noise
bmaj = hd['bmaj'] # beam_major
bmin = hd['bmin'] # beam_minor
cdelt1 = hd['cdelt1'] # delta_x
cdelt2 = hd['cdelt2'] # delta_y
deltav_kms = abs(hd['CDELT3']/1000.) # vel res in kmps
ppb = abs((bmaj*bmin)/(cdelt1*cdelt2)*2*np.pi/(8*np.log(2))) # pixel_per_beam

d = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_dendro.fits')








#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Loading the ppv_catalog
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

cat = Table.read('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_cat.fits')

eta = 1.91 # conversion factor. R = eta * sigma_r
G = 4.302e-3 # units of pc (M_sun)^-1 (km/s)^2
dists = 3800*u.parsec # Distance of the source.
dist_gc = 6598.5452942296305*u.parsec

R12_13 = (6.21*(dist_gc.value/1000.0))+18.71 # C12/C13 abundance ratio (Milam et al. 2005)
R12 = 8.5e-5 # C12/H2 abundance (Frerking et al.)
R13_inv =  R12_13/R12
mu = 2.72 # Average H2 mass including He fraction
mp = 8.4089382e-58*u.solMass # Proton Mass
nu_12CO = 345.79598990*u.GHz
nu_13CO = 330.58796530*u.GHz
delta_nu_12 = 0.0011534512649414282*u.GHz
delta_nu_13 = 0.0011027227552869259*u.GHz

x2 = 18 # x2 = X_CO / 2E20 [cm^-2 / K km s^-1] , where X_CO is a CO-H2 conversion factor.
#x2 = 1 # for 12CO

sigma_v = np.array(cat['v_rms'])
sigma_x = np.array(cat['radius'])*(((1*u.arcsecond).to(u.rad)).value)*dists.value

deltax_pc = abs(np.pi/180.*hd['CDELT1']*dists.value) # delta x in pc
deltay_pc = abs(np.pi/180.*hd['CDELT2']*dists.value) # delta y in pc

sigma_majs = cat['major_sigma']
sigma_mins = cat['minor_sigma']







#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Make additions to the catalog and save them
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

##    Calculate Mass from Dust Maps


#hdu_dust = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/N_Dust-reproject.fits')[0]
#data_dust = hdu_dust.data
#data_dust[data_dust==np.inf]=-np.inf
#w_dust = wcs.WCS(hdu_dust.header)


#cat['mass_dust']=np.zeros(len(cat),dtype=float)
#cat['mass_dust'].unit=u.solMass


#for l in d.all_structures:    ## adding average 8um flux of leaves
#  leaf_mask_3d = l.get_mask()
#  leaf_mask = np.max(leaf_mask_3d,axis=0)
#  masked_dust_map = data_dust[leaf_mask]
#  N_dust = 10**(masked_dust_map)
#  mass = mu*mp*N_dust*deltav_kms*deltax_pc*deltay_pc*(u.parsec.to(u.centimeter))**2
#  total_mass = np.sum(mass)
#  cat[l.idx]['mass_dust']=total_mass.value


##    Save Catalogs after any new modifications   ##
#cat.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_12CO_cat.fits',overwrite=True)






#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Add a tag to each structure whether it is inside or outside the mask.
#&%     [Don't run it if this has been done for a given threshold]
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


cat['Mask_Tag']=np.zeros(len(cat),dtype=object)

cat['Structure Tag'] = np.zeros(len(cat),dtype=object)

cat['HII_Tag']=np.zeros(len(cat),dtype=object)


for l in d.all_structures:

  sr_mask_3d = l.get_mask()   # Make a structure mask
  sr_mask = np.max(sr_mask_3d,axis=0)  # Project it to 2D
  and_mask = sr_mask & glm8_mask  # An overlap mask of the two masks

  # 8um Tag
  if 100.*(np.sum(and_mask)/np.sum(sr_mask))>67:  # > 2/3rd area inside the mask
    cat[l.idx]['Mask_Tag']= "Inside"
  elif 100.*(np.sum(and_mask)/np.sum(sr_mask))<=10:  # > 1/10 area inside the mask
    cat[l.idx]['Mask_Tag']= "Outside"
  else:
    cat[l.idx]['Mask_Tag']= "Partly Inside"

  # HII Tag
  if 100.*(np.sum(and_mask)/np.sum(sr_mask))>50:  # > 1/2 area inside the mask
    cat[l.idx]['HII_Tag']= "Inside"
  else:
    cat[l.idx]['HII_Tag']= "Outside"

  # Leaf or Branch or Trunk
  if l.is_leaf:
    cat[l.idx]['Structure Tag']="leaf"
  elif l.is_branch:
    cat[l.idx]['Structure Tag']="branch"
  else:
    cat[l.idx]['Structure Tag']="trunk"



cat['Trunk Tag'] = np.zeros(len(cat),dtype=object)
cat['Trunk Tag'] = 0
for t in d.trunk:
  cat[t.idx]['Trunk Tag']=1  







#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Viewing and Plotting the Dendrogram
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#v = d.viewer()
#v.show()
#p = d.plotter()

#fig,ax = plt.subplots(nrows = 1, ncols=1)
#plt.rcParams.update({"font.size":6})

#ax = fig.add_subplot(1, 1, 1)
#p.plot_tree(ax,color='seagreen',lw=0.5)

#ax.set_ylabel('$^{13}$CO Peak Intensity')
#ax.set_xlabel('Index of Structure')

#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Dendrogram_G305.eps",dpi=300)






################################
# Plotting Scaling Relations
################################
   
#         Load the catalog variables
# --------------------------------------------------------------------------------------------------

L_CO = cat['luminosity'].astype(np.float64)
radius_pc = cat['radius_pc'].astype(np.float64)
vrms = sigma_v.astype(np.float64)
mass = cat['Mass'].astype(np.float64)
Sigma_exact = cat['Sigma_exact'].astype(np.float64)
Sigma_ellipse = cat['Sigma_ellipse'].astype(np.float64)
alpha = cat['virial_parameter'].astype(np.float64)
mask_tag = cat['Mask_Tag']
HII_tag = cat['HII_Tag']
glm8_flux =  cat['8um Flux Mean'].astype(np.float64)
sr_type = cat['Structure Tag'].astype(str)
trunk_tag = cat['Trunk Tag']
y_var = vrms/np.sqrt(radius_pc)



#         Separate variables for leaves, trunks and trunks
# --------------------------------------------------------------------------------------------------
radius_pc_leaves = radius_pc[(sr_type=="leaf")]
vrms_leaves = vrms[(sr_type=="leaf")]
mask_tag_leaves = mask_tag[(sr_type=="leaf")]
HII_tag_leaves = HII_tag[(sr_type=="leaf")]
mass_leaves = mass[(sr_type=="leaf")]
L_CO_leaves = L_CO[(sr_type=="leaf")]
sigma_leaves = Sigma_ellipse[(sr_type=="leaf")]
y_var_leaves = y_var[(sr_type=="leaf")]
glm8_flux_leaves = glm8_flux[(sr_type=="leaf")]
alpha_leaves = alpha [(sr_type=="leaf")]

radius_pc_branches = radius_pc[(sr_type=="branch")&(trunk_tag==0)]
vrms_branches = vrms[(sr_type=="branch")&(trunk_tag==0)]
mask_tag_branches = mask_tag[(sr_type=="branch")&(trunk_tag==0)]
HII_tag_branches = HII_tag[(sr_type=="branch")&(trunk_tag==0)]
mass_branches = mass[(sr_type=="branch")&(trunk_tag==0)]
L_CO_branches = L_CO[(sr_type=="branch")&(trunk_tag==0)]
sigma_branches = Sigma_ellipse[(sr_type=="branch")&(trunk_tag==0)]
y_var_branches = y_var[(sr_type=="branch")&(trunk_tag==0)]
glm8_flux_branches = glm8_flux[(sr_type=="branch")&(trunk_tag==0)]
alpha_branches = alpha[(sr_type=="branch")&(trunk_tag==0)]


radius_pc_trunks = radius_pc[(sr_type=="branch")&(trunk_tag==1)]
vrms_trunks = vrms[(sr_type=="branch")&(trunk_tag==1)]
mask_tag_trunks = mask_tag[(sr_type=="branch")&(trunk_tag==1)]
HII_tag_trunks = HII_tag[(sr_type=="branch")&(trunk_tag==1)]
mass_trunks = mass[(sr_type=="branch")&(trunk_tag==1)]
L_CO_trunks = L_CO[(sr_type=="branch")&(trunk_tag==1)]
sigma_trunks = Sigma_ellipse[(sr_type=="branch")&(trunk_tag==1)]
y_var_trunks = y_var[(sr_type=="branch")&(trunk_tag==1)]
glm8_flux_trunks = glm8_flux[(sr_type=="branch")&(trunk_tag==1)]
alpha_trunks = alpha[(sr_type=="branch")&(trunk_tag==1)]



#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Plotting Section
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def plot_structure_props (save=False):

  plt.rcParams.update({'font.size':14})
  fig,ax = plt.subplots(ncols=3,nrows=1,figsize=(13,4))
  # Radius
  sns.histplot(np.log10(radius_pc_leaves),\
              ax=ax[0],\
              #cut=0,\
              binwidth=0.125,\
              kde=False,\
              color='tab:orange',\
              alpha=0.5)
  sns.histplot(np.log10(radius_pc_branches),\
              ax=ax[0],\
              #cut=0,\
              binwidth=0.125,\
              kde=False,\
              color='tab:blue',\
              alpha=0.5)

  # Vrms
  sns.histplot(np.log10(vrms_leaves),\
              ax=ax[1],\
              #cut=0,\
              binwidth=0.05,\
              kde=False,\
              color='tab:orange',\
              alpha=0.5)
  sns.histplot(np.log10(vrms_branches),\
              ax=ax[1],\
              #cut=0,\
              binwidth=0.05,\
              kde=False,\
              color='tab:blue',\
              alpha=0.5)

  # Mass
  sns.histplot(np.log10(mass_leaves),\
              ax=ax[2],\
              #cut=0,\
              binwidth=0.25,\
              kde=False,\
              color='tab:orange',\
              alpha=0.5,\
              label="Leaves")
  sns.histplot(np.log10(mass_branches),\
              ax=ax[2],\
              #cut=0,\
              binwidth=0.25,\
              kde=False,\
              color='tab:blue',\
              alpha=0.5,\
              label="Branches")

   # Virial Parameter
  #sns.histplot(np.log10(alpha_leaves),\
  #            ax=ax[3],\
  #            #cut=0,\
  #            kde=True,\
  #            color='tab:orange',\
  #            label="Leaves")
  #sns.histplot(np.log10(alpha_branches),\
  #            ax=ax[3],\
  #            #cut=0,\
  #            kde=True,\
  #            color='tab:blue',\
  #            label="Branches")


  #labels and legends
  ax[0].set_xlabel('Log(Radius) ($pc$)')
  ax[1].set_xlabel(r'Log($\sigma_{\rm{rms}}$) ($km\,s^{-1}$)')
  ax[2].set_xlabel(r'Log(Mass) ($M_{\odot}$)')
  #ax[3].set_xlabel(r'Log($\alpha$)')
  
  ax[0].set_ylabel('Source Counts')
  ax[1].set_ylabel('')
  ax[2].set_ylabel('')
  #ax[3].set_ylabel('')
  
  ax[2].legend(loc='upper right',frameon=False, labelspacing=0.8, fontsize=10)
  
  plt.tight_layout()

  if save:
    plt.savefig(figpath+"Structure_Properties_Comparison.pdf",dpi=800)




#
# Surface Mass Density versus 8um
#

plt.rcParams.update({'font.size':18})
fig = plt.figure(figsize=(9,8))
ax = plt.subplot()

#Fit a power law function to the data points:
from scipy.optimize import curve_fit
def func(x, a, b):
    return (a * (x**b))

x_l = glm8_flux_leaves[(sigma_leaves>0) & (~np.isnan(glm8_flux_leaves).data)]
y_l = sigma_leaves[(sigma_leaves>0) & (~np.isnan(glm8_flux_leaves))]
popt_l, pcov_l = curve_fit(func,x_l,y_l)
popt_b, pcov_b = curve_fit(func, glm8_flux_branches[sigma_branches>0], sigma_branches[sigma_branches>0])

perr_l = np.sqrt(np.diag(pcov_l))
perr_b = np.sqrt(np.diag(pcov_b))

xdata = np.linspace(20, 1500, 1000)
ax.plot(xdata, func(xdata, *popt_l), '-', color='tab:orange',label=r'fit: y = $%2.0f \pm %2.0f \, x^ {%3.2f \pm %3.2f}$' %(popt_l[0],perr_l[0],popt_l[1],perr_l[1]))
ax.plot(xdata, func(xdata, *popt_b), '-', color='tab:blue',label=r'fit: y = $%2.0f \pm %2.0f \, x^ {%3.2f \pm %3.2f}$' %(popt_b[0],perr_b[0],popt_b[1],perr_b[1]))

sns.kdeplot(y=Sigma_ellipse[Sigma_ellipse>0],\
            x=glm8_flux[Sigma_ellipse>0],\
            log_scale=10,\
            hue=cat['Structure Tag'][Sigma_ellipse>0],\
            ax=ax,\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
plt.ylabel(r'log$_{10}\Sigma$  $[M_{\odot}.pc^{-2}]$')
plt.xlabel(r'log$_{10}$ 8$\mu$m Mean Flux [MJy/sr]')


# create legends
labels = ["Leaves","Branches"]
for i, stype in enumerate(["orange", "tab:blue"]):
    plt.scatter([], [], c=stype, alpha=0.6, label=str(labels[i]))
plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Structure Type', loc='lower right')

plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Sigma_vs_8um.pdf",format="pdf", dvi=500,overwrite=True)




#
# Surface Mass Density vs Overlap
#

sig_in = sigma_leaves[(mask_tag_leaves=="Inside")&(sigma_leaves!=0)]
sig_pin = sigma_leaves[(mask_tag_leaves=="Partly Inside")&(sigma_leaves!=0)]
sig_out = sigma_leaves[(mask_tag_leaves=="Outside")&(sigma_leaves!=0)]

sig_br_in = sigma_branches[(mask_tag_branches=="Inside")]
sig_br_pin = sigma_branches[(mask_tag_branches=="Partly Inside")]
sig_br_out = sigma_branches[(mask_tag_branches=="Outside")]

plt.rcParams.update({'font.size':18})
fig,ax = plt.subplots(ncols=1,nrows=2,figsize=(10,15))

sns.histplot(np.log10(sig_in),\
             ax=ax[0],\
             stat='density',\
             element='step',\
             color='firebrick',\
             alpha=0.5,\
             legend = 'False',\
             label = "Mostly Inside")
             
sns.histplot(np.log10(sig_pin),\
             ax=ax[0],\
             stat='density',\
             element='step',\
             color='tab:green',\
             alpha=0.5,\
             legend = 'False',\
             label = "Partly Inside")

sns.histplot(np.log10(sig_out),\
             ax=ax[0],\
             stat='density',\
             element='step',\
             color='tab:blue',\
             alpha=0.5,\
             legend = 'False',\
             label = "Outside")

sns.histplot(np.log10(sig_br_in),\
             ax=ax[1],\
             stat='density',\
             element='step',\
             color='firebrick',\
             alpha=0.5,\
             legend = 'False')
             
sns.histplot(np.log10(sig_br_pin),\
             ax=ax[1],\
             stat='density',\
             element='step',\
             color='tab:green',\
             alpha=0.5,\
             legend = 'False')

sns.histplot(np.log10(sig_br_out),\
             ax=ax[1],\
             stat='density',\
             element='step',\
             color='tab:blue',\
             alpha=0.5,\
             legend = 'False')

ax[0].legend(frameon=False, labelspacing=1, title = 'Feedback Zone Overlap', loc= 'upper left')

ax[1].set_xlabel(r'log$_{10}\Sigma$  $[M_{\odot}.pc^{-2}]$')
ax[1].set_ylabel('Probability Density')
ax[0].set_ylabel('')

ax[0].set_title('Leaves')
ax[1].set_title('Branches')

plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Sigma_Tags.pdf", format='pdf')  
  
  
#
#  Individual Properties vs 8um
#

plt.rcParams.update({'font.size':12})
fig,ax = plt.subplots(ncols=4,nrows=1,figsize=(21,5))

sns.kdeplot(y=np.log10(Sigma_ellipse[Sigma_ellipse>0]),\
            x=np.log10(glm8_flux[Sigma_ellipse>0]),\
            hue=cat['Structure Tag'][Sigma_ellipse>0],\
            ax=ax[0],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[0].axvline(np.log10(60),linestyle='dotted',color='k')
#ax[0].axhline(np.log10(60),linestyle='dotted',color='k')
ax[0].set_ylabel(r'log$_{10}\Sigma [M_{\odot}.pc^{-2}]$')
ax[0].set_xlabel(r'log$_{10}$ 8$\mu$m Mean Flux [MJy/sr]')

sns.kdeplot(y=np.log10(alpha[alpha!=np.inf]),\
            x=np.log10(glm8_flux[alpha!=np.inf]),\
            hue=cat['Structure Tag'][alpha!=np.inf],\
            ax=ax[1],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[1].axhline(np.log10(35),linestyle='dotted',color='k')
ax[1].axhline(np.log10(2.0),linestyle='dashed',color='r', linewidth=0.5)
ax[1].text(3.1,np.log10(3),r'$\alpha_{vir}$=2',va='center',fontsize=12,fontstyle='oblique')
#ax[1].axvline(np.log10(60),linestyle='dotted',color='k')
ax[1].set_ylabel(r'log$_{10}(\alpha_{vir})$')
ax[1].set_xlabel('')

sns.kdeplot(y=np.log10(vrms),\
            x=np.log10(glm8_flux),\
            hue=cat['Structure Tag'],\
            ax=ax[2],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[2].axvline(np.log10(60),linestyle='dotted',color='k')
ax[2].set_ylabel(r'log$_{10}\,v_{rms}\, [km.s^{-1}]$')
ax[2].set_xlabel('')

sns.kdeplot(y=np.log10(y_var[y_var>0]),\
            x=np.log10(glm8_flux[y_var>0]),\
            hue=cat['Structure Tag'][y_var>0],\
            ax=ax[3],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[3].axvline(np.log10(60),linestyle='dotted',color='k')
ax[3].set_ylabel(r'$\sigma_v / R^{0.5} \, [km.s^{-1}.pc^{-0.5}]$')
ax[3].set_xlabel('')

# create legends
labels = ["Leaves","Branches"]
for i, stype in enumerate(["orange", "tab:blue"]):
    ax[0].scatter([], [], c=stype, alpha=0.6, label=str(labels[i]))
ax[0].legend(scatterpoints=1, frameon=False, labelspacing=1, title='Structure Type', loc='lower right')

plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/all_str_props_vs_8um.pdf",format="pdf", dvi=1000,overwrite=True)



#
#  Velocity Properties vs 8um
#


plt.rcParams.update({'font.size':14})
fig,ax = plt.subplots(ncols=1,nrows=2,figsize=(6,10))
sns.kdeplot(y=np.log10(vrms),\
            x=np.log10(glm8_flux),\
            hue=cat['Structure Tag'],\
            ax=ax[0],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[2].axvline(np.log10(60),linestyle='dotted',color='k')
ax[0].set_ylabel(r'log$_{10}$($\sigma_v$) $[km.s^{-1}]$')
ax[0].set_xlabel('')

sns.kdeplot(y=np.log10(y_var[y_var>0]),\
            x=np.log10(glm8_flux[y_var>0]),\
            hue=cat['Structure Tag'][y_var>0],\
            ax=ax[1],\
            bw_adjust=1.5,\
            fill=1,\
            alpha=0.6,\
            legend=False)
#ax[3].axvline(np.log10(60),linestyle='dotted',color='k')
ax[1].set_ylabel(r'log$_{10}$($\sigma_v / R^{0.5}$) $[km.s^{-1}.pc^{-0.5}]$')
ax[1].set_xlabel(r'log$_{10}$ 8$\mu$m Mean Flux [MJy/sr]')

# create legends
labels = ["Leaves","Branches"]
for i, stype in enumerate(["orange", "tab:blue"]):
    ax[1].scatter([], [], c=stype, alpha=0.6, label=str(labels[i]))
ax[1].legend(scatterpoints=1, frameon=False, labelspacing=1, title='Structure Type', loc='lower right')

plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/vel_vs_8um.pdf",format="pdf", dvi=1000,overwrite=True)




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
AGAL_data = AGAL_data[(AGAL_data['Dist']>3.5) & (AGAL_data['Dist']<4.5)]
AGAL_mass = AGAL_data['logMclump'].astype(np.float64)
AGAL_alpha = AGAL_data['alpha'].astype(np.float64)
AGAL_radius = AGAL_data['Rad'].astype(np.float64)
AGAL_vrms = AGAL_data['NH3LW'].astype(np.float64)







#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% CHIMPS and ATLASGAL Data Compare
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

plt.rcParams.update({'font.size':12})
fig,ax = plt.subplots(ncols=4,nrows=1,figsize=(21,5))


sns.kdeplot(data=np.log10(radius_pc_leaves),\
            ax=ax[0],c='goldenrod',fill=0, alpha=0.5,label='G305')
#sns.kdeplot(data=np.log10(cmp_radius[(cmp_flag==3)]),\
#             ax=ax[0],c='blue',fill=1,alpha=0.5,label='CHIMPS')
sns.kdeplot(data=np.log10(AGAL_radius[AGAL_radius>0]),\
            ax=ax[0],c='maroon',fill=0,alpha=0.5,label='AGAL')

ax[0].set_xlabel(r'log$_{10}$(R$_{eq}$/pc)')
ax[0].set_ylabel('')
ax[0].legend(loc='upper right')

sns.kdeplot(data=np.log10(vrms_leaves*(2*np.sqrt(2))),\
            ax=ax[1],c='goldenrod',fill=0,alpha=0.5)
#sns.kdeplot(data=np.log10(cmp_vrms[(cmp_flag==3)]),\
#             ax=ax[1],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_vrms[AGAL_vrms>0]),\
            ax=ax[1],c='maroon',fill=0,alpha=0.5)

ax[1].set_xlabel(r'log$_{10}\,v_{rms}$')
ax[1].set_ylabel('')

sns.kdeplot(data=np.log10(mass_leaves[mass_leaves>0]),\
            ax=ax[2],c='goldenrod',fill=0,alpha=0.5)
#sns.kdeplot(data=np.log10(cmp_mass[(cmp_flag==3)]),\
#             ax=ax[2],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=AGAL_mass[AGAL_mass>0],\
            ax=ax[2],c='maroon',fill=0,alpha=0.5)

ax[2].set_xlabel(r'log$_{10}$(M/M$_{\odot}$)')

sns.kdeplot(data=np.log10(alpha_leaves[alpha_leaves!=np.inf]),\
            ax=ax[3],c='goldenrod',fill=0,alpha=0.5)
#sns.kdeplot(data=np.log10(cmp_alpha[(cmp_flag==3)]),\
#             ax=ax[3],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_alpha[AGAL_alpha>0]),\
            ax=ax[3],c='maroon',fill=0,alpha=0.5)

ax[3].set_xlabel(r'log$_{10}(\alpha_{vir})$')
ax[3].set_ylabel('')






#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&%    R-V Scaling Relation
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#
#         Fit a power law to the branches and trunks    
# 

from scipy.optimize import curve_fit
def func(x,a,b):
    return a*pow(x,b)

popt_out, pcov_out = curve_fit(func,radius_pc_branches[mask_tag_branches=="Outside"],\
                               vrms_branches[mask_tag_branches=="Outside"])  # Fit outside mask points
popt_pin, pcov_pin = curve_fit(func,radius_pc_branches[mask_tag_branches=="Partly Inside"],\
                               vrms_branches[mask_tag_branches=="Partly Inside"])  # Fit to outside mask points
popt_in, pcov_in = curve_fit(func,radius_pc_branches[mask_tag_branches=="Inside"],\
                             vrms_branches[mask_tag_branches=="Inside"])  # Fit to outside mask points
perr_out = np.sqrt(np.diag(pcov_out))  # error to fits
perr_pin = np.sqrt(np.diag(pcov_pin))
perr_in = np.sqrt(np.diag(pcov_in))
print(perr_out, perr_in)


#
#         kde plots for the leaves
# 

plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))

sns.kdeplot(x=radius_pc_leaves[HII_tag_leaves=="Inside"],\
            y=vrms_leaves[HII_tag_leaves=="Inside"],\
            levels=4,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=2,\
            color='firebrick',\
            fill=0,\
            alpha=1)

sns.kdeplot(x=radius_pc_leaves[HII_tag_leaves=="Outside"],\
            y=vrms_leaves[HII_tag_leaves=="Outside"],\
            levels=4,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=2,\
            color='tab:blue',\
            fill=0,\
            alpha=1)



#             Scatter plot of branches
# (inside and outside the threshold mask separate)
# 

plt.scatter(radius_pc_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Inside")],\
            vrms_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Inside")],\
            color='firebrick',\
            marker='o',\
            s=10,\
            alpha=1.0)

plt.scatter(radius_pc_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Outside")],\
            vrms_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Outside")],\
            color='firebrick',\
            marker='o',\
            facecolor='none',\
            s=10,\
            alpha=1.0)

plt.scatter(radius_pc_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Inside")],\
            vrms_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Inside")],\
            color='tab:green',\
            marker='o',\
            s=10,\
            alpha=1.0)

plt.scatter(radius_pc_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Outside")],\
            vrms_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Outside")],\
            color='tab:green',\
            marker='o',\
            facecolor='none',\
            s=10,\
            alpha=1.0)

plt.scatter(radius_pc_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Inside")],\
            vrms_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Inside")],\
            color='tab:blue',\
            marker='o',\
            s=10,\
            alpha=1.0)

plt.scatter(radius_pc_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Outside")],\
            vrms_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Outside")],\
            color='tab:blue',\
            marker='o',\
            facecolor='none',\
            s=10,\
            alpha=1.0)


#             Scatter plot of trunks 
# (inside and outside the threshold mask separate)
# 

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Inside")],\
            vrms_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Inside")],\
            color='firebrick',\
            marker='*',\
            s=50)

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Outside")],\
            vrms_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Outside")],\
            color='firebrick',\
            marker='*',\
            facecolor='none',\
            s=50)

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Inside")],\
            vrms_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Inside")],\
            color='tab:green',\
            marker='*',\
            s=50)

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Outside")],\
            vrms_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Outside")],\
            color='tab:green',\
            marker='*',\
            facecolor='none',\
            s=50)

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Inside")],\
            vrms_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Inside")],\
            color='tab:blue',\
            marker='*',\
            s=50)

plt.scatter(radius_pc_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Outside")],\
            vrms_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Outside")],\
            color='tab:blue',\
            marker='*',\
            facecolor='none',\
            s=50)


#
#         Plotting the fitted power law function
#

xdata = np.logspace(-1,1.3,100)

plt.plot(xdata,\
         func(xdata,popt_out[0],popt_out[1]),\
         color='tab:blue',\
         label='Out:%5.3f x$^{%5.3f\pm%5.3f}$'\
               %(popt_out[0],popt_out[1],perr_out[1]))
                
#plt.plot(xdata,\
#         func(xdata,popt_pin[0],popt_pin[1]),\
#         color='tab:green',\
#         label='Partly In:%5.3f x$^{%5.3f\pm%5.3f}$'\
               %(popt_pin[0],popt_pin[1],perr_pin[1]))

plt.plot(xdata,\
         func(xdata,popt_in[0],popt_in[1]),\
         color='firebrick',\
         label='In:%5.3f x$^{%5.3f\pm%5.3f}$'\
               %(popt_in[0],popt_in[1],perr_in[1]))
#plt.xscale('log')
#plt.yscale('log')



#         Labels and Legends and Aesthetic Settings
# 

plt.xlabel(r'$R\,\,[pc]$')
plt.ylabel(r'$\sigma_v\,\,[km.s^{-1}]$')

plt.legend(framealpha=0.2,loc='upper left')

plt.tight_layout()

plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Scaling_Relation_V_R_vs_8um.pdf",format="pdf")









#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&%    Heyer et al. Plot
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


plt.rcParams.update({'font.size': 14})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,7))

#         kde plots for the leaves
#

sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="Inside")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="Inside")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='firebrick',\
            fill=0,\
            alpha=1,\
            label='Mostly Inside')
            
sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="Partly Inside")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="Partly Inside")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='tab:green',\
            fill=0,\
            alpha=1,\
            label='Partly Inside')
            
sns.kdeplot(x=sigma_leaves[(mask_tag_leaves=="Outside")&(sigma_leaves!=0)],\
            y=y_var_leaves[(mask_tag_leaves=="Outside")&(sigma_leaves!=0)],\
            levels=5,\
            log_scale=10,\
            ax=ax,\
            bw_adjust=1.5,\
            color='tab:blue',\
            fill=0,\
            alpha=1,\
            label='Outside')
            

#         Scatter plot of branches based on inside and outside the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(sigma_branches[(mask_tag_branches=="Inside")], y_var_branches[(mask_tag_branches=="Inside")],\
                         color='firebrick',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_branches[(mask_tag_branches=="Partly Inside")], y_var_branches[(mask_tag_branches=="Partly Inside")],\
                         color='tab:green',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_branches[(mask_tag_branches=="Outside")], y_var_branches[(mask_tag_branches=="Outside")],\
                         color='tab:blue',marker='o',s=10,alpha=0.7)        # for column density dendrogram


#         Scatter plot of trunks based on inside and outside the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(sigma_trunks[(mask_tag_trunks=="Inside")], y_var_trunks[(mask_tag_trunks=="Inside")],\
                         color='firebrick',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_trunks[(mask_tag_trunks=="Partly Inside")], y_var_trunks[(mask_tag_trunks=="Partly Inside")],\
                         color='tab:green',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(sigma_trunks[(mask_tag_trunks=="Outside")], y_var_trunks[(mask_tag_trunks=="Outside")],\
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
plt.xlim(left=1,right=10**4.5)
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
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_Dynamical_State_vs_8um.pdf", format="pdf",dpi=800)


##################################
# Heyer et al. Plot for Dust Mass
##################################
plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))

#         kde plots for the leaves
# --------------------------------------------------------------------------------------------------
sns.kdeplot(x=Sigma_dust_leaves[mask_tag_leaves=="Inside"],y=y_var_leaves[mask_tag_leaves=="Inside"],\
            levels=4,log_scale=10,ax=ax,bw_adjust=1.5,color='firebrick',fill=0,alpha=1)
#sns.kdeplot(x=Sigma_dust_leaves[mask_tag_leaves=="Outside"],y=y_var_leaves[mask_tag_leaves=="Outside"],\
#            levels=2,log_scale=10,ax=ax,bw_adjust=1.5,color='tab:blue',fill=0)
import scipy.stats as st
x = Sigma_dust_leaves[(mask_tag_leaves=="Outside")&(~np.isnan(Sigma_dust_leaves))]
y = y_var_leaves[(mask_tag_leaves=="Outside")&(~np.isnan(Sigma_dust_leaves))]
xx, yy = np.mgrid[0:400:100j, 0.1:3:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
plt.contour(xx, yy, f, levels=4, colors='tab:blue',alpha=1)

#         Scatter plot of branches based on inside and outside the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Inside")], y_var_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Inside")],\
                         color='firebrick',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Outside")], y_var_branches[(mask_tag_branches=="Inside")&(HII_tag_branches=="Outside")],\
                         color='firebrick',marker='o',facecolor='none',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Inside")], y_var_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Inside")],\
                         color='tab:green',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Outside")], y_var_branches[(mask_tag_branches=="Partly Inside")&(HII_tag_branches=="Outside")],\
                         color='tab:green',marker='o',facecolor='none',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Inside")], y_var_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Inside")],\
                         color='tab:blue',marker='o',s=10,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Outside")], y_var_branches[(mask_tag_branches=="Outside")&(HII_tag_branches=="Outside")],\
                         color='tab:blue',marker='o',facecolor='none',s=10,alpha=0.7)        # for column density dendrogram


#         Scatter plot of trunks based on inside and outside the threshold mask.
# --------------------------------------------------------------------------------------------------
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Inside")], y_var_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Inside")],\
                         color='firebrick',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Outside")], y_var_trunks[(mask_tag_trunks=="Inside")&(HII_tag_trunks=="Outside")],\
                         color='firebrick',marker='*',facecolor='none',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Inside")], y_var_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Inside")],\
                         color='tab:green',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Outside")], y_var_trunks[(mask_tag_trunks=="Partly Inside")&(HII_tag_trunks=="Outside")],\
                         color='tab:green',marker='*',facecolor='none',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Inside")], y_var_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Inside")],\
                         color='tab:blue',marker='*',s=50,alpha=0.7)        # for column density dendrogram
plt.scatter(Sigma_dust_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Outside")], y_var_trunks[(mask_tag_trunks=="Outside")&(HII_tag_trunks=="Outside")],\
                         color='tab:blue',marker='*',facecolor='none',s=50,alpha=0.7)        # for column density dendrogram


def func2(x,vir):
    return np.sqrt(np.pi*G*x*vir/5)

x2data = np.linspace(1.0,2000,5000)
y2data = np.linspace(0.3,7,5000)
X2, Y2 = np.meshgrid(x2data, y2data)
def balance (x,y,press):
    return y**2-((4/3)*press/x)-(np.pi*G*x/5)
plt.plot(x2data,func2(x2data,vir=1),'-',color='#014182',label=r'$\alpha_{vir}=1$')
plt.plot(x2data,func2(x2data,vir=2),'-.',color='#014182', label=r'$\alpha_{vir}=2$')

con1=plt.contour(X2, Y2, balance(X2,Y2,press=10), levels=[0], colors=['black'])
fmt = {}
strs = [r"$P_{ext}=10$"]
for l, s in zip(con1.levels, strs):
    fmt[l] = s
manual_location=[(2,1)]
plt.clabel(con1, fmt=fmt, inline=True, fontsize=7, manual=manual_location)
con2=plt.contour(X2, Y2, balance(X2,Y2,press=100), levels=[0], colors=['gray'])
fmt = {}
strs = [r'$P_{ext}=100$']
for l, s in zip(con2.levels, strs):
    fmt[l] = s
manual_location=[(15,3)]
plt.clabel(con2,fmt=fmt,inline=1,fontsize=7, manual=manual_location)
con3=plt.contour(X2, Y2, balance(X2,Y2,press=500), levels=[0], colors=['silver'])
fmt = {}
strs = [r"$P_{ext}=500$"]
for l, s in zip(con3.levels, strs):
    fmt[l] = s
manual_location=[(40,4)]
plt.clabel(con3, fmt=fmt, inline=True, fontsize=7, manual=manual_location)

plt.ylim(bottom=0.2,top=7)
plt.xlim(left=1,right=2000)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\sigma_v / R^{0.5} \,\,\, [km.s^{-1}.pc^{-0.5}]$')
plt.xlabel(r'$\Sigma \,\,\, [M_{\odot}.pc^{-2}]$')
plt.legend(loc='lower right',ncol=2,framealpha=0.5)
plt.show()
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_Dynamical_State_Dust_vs_8um.eps", format="eps",overwrite=True)



#        Plot the 8um Threshold Mask
# --------------------------------------------------------------------------------------------------

import matplotlib
from matplotlib.colors import LogNorm

#image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
#image_data = image_hdu.data
#image_wc = wcs.WCS(image_hdu.header) 
#plt.rcParams.update({"font.size":12})
#fig = plt.figure(figsize=(12,12))
#ax = fig.add_subplot(1,1,1,projection=image_wc)

#current_cmap = matplotlib.cm.get_cmap(plt.cm.YlOrBr)
#current_cmap.set_bad(color='white')
#im = ax.imshow(glm8_mask, origin='lower')
#ax.set_xlabel("Galactic Longitude ($^o$)")
#ax.set_ylabel("Galactic Latitude ($^o$)")



from matplotlib.colors import LogNorm
import matplotlib

#image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_12CO_mom_0.fits')[0]
image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/N_Dust-reproject.fits')[0]
image_data = image_hdu.data
#image_data = RC_mask_sum

image_wc = wcs.WCS(image_hdu.header) 
plt.rcParams.update({"font.size":12})
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1,projection=image_wc)

current_cmap = matplotlib.cm.get_cmap(plt.cm.Greys)
#current_cmap = matplotlib.cm.get_cmap(plt.cm.nipy_spectral)
current_cmap.set_bad(color='white')
#im=ax.imshow(image_data, origin='lower', interpolation='nearest',cmap=current_cmap, vmax=100,norm=LogNorm())
im=ax.imshow(image_data, origin='lower', interpolation='nearest',cmap=current_cmap, norm=LogNorm())
ax.set_xlabel("Galactic Longitude ($^o$)")
ax.set_ylabel("Galactic Latitude ($^o$)")
cax = plt.axes([0.9, 0.125, 0.015, 0.755])
cbar=plt.colorbar(im,cax=cax)
cbar.set_label('[K.km.s$^{-1}$]')
for RC_region in RC_regions_pix:
  ax.add_artist(RC_region.as_artist(facecolor='none', edgecolor='orange'))
#ax.contour(glm8_mask,levels=[1.0],linewidth=2,color='red')
#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_Masks.eps")

p = d.plotter()

counter = 0
plt.rcParams.update({"lines.linewidth":0.5})
trunk_id = cat[cat['Structure Tag']=="branch"]['_idx']
for s in d.all_structures:
  if s.idx in trunk_id:
#    print(s)    
    if cat[s.idx]['Mask_Tag'] == "Inside":
      p.plot_contour(ax, structure=s, colors='firebrick', linewidth=0.5,smooth=2)
    elif cat[s.idx]['Mask_Tag'] == "Partly Inside":
      p.plot_contour(ax, structure=s, colors='tab:green', linewidth=0.5,smooth=2)#  counter += 1 
    elif cat[s.idx]['Mask_Tag'] == "Outside":
      p.plot_contour(ax, structure=s, colors='tab:blue', linewidth=0.5,smooth=2)#  counter += 1 
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_branches.eps",overwrite=True)
