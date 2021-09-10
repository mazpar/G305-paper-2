from astropy.table import Table
from astrodendro import Dendrogram, ppv_catalog, structure
from astropy import units as u
from astropy import wcs
#from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astrodendro.pruning import all_true, min_vchan, min_delta, min_area
from astropy import constants as const
from spectral_cube import SpectralCube
#import radio_beam
#import aplpy
import seaborn as sns
#import scipy.stats as sst

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

data = masked_cube.hdu.data
hd = masked_cube.hdu.header
wc = wcs.WCS(hd)
data_cd = cube_cd.hdu.data
hd_cd = cube_cd.hdu.header
wc_cd = wcs.WCS(hd_cd)

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
#&%
#&% Removing the 4th dimension
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&


#def remove_hd_dim(hd, dim=3):

#    """
 #   Keep information from datacube header
 #   only related to dimension to keep
 #   
 #   hd: astropy.header instance
 #   
 #   dim: int
 #   dimensions to keep, default=3
    
 #  """
    #if dim == 2:
      #  nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS2'],hd['NAXIS1']])).header    
     #   dimlist = ['1','2']
    #if dim == 3:
      #  nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS3'],hd['NAXIS2'],hd['NAXIS1']])).header
     #   dimlist = ['1','2','3']        
    #    
   # for i in dimlist:
  #      for t in ['CRVAL','CRPIX','CDELT','CTYPE','CROTA','CUNIT']:
 #           if hd.get(t+i) != None:
#                nhd[t+i] = hd[t+i]

#    for t in ['BUNIT','BMAJ','BMIN','BPA','RESTFRQ']:
#        if hd.get(t) != None:
#            nhd[t] = hd[t]

#    return nhd

#data = data.squeeze()  # remove 4th dimension
#hd = remove_hd_dim(hd) # change the header
#hdu.data = data
#hdu.header = hd
#hdu.writeto('13CO_squeeze.fits', overwrite=1)


##    Custom Definitions for the Dendrogram    ##
rms = 0.12 # rms noise
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


is_independent = all_true((min_delta(5*rms), min_area(1*ppb), min_vchan(6)))

d = Dendrogram.compute(data, min_value=5*rms, wcs=wc, is_independent = is_independent, verbose=1)
d.save_to('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_dendro.fits')





#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% START HERE IF DENDROGRAM ALREADY RUN
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

d = Dendrogram.load_from('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_dendro.fits')

leaf_id=np.array([])
for s in d.all_structures:
    if s.is_leaf:
        leaf_id=np.append(leaf_id,s.idx)
print(leaf_id)
leaf_id.sort()




#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Viewing and Plotting the Dendrogram
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#v = d.viewer()
#v.show()
p = d.plotter()


fig,ax = plt.subplots(nrows = 2, ncols=1)
plt.rcParams.update({"font.size":6})

p.plot_tree(ax[0],color='seagreen',lw=0.5)

ax[0].set_ylabel('$^{13}$CO Peak Intensity')
ax[0].set_xlabel('Index of Structure')

plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Dendrogram_G305.eps",dpi=300)





############################
##    Read the 8um Map    ##
############################

hdu_glm8 = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
data_glm8 = hdu_glm8.data
w_glm8 = wcs.WCS(hdu_glm8.header)


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Creating the ppv_catalog
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

metadata = {}
metadata['data_unit'] = u.Jy
#metadata['data_unit'] = u.cm**-2
#metadata ['spatial_scale'] = 9 * u.arcsecond # Pixel Size
metadata['beam_major'] = (bmaj * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['beam_minor'] = (bmin * u.deg .to(u.arcsecond))*u.arcsecond # FWHM
metadata['velocity_scale'] = 0.5 * u.km/u.s # v_res

cat = ppv_catalog(d,metadata)




##
##   Additions to Integrated Intensity Catalog
##

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



##    adding location of clumps
cat['x_c'] = np.zeros(len(cat),dtype=float)
cat['x_c'] = hd['crval1']+((cat['x_cen']-hd['crpix1'])*hd['CDELT1'])
cat['x_c'].unit = u.degree

cat['y_c'] = np.zeros(len(cat),dtype=float)
cat['y_c'] = hd['crval2']+((cat['y_cen']-hd['crpix2'])*hd['CDELT2'])
cat['y_c'].unit = u.degree


##    adding a radius column to the catalog
cat['radius_pc'] = np.zeros(len(cat),dtype=float)
cat['radius_pc'] = eta*np.sqrt((sigma_majs*deltax_pc)*(sigma_mins*deltay_pc))
cat['radius_pc'].unit = u.parsec


##    adding a luminosity column to the catalog
cat['luminosity']=np.zeros(len(cat),dtype=float)
cat['luminosity'] = cat['flux']*deltav_kms*deltax_pc*deltay_pc
cat['luminosity'].unit = u.K * u.km / u.s * u.pc * u.pc


##    adding a luminosity mass column to the catalog
#cat['M_lum']=np.zeros(len(cat),dtype=float)
#cat['M_lum'] = cat['luminosity']*4.4*x2
#cat['M_lum'].unit = u.solMass

##   adding a mass column
cat['Mass']=np.zeros(len(cat),dtype=float)
cat['Mass'].unit = u.solMass
for s in d.all_structures:
  str_mask = s.get_mask()
  str_cube = data_cd[str_mask]
  total_cd = np.nansum(str_cube)
  mass = mu*mp.value*R13_inv*total_cd*deltax_pc*deltay_pc*(u.parsec.to(u.centimeter))**2
  cat[s.idx]['Mass'] = mass


##    calculating virial parameter alpha
cat['virial_parameter'] = np.zeros(len(cat),dtype=float)
cat['virial_parameter'] = (3*((sigma_v)**2)*cat['radius_pc'])/(cat['Mass']*G)
cat['virial_parameter'].unit = ''


##    adding a surface density column to the catalog
cat['Sigma_exact']=np.zeros(len(cat),dtype=float)
cat['Sigma_exact'] = cat['Mass']/(cat['area_exact']*deltax_pc*deltay_pc)
cat['Sigma_exact'].unit = u.solMass/(u.pc*u.pc)

cat['Sigma_ellipse'] = np.zeros(len(cat),dtype=float)
cat['Sigma_ellipse'] = cat['Mass']/(np.pi*cat['radius_pc']**2)
cat['Sigma_ellipse'].unit = u.solMass/(u.pc*u.pc)


##    adding a surface density column to the catalog
cat['n_H2'] = np.zeros(len(cat),dtype=float)
cat['n_H2'] = (cat['Mass']/(mu*mp*(4/3)*np.pi*cat['radius_pc']*cat['radius_pc']*cat['radius_pc'])).to(1/u.cm**3)
cat['n_H2'].unit = 1/(u.cm**3)



##  Structure tags 

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




##    Calculate 8um flux over masks
cat['8um Flux Mean']=np.zeros(len(cat),dtype=float)

for l in d.all_structures:    ## adding average 8um flux of leaves
  leaf_mask_3d = l.get_mask()
  leaf_mask = np.max(leaf_mask_3d,axis=0)
  masked_8um_map = data_glm8[leaf_mask]
  mean_8um_flux = np.nanmean(masked_8um_map)
  cat[l.idx]['8um Flux Mean']=mean_8um_flux



##    Save Catalogs   ##
cat.write('/home/pmazumdar/Documents/LASMA/Dendrogram/fits_files/G305_13CO_cat.fits',overwrite=True)


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
mask_tag = np.array([])
glm8_flux = np.array([])


for i in leaf_id:
    L_CO = np.append(L_CO, cat[int(i)]['luminosity'])
    radius_pc = np.append(radius_pc,cat[int(i)]['radius_pc'])
    vrms = np.append(vrms,sigma_v[int(i)])
    mass = np.append(mass, cat[int(i)]['M_lum'])
    #mass = np.append(mass, cat[i]['Mass'])
    Sigma_exact = np.append(Sigma_exact, cat[int(i)]['Sigma_exact'])
    Sigma_ellipse = np.append(Sigma_ellipse,cat[int(i)]['Sigma_ellipse'])
    alpha = np.append(alpha,cat[int(i)]['virial_parameter'])
    mask_tag = np.append(mask_tag,cat[int(i)]['mask_ID'])
    glm8_flux =  np.append(glm8_flux,cat[int(i)]['8um Flux Mean'])


#
# Variables for plotting all structures
#

L_CO = np.array(cat['luminosity'])
radius_pc = np.array(cat['radius_pc'])
vrms = np.array(cat['v_rms']*deltav_kms)
mass = np.array(cat['M_lum'])
Sigma_exact = np.array(cat['Sigma_exact'])
Sigma_ellipse = np.array(cat['Sigma_ellipse'])
alpha = np.array(cat['virial_parameter'])
glm8_flux =np.array(cat['8um Flux Mean'])


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Plotting Individual Properties
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

#from astropy.visualization import hist
plt.rcParams.update({'font.size':10})
fig,ax = plt.subplots(ncols=3,nrows=1,figsize=(12,4))
#Mlum_low_hist = hist(np.log10(mass[mask_tag==0]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[0], color='orange',label='Lum_out')
#Mlum_mid_hist = hist(np.log10(mass[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[0], color='gray',label='Lum_mid')
#Mlum_high_hist = hist(np.log10(mass[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[0], color='green',label='Lum_in')
#Mcd_hist_low = hist(np.log10(mass_cd[mask_tag_cd==0]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed', density=True, ax=ax[0], color='orange', label='Ncd_low')
#Mcd_hist_mid = hist(np.log10(mass_cd[mask_tag_cd==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, linestyle='dashed', density=True, ax=ax[0], color='gray', label='Ncd_mid')
#Mcd_hist_high = hist(np.log10(mass_cd[mask_tag_cd==1]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed', density=True, ax=ax[0], color='green', label='Ncd_high')
#sns.kdeplot(data=np.log10(mass_cd),hue=mask_tag_cd,ax=ax[0],fill=1,alpha=0.25,label='Outside_Mask')
#sns.kdeplot(x=np.log10(mass),y=np.log10(glm8_flux),ax=ax[0])
#sns.kdeplot(x=np.log10(mass_cd),y=np.log10(glm8_flux_cd),ax=ax[0],levels=6,fill=1,color='tab:blue')
sns.kdeplot(y=np.log10(Sigma_ellipse),x=np.log10(glm8_flux),ax=ax[0],bw_adjust=2,fill=1,color='tab:blue')
#ks_stat_mass = sst.ks_2samp(mass_cd[mask_tag_cd==0],mass_cd[mask_tag_cd==1])

#ax[0].set_xlabel(r'log$_{10}$(M/M$_{\odot}$)')
ax[0].set_ylabel(r'log$_{10}\Sigma (M_{\odot}.pc^{-2})$')
ax[0].set_xlabel(r'log$_{10}$ 8$\mu$m Mean Flux (MJy/sr)')
#ax[0].legend()

#Req_hist_low = hist(np.log10(radius_pc[mask_tag==0]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[1], color='orange')
#Req_hist_mid = hist(np.log10(radius_pc[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[1], color='gray')
#Req_hist_high = hist(np.log10(radius_pc[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[1], color='green')
#Req_cd_hist_low =  hist(np.log10(radius_pc_cd[mask_tag_cd==0]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed',density=True, ax=ax[1], color='orange', label='Ncd_out')
#Req_cd_hist_mid =  hist(np.log10(radius_pc_cd[mask_tag_cd==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, linestyle='dashed',density=True, ax=ax[1], color='gray', label='Ncd_mid')
#Req_cd_hist_high =  hist(np.log10(radius_pc_cd[mask_tag_cd==1]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed',density=True, ax=ax[1], color='green', label='Ncd_in')
#sns.kdeplot(data=np.log10(radius_pc_cd),hue=mask_tag_cd,fill=1,alpha=0.25,ax=ax[1])
#sns.kdeplot(x=np.log10(radius_pc_cd),y=np.log10(glm8_flux_cd),ax=ax[1],levels=6,fill=1,color='tab:orange')
#ks_stat_Req = sst.ks_2samp(radius_pc_cd[mask_tag_cd==0],radius_pc_cd[mask_tag_cd==1])

#ax[1].set_xlabel(r'log$_{10}$(R$_{eq}$/pc)')
#ax[1].set_ylabel('')

#vir_hist_low = hist(np.log10(alpha[mask_tag==0]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[2], color='orange')
#vir_hist_mid = hist(np.log10(alpha[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[2], color='gray')
#vir_hist_high = hist(np.log10(alpha[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[2], color='green')
#vir_cd_hist_low =  hist(np.log10(alpha_cd[mask_tag_cd==0]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed',density=True, ax=ax[2], color='orange')
#vir_cd_hist_mid =  hist(np.log10(alpha_cd[mask_tag_cd==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, linestyle='dashed',density=True, ax=ax[2], color='gray')
#vir_cd_hist_high =  hist(np.log10(alpha_cd[mask_tag_cd==1]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed',density=True, ax=ax[2], color='green')
#sns.kdeplot(data=np.log10(alpha_cd),hue=mask_tag_cd,fill=1,alpha=0.25,ax=ax[2])
#sns.kdeplot(x=np.log10(alpha_cd),y=np.log10(glm8_flux_cd),ax=ax[2],levels=6,fill=1, color='yellow')
sns.kdeplot(y=np.log10(alpha),x=np.log10(glm8_flux),ax=ax[1],bw_adjust=2,fill=1, color='yellow')
#ks_stat_vir = sst.ks_2samp(alpha_cd[mask_tag_cd==0],alpha_cd[mask_tag_cd==1])

ax[1].set_ylabel(r'log$_{10}(\alpha_{vir})$')
ax[1].set_xlabel('')

#vrms_hist_low = hist(np.log10(vrms[mask_tag==0]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[3], color='orange')
#vrms_hist_mid = hist(np.log10(vrms[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[3], color='gray')
#vrms_hist_high = hist(np.log10(vrms[mask_tag==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, density=True, ax=ax[3], color='green')
#vrms_cd_hist_low =  hist(np.log10(vrms_cd[mask_tag_cd==0]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed', density=True, ax=ax[3], color='orange')
#vrms_cd_hist_mid =  hist(np.log10(vrms_cd[mask_tag_cd==1]), bins='knuth', histtype='step', alpha=1, linewidth=2, linestyle='dashed', density=True, ax=ax[3], color='gray')
#vrms_cd_hist_high =  hist(np.log10(vrms_cd[mask_tag_cd==1]), bins='knuth', histtype='stepfilled', alpha=0.4, linewidth=2, linestyle='dashed', density=True, ax=ax[3], color='green')
#sns.kdeplot(data=np.log10(vrms_cd,hue=mask_tag_cd,fill=1,alpha=0.25,ax=ax[3]))
sns.kdeplot(y=np.log10(vrms),x=np.log10(glm8_flux),ax=ax[2],bw_adjust=2,fill=1,color='red')
#ks_stat_vrms = sst.ks_2samp(vrms_cd[mask_tag_cd==0],vrms_cd[mask_tag_cd==1])

ax[2].set_ylabel(r'log$_{10}\,v_{rms}$')
ax[2].set_xlabel('')
plt.tight_layout()
#plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/histo_clump_props.png",format="png", overwrite=True)
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/all_str_props_vs_8um.png",format="png", dvi=300,overwrite=True)

#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Scaling relation tests
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))

#plt.scatter(radius_pc[mask_tag==0], vrms[mask_tag==0],color='orange',marker='x',s=20,label='Lum_Out')             # for intensity dendrogram
#plt.scatter(radius_pc[mask_tag==1], vrms[mask_tag==1],color='orange',marker='o',facecolor='none',s=10,label='Lum_In')             # for intensity dendrogram
plt.scatter(radius_pc_cd[mask_tag_cd==0], vrms_cd[mask_tag_cd==0],color='k',marker='x',s=20, label='Ncd_Out')        # for column density dendrogram
plt.scatter(radius_pc_cd[mask_tag_cd==1], vrms_cd[mask_tag_cd==1],color='k',marker='o',facecolor='none',s=10, label='Ncd_In')        # for column density dendrogram
#sns.kdeplot(x=radius_pc_cd,y=vrms_cd,hue=mask_tag_cd,levels=7,log_scale=10,ax=ax)
sns.kdeplot(x=radius_pc_cd,y=vrms_cd,levels=7,log_scale=10,ax=ax,color='k')
#sns.kdeplot(x=radius_pc_cd[mask_tag_cd==1],y=vrms_cd[mask_tag_cd==1],ax=ax)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0.1,2)
#plt.ylim(0.1,2)
from scipy.optimize import curve_fit
def func(x,a,b):
    return a*pow(x,b)
#popt, pcov = curve_fit(func,radius_pc,vrms)              # for intensity dendrogram
popt_cd_out, pcov_cd_out = curve_fit(func,radius_pc_cd[mask_tag_cd==0],vrms_cd[mask_tag_cd==0])         # for column density dendrogram
popt_cd_in, pcov_cd_in = curve_fit(func,radius_pc_cd[mask_tag_cd==1],vrms_cd[mask_tag_cd==1])         # for column density dendrogram
#perr = np.sqrt(np.diag(pcov))
perr_cd_out = np.sqrt(np.diag(pcov_cd_out))
perr_cd_in = np.sqrt(np.diag(pcov_cd_in))
#print(perr)
print(perr_cd_out, perr_cd_in)
xdata = np.logspace(-1,0.3,100)
#plt.plot(xdata,func(xdata,popt[0],popt[1]), 'k-', label='fit:a=%5.3f, b=%5.3f' %tuple(popt))
plt.plot(xdata,func(xdata,popt_cd_out[0],popt_cd_out[1]), color='tab:blue', label='outside:%5.3f x$^{%5.3f}$' %tuple(popt_cd_out))
plt.plot(xdata,func(xdata,popt_cd_in[0],popt_cd_in[1]), color='tab:orange', label='inside:%5.3f x$^{%5.3f}$' %tuple(popt_cd_in))
plt.xlabel(r'$R\,\,[pc]$')
plt.ylabel(r'$\sigma_v\,\,[km.s^{-1}]$')
plt.legend(framealpha=0.5)
plt.tight_layout()
plt.show()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Scaling_Relation_V_R_cd.eps",dpi=300,format="eps")

plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(5,5))
#plt.scatter(radius_pc[mask_tag==0],L_CO[mask_tag==0], color='k', marker='x',s=20,label='Outside_Mask')
#plt.scatter(radius_pc[mask_tag==1],L_CO[mask_tag==1], color='k', marker='o',s=10,facecolor='none',label='Inside Mask')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0.1,2)
#plt.ylim(0.01,500)
sns.kdeplot(x=radius_pc,y=L_CO,hue=mask_tag,levels=7,log_scale=10,ax=ax)

def func1(x,a,b):
    return a*pow((x*4.4*5),b)
popt_out, pcov_out = curve_fit(func1,radius_pc[mask_tag==0],L_CO[mask_tag==0])
perr_out = np.sqrt(np.diag(pcov_out))
popt_in, pcov_in = curve_fit(func1,radius_pc[mask_tag==1],L_CO[mask_tag==1])
perr_in = np.sqrt(np.diag(pcov_in))
print(perr_out, perr_in)
plt.plot(xdata,func1(xdata,popt_out[0],popt_out[1]), color='tab:blue', label='outside:%5.3f x$^{%5.3f}$' %tuple(popt_out))
plt.plot(xdata,func1(xdata,popt_in[0],popt_in[1]), color='tab:orange', label='inside:%5.3f x$^{%5.3f}$' %tuple(popt_in))
plt.xlabel(r'$R\,\,[pc]$')
plt.ylabel(r'$L_{^{13}CO}\,\,[K.km.s^{-1}]$')
plt.legend()
plt.show()
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/Scaling_Relation_L_R_leaf.eps",dpi=300,format="eps")

#
# Heyer et al. Plot
#
plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
def func2(x,vir):
    return np.sqrt(np.pi*G*x*vir/5)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(20,10000)
#plt.ylim(0.05,6)
#plt.scatter(Sigma_exact,vrms/np.sqrt(radius_pc),c='red',label=r'$\Sigma_{exact}$')
#plt.scatter(Sigma_exact_cd,vrms_cd/np.sqrt(radius_pc_cd),c='red',label=r'$\Sigma_{exact}$')
#plt.scatter(Sigma_ellipse_cd[mask_tag_cd==0],(vrms_cd/np.sqrt(radius_pc_cd))[mask_tag_cd==0],s=20,marker='o',c='tab:blue')   # for col. density dendrogram
#plt.scatter(Sigma_ellipse_cd[mask_tag_cd==1],(vrms_cd/np.sqrt(radius_pc_cd))[mask_tag_cd==1],s=20,marker='o',c='tab:orange')   # for col. density dendrogram
#plt.scatter(Sigma_ellipse[mask_tag==0],(vrms/np.sqrt(radius_pc))[mask_tag==0],s=20,marker='+',c='tab:blue',label=r'$\Sigma_{Intensity}$ Outside Mask')           # for intensity dendrogram
#plt.scatter(Sigma_ellipse[mask_tag==1],(vrms/np.sqrt(radius_pc))[mask_tag==1],s=10,marker='o',c='tab:blue',facecolor='none',label=r'$\Sigma_{Intensity}$ Inside Mask')           # for intensity dendrogram
#sns.kdeplot(x=Sigma_ellipse_cd,y=vrms_cd/np.sqrt(radius_pc_cd),hue=mask_tag_cd,log_scale=10,levels=7,legend=0,label=['outside','inside'])
cmap=plt.nipy_spectral()
plt.scatter(Sigma_ellipse_cd,(vrms_cd/np.sqrt(radius_pc_cd)),s=10,marker='o',c=np.log10(glm8_flux_cd),cmap=plt.cm.Greens)   # for col. density dendrogram
cbar = plt.colorbar(aspect=50,pad=0.01)
cbar.set_label(r'log(Mean 8$\mu$m Flux)')
sns.kdeplot(x=Sigma_ellipse_cd,y=vrms_cd/np.sqrt(radius_pc_cd),log_scale=10,levels=7,color='red',fill=1,alpha=0.25)

x2data = np.linspace(30,10000,5000)
y2data = np.linspace(0.08,7,5000)
#x2data = np.linspace(0.8,60,500)
#y2data = np.linspace(0.1,4,500)
X2, Y2 = np.meshgrid(x2data, y2data)
def balance (x,y,press):
    return y**2-((4/3)*press/x)-(np.pi*G*x/5)
plt.plot(x2data,func2(x2data,vir=1),'-',color='#014182',label=r'$\alpha_{vir}=1$')
plt.plot(x2data,func2(x2data,vir=2),'-.',color='#014182', label=r'$\alpha_{vir}=2$')
#plt.contour(X2,Y2,press,levels=1, color='k')

con1=plt.contour(X2, Y2, balance(X2,Y2,press=10), levels=[0], colors=['black'])
fmt = {}
strs = [r"$P_{ext}=10$"]
for l, s in zip(con1.levels, strs):
    fmt[l] = s
manual_location=[(40,0.5)]
plt.clabel(con1, fmt=fmt, inline=True, fontsize=7, manual=manual_location)
con2=plt.contour(X2, Y2, balance(X2,Y2,press=100), levels=[0], colors=['gray'])
fmt = {}
strs = [r'$P_{ext}=100$']
for l, s in zip(con2.levels, strs):
    fmt[l] = s
manual_location=[(60,2)]
plt.clabel(con2,fmt=fmt,inline=1,fontsize=7, manual=manual_location)
con3=plt.contour(X2, Y2, balance(X2,Y2,press=500), levels=[0], colors=['silver'])
fmt = {}
strs = [r"$P_{ext}=500$"]
for l, s in zip(con3.levels, strs):
    fmt[l] = s
manual_location=[(100,3)]
plt.clabel(con3, fmt=fmt, inline=True, fontsize=7, manual=manual_location)

plt.ylabel(r'$\sigma_v / R^{0.5} \,\,\, [km.s^{-1}.pc^{-0.5}]$')
plt.xlabel(r'$\Sigma \,\,\, [M_{\odot}.pc^{-2}]$')
plt.legend(loc='lower right',ncol=2,framealpha=0.5)
plt.show()
plt.tight_layout()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_Dynamical_State_leaf.png", format="png",overwrite=True)


#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% Plotting virialized structures
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
from matplotlib.colors import LogNorm
import matplotlib

image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/G305_13CO_mom_0.fits')[0]
#image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/8um-reproject.fits')[0]
image_data = image_hdu.data
image_wc = wcs.WCS(image_hdu.header) 
plt.rcParams.update({"font.size":12})
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1,projection=image_wc)

current_cmap = matplotlib.cm.get_cmap(plt.cm.Blues)
#current_cmap = matplotlib.cm.get_cmap(plt.cm.nipy_spectral)
current_cmap.set_bad(color='white')
#im=ax.imshow(image_data, origin='lower', interpolation='nearest',cmap=current_cmap, vmin=40, vmax=700,norm=LogNorm())
im=ax.imshow(image_data, origin='lower', interpolation='nearest',cmap=current_cmap)
ax.set_xlabel("Galactic Longitude ($^o$)")
ax.set_ylabel("Galactic Latitude ($^o$)")
cax = plt.axes([0.9, 0.125, 0.015, 0.755])
cbar=plt.colorbar(im,cax=cax)
cbar.set_label('[K.km.s$^{-1}$]')
#cbar.set_label('[MJy.sr$^{-1}$]')

counter = 0
plt.rcParams.update({"lines.linewidth":0.5})

for s in d.all_structures:
    if s.is_leaf:
        print(s.idx)    
        #if cat[s.idx]['virial_parameter']<=2:
        p.plot_contour(ax, structure=s, colors='black', linewidth=0.8)
#        else:
#          p.plot_contour(ax, structure=s, colors='red', linewidth=0.8)
#        par = s.parent
#        if par is not None:
#            p.plot_contour(ax, structure=par, colors='orange', linewidth=0.5)
        counter += 1 

import matplotlib.markers as markers
xhii = np.array([305.353,305.195,305.270,305.320,305.254,305.348,305.551,305.532])
yhii = np.array([0.193,0.033,-0.007,0.070,0.204,0.223,0.014,0.348])
hii_marker = markers.MarkerStyle(marker=(3,0,0),fillstyle="none")
#ax.scatter(xhii,yhii,c='k',marker=hii_marker,s=70,transform=ax.get_transform('galactic'),label='HII')

xuchii = np.array([305.362,305.368,305.562,305.55,305.200])
yuchii = np.array([0.150,0.213,0.013,-0.012,0.019])
uchii_marker = markers.MarkerStyle(marker=(3,0,0),fillstyle="none")
#ax.scatter(xuchii,yuchii,c='k',marker=uchii_marker,s=50,transform=ax.get_transform('galactic'),label='UC HII')

xbrc = np.array([305.244])
ybrc = np.array([0.224])
brc_marker = markers.MarkerStyle(marker=(4,1,0),fillstyle="none")
#ax.scatter(xbrc,ybrc,c='k',marker=brc_marker,s=50,transform=ax.get_transform('galactic'),label='BRC')

xh2o = np.array([305.22,305.21,305.41,305.35,305.35,305.33,305.09,305.13,305.20,305.26,305.72,305.89,305.75,305.81,305.80])
yh2o = np.array([0.28,0.21,0.25,0.20,0.15,0.07,0.1,0.08,0.00,-0.01,0.09,0.03,-0.08,-0.11,-0.24])
h2o_marker = markers.MarkerStyle(marker=(4,0,0), fillstyle="none")
#ax.scatter(xh2o,yh2o,c='purple',marker=h2o_marker,s=50,transform=ax.get_transform('galactic'),label='H$_2$O Maser')


xmmb = np.array([305.199,305.200,305.202,305.208,305.248,305.362,305.366,305.475,305.563,305.615,305.799,305.822,305.887])
ymmb = np.array([0.005,0.019,0.208,0.206,0.245,0.150,0.184,0.096,0.013,-0.344,-0.245,-0.115,0.017])
mmb_marker = markers.MarkerStyle(marker=(4,0,45), fillstyle="none")
#ax.scatter(xmmb,ymmb,c='orange',marker=mmb_marker,s=50,transform=ax.get_transform('galactic'),label='MMB Maser')
#ax.legend(ncol=2)

#for s in d.leaves:
#    if cat[s.idx]['virial_parameter']<=1:
        #mask = mask | s.get_mask()
 #       print(s.idx)        
#        p.plot_contour(ax, structure=s.idx, colors='red', linewidth=0.5)
#        counter += 1 
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_13CO_dendrogram_leaves_on_8um_map.eps",format="eps", overwrite=True)


# Column Density Plot
image_hdu = fits.open('/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/ntotal_sum_new.fits')[0]
image_data = image_hdu.data
image_wc = wcs.WCS(image_hdu.header) 
plt.rcParams.update({"font.size":12})
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(1,1,1,projection=image_wc)

current_cmap = matplotlib.cm.get_cmap(plt.cm.YlOrBr)
current_cmap.set_bad(color='white')
im = ax.imshow(image_data, origin='lower', interpolation='nearest',cmap=current_cmap,vmin=8e15,vmax=7e16,norm=LogNorm())
ax.set_xlabel("Galactic Longitude ($^o$)")
ax.set_ylabel("Galactic Latitude ($^o$)")
cax = plt.axes([0.9, 0.125, 0.015, 0.755])
cbar = plt.colorbar(im,cax=cax)
cbar.set_label('[cm$^{-2}$]')
counter = 0
plt.rcParams.update({"lines.linewidth":0.8})

for s in d.all_structures:
    if s.is_leaf:
        print(s.idx)        
        if cat[s.idx]['virial_parameter']<=2:
          p.plot_contour(ax, structure=s, colors='black', linewidth=0.5)
#        else:
#          p.plot_contour(ax, structure=s, colors='red', linewidth=0.5)
#        par = s.parent
#        if par is not None:
#            p.plot_contour(ax, structure=par, colors='orange', linewidth=0.5)
        counter += 1 

ax.scatter(xhii,yhii,c='k',marker=hii_marker,s=40,transform=ax.get_transform('galactic'),label='HII')
ax.scatter(xuchii,yuchii,c='k',marker=uchii_marker,s=70,facecolors='none',transform=ax.get_transform('galactic'),label='UC HII')
ax.scatter(xbrc,ybrc,c='k',marker=brc_marker,s=60,transform=ax.get_transform('galactic'),label='BRC')
ax.scatter(xh2o,yh2o,c='royalblue',marker=h2o_marker,s=40,transform=ax.get_transform('galactic'),label='H$_2$O Maser')
ax.scatter(xmmb,ymmb,c='lime',marker=mmb_marker,s=40,transform=ax.get_transform('galactic'),label='MMB Maser')
ax.legend()
plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/G305_ncd_dendrogram_leaves.eps",format="eps", overwrite=True)






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
sns.kdeplot(data=np.log10(mass_cd),ax=ax[0],c='gray',fill=0,alpha=0.5,label='N$_{cd}$')
sns.kdeplot(data=np.log10(cmp_mass[(cmp_flag==3)]),ax=ax[0],c='blue',fill=1,alpha=0.5,label='CHIMPS')
sns.kdeplot(data=AGAL_mass,ax=ax[0],c='maroon',fill=0,alpha=0.5,linestyle='dashed',label='AGAL')
ax[0].set_xlabel(r'log$_{10}$(M/M$_{\odot}$)')
ax[0].legend(loc='upper left')

sns.kdeplot(data=np.log10(radius_pc),ax=ax[1],c='goldenrod',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(radius_pc_cd),ax=ax[1],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_radius[(cmp_flag==3)]),ax=ax[1],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_radius[AGAL_radius>0]),ax=ax[1],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[1].set_xlabel(r'log$_{10}$(R$_{eq}$/pc)')
ax[1].set_ylabel('')

sns.kdeplot(data=np.log10(alpha),ax=ax[2],c='goldenrod',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(alpha_cd),ax=ax[2],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_alpha[(cmp_flag==3)]),ax=ax[2],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_alpha[AGAL_alpha>0]),ax=ax[2],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[2].set_xlabel(r'log$_{10}(\alpha_{vir})$')
ax[2].set_ylabel('')

sns.kdeplot(data=np.log10(vrms),ax=ax[3],c='goldenrod',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(vrms_cd),ax=ax[3],c='gray',fill=0,alpha=0.5)
sns.kdeplot(data=np.log10(cmp_vrms[(cmp_flag==3)]),ax=ax[3],c='blue',fill=1,alpha=0.5)
sns.kdeplot(data=np.log10(AGAL_vrms[AGAL_vrms>0]),ax=ax[3],c='maroon',fill=0,alpha=0.5,linestyle='dashed')
ax[3].set_xlabel(r'log$_{10}\,v_{rms}$')
ax[3].set_ylabel('')

plt.savefig("/home/pmazumdar/Documents/LASMA/Reduction/class_maps/temp/plots/clump_properties_comparison.png",format="png", overwrite=True)






#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
#&%
#&% ATLASGAL CSC Data Load
#&%
#&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


