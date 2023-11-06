import numpy as np
import glob
import os
from tqdm import tqdm
from astropy.table import Table
from astropy.io import fits
import astropy

ver = 'v0'

outdir = '/n/holyscratch01/conroy_lab/vchandra/mage/data/reduced/%s/' % ver
plotdir1 = '/n/holyscratch01/conroy_lab/vchandra/mage/plots/%s/'% ver
plotdir2 = '/n/holyscratch01/conroy_lab/vchandra/mage/plots/%s/reduction/' % ver

################################
# COLLATE SPECTRA
################################

try:
    os.mkdir(outdir)
except:
    print('dir exists')

try:
    os.mkdir(plotdir1)
except:
    print('dir exists')

try:
    os.mkdir(plotdir2)
except:
    print('dir exists')

specfiles = glob.glob('/n/holyscratch01/conroy_lab/vchandra/mage/data/202*/reduced_%s/magellan_mage_A/Science/coadd/*.fits' % ver)

print('there are %i co-added spectra' % len(specfiles))

for file in (specfiles):
    name = file.split('/')[-1].split('_')[0].replace('.fits', '')
    date = file.split('/')[-6]

    if name == 'j2035m2245': # mis-named file in header
        name = 'j2035m2445'
    
    outname = date + '_' + name + '.fits'
    
    cmd = 'rsync -vzr %s %s%s' % (file, outdir, outname)
    
    os.system(cmd)

################################
# COLLATE PLOTS
################################

# plotfiles = glob.glob('/n/holyscratch01/conroy_lab/vchandra/mage/data/202*/reduced_%s/magellan_mage_A/plots/*.png' % ver)

# for file in (plotfiles):
#     name = file.split('/')[-1].replace('.png', '')
#     date = file.split('/')[-5]

#     if name == 'j2035m2245': # mis-named file in header
#         name = 'j2035m2445'
    
#     outname = date + '_' + name + '.fits'
    
#     cmd = 'rsync -vzr %s %s%s' % (file, plotdir2, outname)
    
#     os.system(cmd)

################################
# MAKE SPALL
################################

specfiles = glob.glob(outdir + '*.fits')

keys = ['RA', 'DEC', 'TARGET', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 
       'EXPTIME']

spall_dict = dict(

    name = [],
    date = [],
    specfile = [],

)

for ii,file in enumerate(specfiles):
    name = file.split('/')[-1].split('_')[-1].replace('.fits', '')
    date = '_'.join(file.split('/')[-1].split('_')[-4:-1])
    
    spall_dict['name'].append(name)
    spall_dict['date'].append(date)
    spall_dict['specfile'].append(file)
    
    with fits.open(file) as f:
        
        if ii == 0:
            for key in keys:
                spall_dict['mage_' + key.lower()] = [];
        
        for key in keys:
            spall_dict['mage_' + key.lower()].append(f[0].header[key])


spall = Table(spall_dict)
print('there are %i observations in spall, for %i unique targets...' % (len(spall), len(np.unique(spall['name']))))

################################
# MATCH TO TARGETDB
################################

tdb = Table.read('/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/tdb/targetdb_2022b.fits')
tdb23a = Table.read('/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/tdb/targetdb_2023a.fits')
tdb23b = Table.read('/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/tdb/targetdb_2023b.fits')


tdb = astropy.table.unique(astropy.table.vstack((tdb, tdb23a, tdb23b)), keys = 'name')

for key in list(tdb.columns):
    if key == 'name':
        continue
    
    tdb.rename_column(key, 'tdb_' + key)

spall_tdb = astropy.table.join(spall, tdb, keys = 'name', join_type = 'left')
isnan = spall_tdb['tdb_ra'].mask

print('there are %i rows after matching to targetDB...' % len(spall_tdb))
print('there are %i NaN coordinates, removing from spall:' % np.sum(isnan))
print(list(spall_tdb[isnan]['name']))
spall_tdb= spall_tdb[~isnan]
print('there are %i rows after removing nans...' % len(spall_tdb))

spall_tdb.write('/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/spall.fits', overwrite = True)