# xmatches supplied clean KTABLE to all-band photometry on Cannon

import numpy as np
from astropy.io import ascii
from astropy.table import Table
import glob
import astropy
import os

print('starting script')

infile = '/n/holyscratch01/conroy_lab/vchandra/mage/MS/ms_input/source_ids.txt'
outfolder = '/n/holyscratch01/conroy_lab/vchandra/mage/MS/sams_input/'

source_ids = np.loadtxt(infile)

# XMATCH Charlie's Gaia EDR3 photometric tables to my Gaia giants sample
# Loop over Charlie's files, match my file, save as h5

charlie_cats = '/n/holystore01/LABS/conroy_lab/Lab/gaia/edr3/gall2/catalogs/'
charlie_cats = glob.glob(charlie_cats + '*')

#charlie_cats = charlie_cats[-1:] # for testing, REMOVE

#print(charlie_cats)

tab = Table();

# try:
#     os.mkdir(outfolder)
# except:
#     print('output folder exists')

for bigcat in charlie_cats:
    print('-----------')
    print('Big Catalog: %s' % bigcat.split('/')[-1])
    
    bs = bigcat.split('/')[-1].split('b')[-1].split('_')
    b1 = float(bs[0])
    b2 = float(bs[1][:5])
    bs = np.arange(b1, b2 + 1).astype(int)

    print(b1, b2)
    
    if b1 >= -20 and b1 < 20:
        print('skipping this')
        continue;
	
    bigtable = Table.read(bigcat)
 
    sel = np.isin(bigtable['GAIAEDR3_ID'], source_ids)

    tab = astropy.table.vstack((tab, bigtable[sel]))

print(len(tab))

print('writing output!')
tab.write(outfolder + 'photometry.fits', overwrite = True)
