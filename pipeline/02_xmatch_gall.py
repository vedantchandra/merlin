# xmatches supplied clean table to all-band photometry on Cannon

import numpy as np
from astropy.io import ascii
from astropy.table import Table
import glob
import astropy
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
import sys

print('starting script')

#infile = '/n/holyscratch01/conroy_lab/vchandra/mage/MS/ms_input/source_ids.txt'

infile = '/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/spall.fits'
outfolder = '/n/holyscratch01/conroy_lab/vchandra/mage/catalogs/xgall/'

try:
	os.mkdir(outfolder)
except:
	print('folder exists!')

# XMATCH Charlie's Gaia EDR3 photometric tables to my magE stars
# Loop over Charlie's files, match my file, save as h5

charlie_cats = '/n/holystore01/LABS/conroy_lab/Lab/gaia/edr3/gall2/catalogs/'
charlie_cats = glob.glob(charlie_cats + '*')

galidx = int(sys.argv[1])
bigcat = charlie_cats[galidx]
table = Table.read(infile)

#outtab = Table();

c = SkyCoord(ra=table['tdb_ra'] * u.degree, dec = table['tdb_dec']*u.degree)

# for bigcat in charlie_cats:
print('-----------')
print('Big Catalog: %s' % bigcat.split('/')[-1])

bs = bigcat.split('/')[-1].split('b')[-1].split('_')
b1 = float(bs[0])
b2 = float(bs[1][:5])
bs = np.arange(b1, b2 + 1).astype(int)

print(b1, b2)

bigtable = Table.read(bigcat)

max_sep = np.repeat(3.0, len(table)) * u.arcsec
notgaia = (table['tdb_selection'] == 'rvs') | (table['tdb_selection'] == 'tell')
max_sep[notgaia] = 20 * u.arcsec

catalog = SkyCoord(ra=bigtable['RA'] * u.degree, dec = bigtable['DEC'] * u.degree)
idx, d2d, d3d = c.match_to_catalog_sky(catalog)
sep_constraint = d2d < max_sep
c_matches = table[sep_constraint]
cat_matches = bigtable[idx[sep_constraint]]

outtab = astropy.table.hstack((c_matches, cat_matches))

#outtab = astropy.table.vstack((outtab, tab))

del bigtable

print('output table has %i rows' % len(outtab))
print('writing output!')
outtab.write(outfolder + 'mage_xgall_%i.fits' % galidx, overwrite = True)