from astropy.table import Table
import astropy
import glob
import os
from tqdm import tqdm
import numpy as np
import glob

datadir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/'


infiles = glob.glob(datadir + 'catalogs/xgall/*.fits')

# with open('/n/home03/vchandra/outerhalo/08_mage/pipeline/control/redux.txt', 'r') as file:
# 	redux = file.read().replace('\n','')

redux = 'v0'

acat = Table()

for file in infiles:
	tab = Table.read(file)
	acat = astropy.table.vstack((acat, tab))
	print(len(acat))


# Reconcile missing spectra

# print('checking how many ACAT stars have spectra...')
# datadir = '/n/holyscratch01/conroy_lab/vchandra/mage/data/reduced/%s/' % redux
# specfiles = glob.glob(datadir + 'spectra/%s/*/*.fits' % redux)
# acatfiles = [datadir + 'spectra/%s/%s/%s' % (redux, row['carton'], row['SPEC_FILE'].strip()) for row in acat]
# has_spec = np.isin(specfiles, acatfiles)
# print('%i stars out of %i in ACAT have downloaded spectra...' % (np.sum(has_spec), len(acat)))

# remove bad columns

acat.remove_columns(['UNWISE_FRACFLUX', 'UNWISE_FLAGS', 'UNWISE_INFO_FLAGS'])

acat = acat[np.argsort(acat['mage_mjd'])]

acat['ACAT_ID'] = np.arange(len(acat))

print('there are %i stars in the ACAT' % len(acat))

acat.write(datadir + 'catalogs/mage_acat.fits', overwrite = True)