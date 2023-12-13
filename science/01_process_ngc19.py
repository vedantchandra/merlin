#%%

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.time import Time
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from chandra import outerhalo as oh
import os
import glob
import pandas as pd
from scipy import stats

plt.style.use('~/vedant.mplstyle')


# %%

datadir = '/n/holystore01/LABS/conroy_lab/Lab/h3/mocks/lmc/'
files = np.sort(glob.glob(datadir + 'MWLMC5*icrs_nocuts.fits'))
print(files)

outdir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/ngc19/'

# %%
file = files[1]
# %%
tab = Table.read(file)


# %%

sel = (

    (tab['distance'] > 20)

)
# %%
np.sum(sel) / len(sel)
# %%

selcol = ['X_gal', 'Y_gal', 'Z_gal', 
          'Vx_gal', 'Vy_gal', 'Vz_gal',
            'l', 'b', 'distance',
            'V_gsr']

# %%

Table(np.array(tab[selcol])).write(outdir + file.split('/')[-1].split('.')[0] + '_trimmed.h5')

# %%
