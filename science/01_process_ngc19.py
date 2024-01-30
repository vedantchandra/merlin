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
import random

plt.style.use('~/vedant.mplstyle')


# %%

datadir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/ngc19/tracers/'
files = np.sort(glob.glob(datadir + '*.txt'))
print(files)

outdir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/ngc19/tracers/processed/'

# %%
for file in files:
    print('reading %s' % file)

    p = 0.01  # 1% of the lines

    df = pd.read_csv(
            file,
            names = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'm'],
            skiprows=lambda i: i>0 and random.random() > p,
            delimiter=' ',
    )

    print('%i rows' % len(df))

    coo = SkyCoord(x = df['x'] * u.kpc, y = df['y'] * u.kpc, z = df['z'] * u.kpc,
                v_x = df['vx'] * u.km / u.s, v_y = df['vy'] * u.km / u.s, v_z = df['vz'] * u.km / u.s,
                frame = 'galactocentric')


    cooeq = coo.icrs
    coogal = coo.galactic

    df['ra'] = cooeq.ra.value
    df['dec'] = cooeq.dec.value
    df['pmra'] = cooeq.pm_ra_cosdec.value
    df['pmdec'] = cooeq.pm_dec.value
    df['distance'] = cooeq.distance.value
    df['rv'] = cooeq.radial_velocity.value

    df['l'] = coogal.l.value
    df['b'] = coogal.b.value

    # # %%
    # plt.scatter(df['ra'], df['dec'])
    # # %%

    # plt.hist(df['distance'])

    # # %%
    # rgal = np.sqrt(df['x']**2 + df['y']**2 + df['z']**2)
    # # %%
    # plt.hist(rgal)
    # #%%%

    name = file.split('/')[-1].replace('.txt', '')
    outname = outdir + name + '_observed.h5'

    print('writing %s' % outname)

    df.to_hdf(outname, key = 'table', mode = 'w', format = 'fixed')
# %%
