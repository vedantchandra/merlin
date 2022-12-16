import numpy as np
from astropy.table import Table,hstack,vstack,join
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob

from scipy import constants
speedoflight = constants.c / 1000.0


datadir = '/n/holyscratch01/conroy_lab/vchandra/mage/'

with open('/n/home03/vchandra/outerhalo/08_mage/pipeline/control/redux.txt', 'r') as file:
    redux = file.read().replace('\n','')

def getdata(GaiaID = None, acat_id = None, date = None,
                mask_hbeta = False):    
    
    acat = Table.read(datadir + 'catalogs/mage_acat.fits')
    acat = acat.filled(99.0)
    
    if acat_id is not None:
        row = acat[acat['ACAT_ID'] == acat_id]
    elif GaiaID is not None and date is not None:
        row = acat[(acat['GAIAEDR3_ID'] == GaiaID)&(acat['date'] == date)]
    # elif index is not None:
    #     row = acat[acat['ACAT_ID'] == index][0]
    else:
        print('must pass either GaiaID+date or ACAT index!!')
        raise

    if len(row) == 1:
        row = row[0]
    elif len(row) == 0:
        print('no matches in table! cannot get data...')
        raise
    elif len(row) > 1:
        print('getdata query does not return unique row in ACAT')
        raise
    
    filepath = row['specfile']

    with fits.open(filepath) as f:
    
        header = f[0].header

        phot = dict(row)

        # pull individual spectrum parts
        wave = 10**f[1].data['LOGLAM']
        flux = f[1].data['FLUX']
        ivar = f[1].data['IVAR']
        andmask = f[1].data['AND_MASK']
        ormask  = f[1].data['OR_MASK']
        #lsf = f[1].data['WRESL'] * lsf_fudge
        wresl = (np.log(10) * header['CD1_1'] * f[1].data['WDISP']) * wave

    ## QUALITY CUTS TO SPECTRUM

    if mask_hbeta:
        hbsel = (wave > 4830.0) & (wave < 4890)
        ivar[hbsel] = 0.0


    cond = (
        np.isfinite(flux) & 
        (ivar > 0.0) & 
        (wresl> 0.0) &
        (wave > 4750) &
        (wave < 5550.0) 
    )

    wave   = wave[cond]
    flux   = flux[cond]
    ivar = ivar[cond]
    wresl = wresl[cond]

    medflux = np.median(flux)
    flux /= medflux
    ivar *= medflux**2

    return {'phot':phot,'spec':[wave,flux,ivar,andmask,ormask,wresl]}
