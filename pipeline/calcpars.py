from astropy.table import Table
import numpy as np
import sys,os,argparse,glob
from quantiles import quantile

from scipy import constants
speedoflight = constants.c / 1000.0

datadir = '/n/holyscratch01/conroy_lab/vchandra/sdss5/'
specNN = datadir + 'ms/NN/modV0_spec_LinNet_R5K_WL445_565.h5'
contNN = datadir + 'ms/NN/modV0_cont_LinNet_R12K_WL445_565.h5'
photNN = datadir + 'ms/VARRV/'
SBlib = datadir + 'ms/CKC/ckc_R500.h5'
MISTgrid = datadir + 'ms/MIST_2.0_spot_EEPtrk_small.h5'
outdir = datadir
NNtype = 'LinNet'

from minesweeper import genmod
from minesweeper.fastMISTmod import GenMIST 

from getdata import getdata
from phaseafy import phaseafy

fwhm_sigma = 2.0*np.sqrt(2.0*np.log(2.0))

photbands = ({
    'GAIAEDR3_G':'GaiaEDR3_G',
    'GAIAEDR3_BP':'GaiaEDR3_BP',
    'GAIAEDR3_RP':'GaiaEDR3_RP',
    'PS_G':'PS_g',
    'PS_R':'PS_r',
    'PS_I':'PS_i',
    'PS_Z':'PS_z',
    'PS_Y':'PS_y',
    'TMASS_J':'2MASS_J',
    'TMASS_H':'2MASS_H',
    'TMASS_K':'2MASS_Ks',
    'UNWISE_W1':'WISE_W1',
    'UNWISE_W2':'WISE_W2',
    'SDSS_U':'SDSS_u',
    'SDSS_G':'SDSS_g',
    'SDSS_R':'SDSS_r',
    'SDSS_I':'SDSS_i',
    'SDSS_Z':'SDSS_z',
    })

usedbands_i = list(photbands.keys())
usedbands = [photbands[x] for x in usedbands_i]

GM = genmod.GenMod() 
GM._initspecnn(
    nnpath=specNN,
    NNtype=NNtype,
    Cnnpath=contNN)
GM._initphotnn(
    usedbands,
    nnpath=photNN,
    )

GMIST = GenMIST(
    MISTpath=MISTgrid,
    ageweight=False)

def bfspec(specdata,bf):
    spec = {}
    cond = specdata[2] != 0.0
    spec['WAVE']   = specdata[0][cond]
    spec['FLUX']   = specdata[1][cond]
    spec['E_FLUX'] = 1.0/np.sqrt(specdata[2][cond])
    spec['LSF']  = specdata[-1][cond]

    cond = np.isfinite(spec['FLUX']) & np.isfinite(spec['E_FLUX']) & (spec['LSF'] > 0.0)
    spec['WAVE']   = spec['WAVE'][cond]
    spec['FLUX']   = spec['FLUX'][cond]
    spec['E_FLUX'] = spec['E_FLUX'][cond]
    spec['LSF']    = spec['LSF'][cond]

    # create the WRESL array
    spec['WRESL'] = (spec['WAVE'] * spec['LSF']) / speedoflight

    # cond = (spec['WAVE'] > 3850.0) & (spec['WAVE'] < 8900.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    # cond = (spec['WAVE'] > 4000.0) & (spec['WAVE'] < 7000.0)
    # cond = (spec['WAVE'] > 4455.0) & (spec['WAVE'] < 5645.0)
    # cond = (spec['WAVE'] > 5000.0) & (spec['WAVE'] < 5500.0)
    cond = (spec['WAVE'] > 4750.0) & (spec['WAVE'] < 5500.0)
    spec['WAVE']   = spec['WAVE'][cond]
    spec['FLUX']   = spec['FLUX'][cond]
    spec['E_FLUX'] = spec['E_FLUX'][cond]
    spec['WRESL']  = spec['WRESL'][cond]

    # # mask out Na doublet due to ISM absorption
    # cond = (spec['WAVE'] < 5850.0) | (spec['WAVE'] > 5950.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    # # mask out telluric features
    # cond = (spec['WAVE'] < 7500.0) | (spec['WAVE'] > 7700.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    medflux = np.median(spec['FLUX'])
    # spec['FLUX']   = spec['FLUX']/medflux
    # spec['E_FLUX'] = spec['E_FLUX']/medflux

    ii = 0
    for x in bf.keys():
        if ('pc_' in x) and ('err' not in x):
            ii += 1

    bfpars = [bf['Teff'],bf['log(g)'],bf['[Fe/H]'],bf['[a/Fe]'],bf['Vrad'],bf['Vrot'],np.nan,spec['WRESL']]
    bfpcpars = [bf['pc_{}'.format(x)] for x in range(ii)]
    bfspecpars = bfpars + bfpcpars

    bfmod = GM.genspec(bfspecpars,outwave=spec['WAVE'],modpoly=True)
    bfflux = bfmod[1]*medflux

    chisq_spec = np.nansum([((m-d)/s)**2.0 for m,d,s in zip(bfflux,spec['FLUX'],spec['E_FLUX'])])
    nspecpix = len(bfflux)

    return chisq_spec,nspecpix

def bfphot(phdata,bf):

    phot = {}
    usedphot = {}
    for pp in usedbands_i:
        if ((phdata[pp] > 5.0) &
            (phdata[pp] < 90.0) & 
            (np.abs(phdata[pp+'_ERR']) < 90.0) & 
            ~np.isnan(phdata[pp]) & 
            ~np.isnan(phdata[pp+'_ERR'])
            ):
            filtersys = photbands[pp].split('_')[0]
            filtername = photbands[pp].split('_')[1]

            if filtersys == 'PS':
                if phdata[pp] <= 14.0:
                    continue
                photdata = phdata[pp]
                photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == '2MASS':
                if phdata[pp] <= 5.0:
                    continue
                photdata = phdata[pp]
                photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
            elif (filtersys == 'WISE') or (filtersys == 'UNWISE'):
                if phdata[pp] <= 8.0:
                    continue
                photdata = phdata[pp]
                photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
            elif filtersys == 'SDSS':
                if phdata[pp] <= 12.0:
                    continue
                if filtername == 'u':
                    photdata = phdata[pp] - 0.04
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
                elif filtername == 'i':
                    photdata = phdata[pp] + 0.02
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.02**2.0)))
                else:
                    photdata = phdata[pp]
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == 'GaiaEDR3':
                if filtername == 'G':
                    photdata = phdata[pp+'_CORRECTED']
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'BP':
                    photdata = phdata[pp]
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'RP':
                    photdata = phdata[pp]
                    photerr = np.sqrt((phdata[pp+'_ERR']**2.0)+((0.05**2.0)))
            else:
                photdata = phdata[pp]
                photerr = phdata[pp+'_ERR']

            usedphot[photbands[pp]] = pp
            phot[photbands[pp]] = [float(photdata),float(photerr)]

    filterbands = list(phot.keys())

    sedpars = ([
        bf['Teff'],
        bf['log(g)'],
        bf['[Fe/H]'],
        bf['[a/Fe]'],
        bf['log(R)'],
        bf['Dist']*1000.0,
        bf['Av'],
        ])
    sed = GM.genphot(sedpars)
    modmag = [sed[kk] for kk in filterbands]

    chisq_phot = np.nansum([
        ((m-d)/s)**2.0 for m,d,s in zip(
            modmag,
            [phot[x][0] for x in filterbands],
            [phot[x][1] for x in filterbands])
        ])
    nbands = len(filterbands)

    return chisq_phot,nbands


def run(index=None,GaiaID=None,version='VX',catalog = None, acat_id = None):
    # read in SEGUE data
    if index is not None:
        data = getdata(index=index)
    if GaiaID is not None:
        data = getdata(GaiaID=GaiaID)
    if acat_id is not None:
        data = getdata(acat_id = acat_id)
    try:
        assert data['phot']
    except AssertionError:
        print('Warning: User did not pass a valid selection criteria')
        raise

    try:
        assert data['phot']
    except AssertionError:
        print('Warning: User did not pass a valid selection criteria')
        raise

    if catalog is None:
        catalog = 'solo'

    photcat = data['phot']

    samplefile = '{OUTDIR}{CATALOG}/{VER}/mwm_gaiaID_{GAIAID}_fieldID_{FIELDID}_mjd_{MJD}_catID_{CATID}_{VER}_samp.dat'.format(
            OUTDIR = outdir + 'samples/',
            FIELDID=data['phot']['FIELD'],
            GAIAID=data['phot']['GAIAEDR3_ID'],
            CATID=data['phot']['CATALOGID'],
            MJD=data['phot']['MJD'],
            VER=version,
            CATALOG=catalog)
    samplefile_gz = samplefile + '.gz'

    try:
        try:
            samp = Table.read(samplefile,format='ascii')
        except:
            samp = Table.read(samplefile_gz,format='ascii')
        if samp['delta(log(z))'][-1] > 0.01:
            print('CALCPARS: delta(log(z)) DID NOT CONVERGE!')
            return
    except:
        print('CALCPARS: Problem with reading sample file!')
        return

    samp['Pr'] = np.exp(samp['log(wt)']-samp['log(z)'][-1])
    samp = samp[samp['Pr'] > 1E-10]

    samp['Dist'] = samp['Dist']/1000.0
    if 'Para' not in samp.keys():
        samp['Para'] = 1.0/samp['Dist']
    if 'Age' not in samp.keys():
        samp['Age'] = 10.0**(samp['log(Age)']-9.0)

    fitpars = ([
        x for x in samp.keys() if x not in 
        ['Iter','Agewgt','log(lk)','log(vol)','log(wt)','h','nc','log(z)','delta(log(z))','Pr']
        ])

    pardict = {}
    for ff in fitpars:
        pardict[ff] = quantile(samp[ff],[0.16,0.5,0.84],weights=samp['Pr'])

    # define an output dictionary
    outdict = {}

    # add photcat info
    for nn in photcat.keys():
        outdict[nn] = photcat[nn]

    # delete OBJID as it is a list
    # del outdict['OBJID']

    # convert any byte string to ascii
    for kk in outdict.keys():
        if isinstance(outdict[kk],bytes):
            outdict[kk] = outdict[kk].decode('ascii')

    # add best fit pars w/ errors
    for ff in pardict.keys():
        ll,med,ul = pardict[ff]
        outdict[ff] = med
        outdict[ff+'_lerr'] = med - ll
        outdict[ff+'_uerr'] = ul - med
        outdict[ff+'_err']  = (ul - ll)/2.0

    # add log(z) log(wt) log(lk)
    outdict['lnZ'] = samp['log(z)'][-1]
    outdict['lnL'] = samp['log(lk)'].max()
    outdict['lnP'] = np.log(samp['Pr'].max())

    # calc chi-sq and add
    outdict['chisq_spec'],outdict['nspecpix'] = bfspec(data['spec'],outdict)
    outdict['chisq_phot'],outdict['nbands'] = bfphot(data['phot'],outdict)

    # doing phase-a-fy
    PH = phaseafy()

    # add phaseafy stuff
    phpar = ([
        'R_gal','X_gal','Y_gal','Z_gal',
        'Vx_gal','Vy_gal','Vz_gal',
        'Vr_gal','Vphi_gal','Vtheta_gal','V_tan','V_gsr',
        'Lx','Ly','Lz','Ltot'])
    for POTN in PH.potentials.keys():
        phpar.append('E_kin_{0}'.format(POTN))
        phpar.append('E_pot_{0}'.format(POTN))
        phpar.append('E_tot_{0}'.format(POTN))
        phpar.append('circLz_{0}'.format(POTN))
        phpar.append('circLtot_{0}'.format(POTN))

    for pp in phpar:
        outdict[pp] = np.nan 
        outdict[pp+'_err'] = np.nan 
    outdict = PH.calcphase(outdict,nsamples=50000,verbose=False)

    outfile = '{OUTDIR}{CATALOG}/{VER}/mwm_gaiaID_{GAIAID}_fieldID_{FIELDID}_mjd_{MJD}_catID_{CATID}_{VER}.pars'.format(
        OUTDIR=outdir + 'pars/',
        FIELDID=data['phot']['FIELD'],
        GAIAID=data['phot']['GAIAEDR3_ID'],
        CATID=data['phot']['CATALOGID'],
        MJD=data['phot']['MJD'],
        VER=version,
        CATALOG = catalog)

    print('... writing file to: {0}'.format(outfile))

    with open(outfile,'w') as ofile:
        parkeys = list(outdict.keys())
        for kk in parkeys:
            ofile.write('{} '.format(kk))
        ofile.write('\n')
        for kk in parkeys:
            ofile.write('{} '.format(outdict[kk]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--survey',help='SEGUE Survey ID',type=str,
    #     choices=['SEGUE','SEGUE_clusters'],default='SEGUE')
    # parser.add_argument('--tileID',help='tile ID number',type=int,default=1660)
    parser.add_argument('--index',help='Index of star in acat',type=int,default=None)
    parser.add_argument('--GaiaID',help='Gaia EDR3 ID of star',type=int,default=None)
    # parser.add_argument('--FiberID',help='Fiber ID of star',type=int,default=None)
    # parser.add_argument('--mjd',help='MJD of plate to run',type=int,default=None)
    parser.add_argument("--version", "-v", help="run version",type=str,default='VX')
    args = parser.parse_args()
    run(
        index=args.index,
        GaiaID=args.GaiaID,
        version=args.version)