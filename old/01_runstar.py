import sys

# sys.path.append('/n/home03/vchandra/software/MINESweeper/')

from minesweeper import fitstar
from astropy.table import Table
import sys,argparse,os
import numpy as np
from numpy.polynomial.polynomial import polyval

from astropy.io import fits

from scipy import constants
speedoflight = constants.c / 1000.0

import socket
hostname = socket.gethostname()
print(hostname)
if hostname[:4] == 'holy': # REWRITE ALL THIS
    specNN = '/n/home03/vchandra/software/MS_files/NN/R12K/modV0_spec_LinNet_R12K_WL445_565.h5' # CHANGE THES
    contNN = '/n/home03/vchandra/software/MS_files/NN/R12K/modV0_cont_LinNet_R12K_WL445_565.h5' #'msdata/lowres/YSTANN_4000_7000_cont.h5' # FIT CONTINUUM_NORMALIZED
    photNN = '/n/home03/vchandra/software/MS_files/VARRV/'
    MISTgrid = '/n/home03/vchandra/software/MS_files/MIST_2.0_spot_EEPtrk_small.h5'
    datadir = '/n/holyscratch01/conroy_lab/vchandra/mage/MS/ms_input/'
    npoints = 250
else:
    specNN = '/Users/vedantchandra/0_research/software/MS_files/NN/R12K/modV0_spec_LinNet_R12K_WL445_565.h5'
    contNN = '/Users/vedantchandra/0_research/software/MS_files/NN/R12K/modV0_cont_LinNet_R12K_WL445_565.h5'
    photNN = '/Users/vedantchandra/0_research/software/MS_files/VARRV/'
    MISTgrid = '/Users/vedantchandra/0_research/software/MS_files/MIST_2.0_spot_EEPtrk_small.h5'
    datadir = '/Users/vedantchandra/0_research/00_outerhalo/05_dr3_giants/mage2022b/DATA/ms_input/'
    npoints = 50

# GET DATA

table = Table.read(datadir + 'photometry_clean.fits')
# table['phot_g_mean_mag_ERR'] = (2.5 / np.log(10)) * (1 / table['phot_g_mean_flux_over_error'])
# table['phot_bp_mean_mag_ERR'] = (2.5 / np.log(10)) * (1 / table['phot_bp_mean_flux_over_error'])
# table['phot_rp_mean_mag_ERR'] = (2.5 / np.log(10)) * (1 / table['phot_rp_mean_flux_over_error'])

#fwhm_sigma = 2.0*np.sqrt(2.0*np.log(2.0))
#wresl_p = np.loadtxt(datadir + 'wresl_model.txt') # polycoef for WRESL sigma from WL Angstrom

sigma_to_fwhm = 2.355

# limit WL range
# w1 = 4470
# w2 = 5635

w1 = 5050
w2 = 5300

def GD(index): # index is row of photometry_clean table

    star = table[index]

    name = str(star['GAIAEDR3_ID'])
    specfile = datadir + name + '.fits'
    print('star is GAIA DR3 %s' % name)
    print('spectrum file is %s' % specfile)

    with fits.open(specfile) as f:
        wl = f[1].data['wave']
        fl = f[1].data['flux']
        ivar = f[1].data['ivar']
        
    # crop spectrum 

    mask = np.ones(len(fl))
    mask[wl < 4470] = 0
    mask[wl > 5635] = 0
    mask_b = mask > 0
    wl = wl[mask_b]
    fl = fl[mask_b]
    ivar = ivar[mask_b]

    spec = {};
    spec['WAVE'] = wl
    spec['FLUX'] = fl
    spec['IVAR'] = ivar
    spec['WRESL'] = np.repeat(7000, len(fl)) # REPLACE WITH PROPER SKY LINE FIT

    print(spec['WRESL'])

    data = {};
    data['mcat'] = star
    data['spec'] = spec

    return data


# SEDpars = Table.read(BOSSdatadir+'SEDpars.dat',format='ascii')

def run(ind,version='V1',overwrite=True):

    print('... Running {}'.format(ind))

    print('    ... Pulling Data')
    data = GD(ind)

    specdata = data['spec']
    mcat = data['mcat']

    # print('    ... Spec Index: {}'.format(mcat['BOSS_INDEX']))

    # print('    ... SDSS INFO:')
    # print('        CLASS: {}'.format(mcat['CLASS']))
    # print('        SUBCLASS: {}'.format(mcat['SUBCLASS']))
    # z = mcat['Z']
    # zerr = mcat['Z_ERR']
    # rv = z * speedoflight
    # rv_err = zerr * speedoflight
    # print('        RV: {0:n} +/- {1:n}'.format(rv,rv_err))

    # print('    ... ELODIE INFO:')
    # print('        Teff: {0:n}'.format(mcat['ELODIE_TEFF']))
    # print('        logg: {0:n}'.format(mcat['ELODIE_LOGG']))
    # print('        Fe/H: {0:n}'.format(mcat['ELODIE_FEH']))

    # SEDpars_i = SEDpars[SEDpars['index'] == ind]

    # print('    ... SED INFO:')
    # print('        Teff: {0:n}'.format(SEDpars_i['Teff'][0]))
    # print('        logg: {0:n}'.format(SEDpars_i['log(g)'][0]))
    # print('        Fe/H: {0:n}'.format(SEDpars_i['[Fe/H]'][0]))
    # print('        a/Fe: {0:n}'.format(SEDpars_i['[a/Fe]'][0]))

    parallax = float(mcat['GAIAEDR3_PARALLAX']) ## APPLY ZP CORRECTION??
    parallax_error = float(mcat['GAIAEDR3_PARALLAX_ERROR'])

    Ebv = float(mcat['EBV'])
    Av = 3.1 * Ebv * 0.86

    print('    ... GaiaEDR3 Parallax = {0:n} +/- {1:n}'.format(parallax,parallax_error))
    print('    ... SFD Av = {0:n}'.format(Av))

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
        # 'SDSS_U':'SDSS_u',
        # 'SDSS_G':'SDSS_g',
        # 'SDSS_R':'SDSS_r',
        # 'SDSS_I':'SDSS_i',
        # 'SDSS_Z':'SDSS_z',
        })

    print('    ... Building Phot')
    phot = {}
    usedphot = {}
    for pp in photbands.keys():
        if ((mcat[pp] > 5.0) &
            (mcat[pp] < 90.0) & 
            (np.abs(mcat[pp+'_ERR']) < 90.0) & 
            ~np.isnan(mcat[pp]) & 
            ~np.isnan(mcat[pp+'_ERR'])
            ):
            filtersys = photbands[pp].split('_')[0]
            filtername = photbands[pp].split('_')[1]

            if filtersys == 'PS':
                if mcat[pp] <= 14.0:
                    continue
                photdata = mcat[pp]
                photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == '2MASS':
                if mcat[pp] <= 5.0:
                    continue
                photdata = mcat[pp]
                photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
            elif (filtersys == 'WISE') or (filtersys == 'UNWISE'):
                if mcat[pp] <= 8.0:
                    continue
                photdata = mcat[pp]
                photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
            elif filtersys == 'SDSS':
                if mcat[pp] <= 12.0:
                    continue
                if filtername == 'u':
                    photdata = mcat[pp] - 0.04
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
                elif filtername == 'i':
                    photdata = mcat[pp] + 0.02
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.02**2.0)))
                else:
                    photdata = mcat[pp]
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == 'GAIAEDR3':
                if filtername == 'G':
                    photdata = mcat[pp]
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'BP':
                    photdata = mcat[pp]
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'RP':
                    photdata = mcat[pp]
                    photerr = np.sqrt((mcat[pp+'_ERR']**2.0)+((0.05**2.0)))
            else:
                photdata = mcat[pp]
                photerr = np.sqrt(mcat[pp+'_ERR']**2 + 0.02**2)

            usedphot[photbands[pp]] = pp
            phot[photbands[pp]] = [float(photdata),float(photerr)]

    print('    ... Building Spec')
    spec = {}
    cond = (specdata['IVAR'] != 0.0) & (specdata['WAVE'] > w1) & (specdata['WAVE'] < w2) # crop to wavelength range
    spec['WAVE']   = specdata['WAVE'][cond]
    spec['FLUX']   = specdata['FLUX'][cond]
    spec['E_FLUX'] = 1.0/np.sqrt(specdata['IVAR'][cond])
    spec['WRESL']  = specdata['WRESL'][cond]

    cond = np.isfinite(spec['FLUX']) & np.isfinite(spec['E_FLUX'])
    spec['WAVE']   = spec['WAVE'][cond]
    spec['FLUX']   = spec['FLUX'][cond]
    spec['E_FLUX'] = spec['E_FLUX'][cond]
    spec['WRESL']  = spec['WRESL'][cond]

    # cond = (spec['WAVE'] > 3850.0) & (spec['WAVE'] < 8900.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    # cond = (spec['WAVE'] > 4000.0) & (spec['WAVE'] < 7000.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    # mask out Na doublet due to ISM absorption
    cond = (spec['WAVE'] < 5850.0) | (spec['WAVE'] > 5950.0)
    spec['WAVE']   = spec['WAVE'][cond]
    spec['FLUX']   = spec['FLUX'][cond]
    spec['E_FLUX'] = spec['E_FLUX'][cond]
    spec['WRESL']  = spec['WRESL'][cond]

    # # mask out telluric features
    # cond = (spec['WAVE'] < 7500.0) | (spec['WAVE'] > 7700.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    medflux = np.median(spec['FLUX'])
    spec['FLUX']   = spec['FLUX']/medflux
    spec['E_FLUX'] = spec['E_FLUX']/medflux
    print('spectrum from %i to %i wavelengths' % (np.min(spec['WAVE']),np.max(spec['WAVE'])))
    print('    ... Building Input Dict')
    # build input dict
    inputdict = {}

    inputdict['specANNpath'] = specNN
    inputdict['contANNpath'] = contNN
    inputdict['NNtype'] = 'LinNet'
    inputdict['photANNpath'] = photNN
    inputdict['MISTpath'] = MISTgrid
    inputdict['isochrone_prior'] = True
    inputdict['ageweight'] = True
    inputdict['mistdif'] = True
    inputdict['Rvfree'] = False

    inputdict['phot'] = phot

    inputdict['spec'] = {}
    inputdict['spec']['obs_wave']  = spec['WAVE']
    inputdict['spec']['obs_flux']  = spec['FLUX']
    inputdict['spec']['obs_eflux'] = spec['E_FLUX']
    inputdict['spec']['modpoly'] = True
    inputdict['spec']['convertair'] = False

    maxdist = 1000.0/(parallax-10.0*parallax_error)
    if maxdist < 0.0:
        maxdist = np.inf
    mindist = 1000.0/(parallax+10.0*parallax_error)

    print('    ... Building Sampler Options')

    # set parameter for sampler
    inputdict['sampler'] = {}
    inputdict['sampler']['samplertype'] = 'Static'
    inputdict['sampler']['samplemethod'] = 'rwalk'
    # inputdict['sampler']['samplemethod'] = 'unif'
    inputdict['sampler']['npoints'] = npoints # 500
    # inputdict['sampler']['npoints'] = 150
    inputdict['sampler']['samplerbounds'] = 'multi'
    inputdict['sampler']['flushnum'] = 250 # 250
    inputdict['sampler']['delta_logz_final'] = 0.01 # 0.01
    inputdict['sampler']['bootstrap'] = 0
    inputdict['sampler']['walks'] = 25
    # inputdict['sampler']['walks'] = 3

    ############## Priors ################
    print('    ... Building Prior Options')
    inputdict['priordict'] = {}

    inputdict['priordict']['EEP'] = {'pv_uniform':[1,808]}
    inputdict['priordict']['initial_Mass'] = {'pv_uniform':[0.5,2.0]}
    inputdict['priordict']['initial_[Fe/H]'] = {'pv_uniform':[-4.0,0.25]}
    inputdict['priordict']['initial_[a/Fe]'] = {'pv_uniform':[-0.2,0.6]}
    #inputdict['priordict']['initial_[a/Fe]'] = {'fixed' : 0}
    inputdict['priordict']['Dist']   = {'pv_uniform':[max([1.0,mindist]),min([maxdist,300000.0])]}
    inputdict['priordict']['Av']     = {'pv_tgaussian':[0.0,5.0*Av,Av,Av*0.1]}

    inputdict['priordict']['Age'] = {'uniform':[6.0,14.0]}
    inputdict['priordict']['Parallax'] =  {'gaussian':[parallax,parallax_error]}

    numpoly = 4
    coeffarr = [[1.5,1.0],[0.0,0.5],[0.0,0.5]] +[[0.0,0.25] for ii in range(numpoly-3)]
    inputdict['priordict']['blaze_coeff'] = coeffarr
    inputdict['priordict']['Vrad'] = {'pv_uniform':[-350.0, 350.0]}
    #inputdict['priordict']['Vrad'] = {'fixed': mcat['mle_rv']} # MAYBE REPLACE WITH GAUSSIAN 
    # inputdict['priordict']['Vrot'] = {'pv_tgaussian':[0.0, 200.0, 25.0, 25.0]}
    # inputdict['priordict']['Vrot'] = {'pv_uniform':[0.0, 100.0]}
    inputdict['priordict']['Vrot'] = {'pv_beta':[1.05,1.5,0.0,250.0]}
    #inputdict['priordict']['Vrot'] = {'fixed': 5}
    #inputdict['priordict']['Inst_R'] = {'fixed' : spec['WRESL']}
    inputdict['priordict']['Inst_R'] = {'fixed':8000.0}

    inputdict['priordict']['IMF'] = {'IMF_type':'Kroupa'}
    # inputdict['priordict']['GAL'] = {'lb_coords':[float(mcat['L']),float(mcat['B'])]}

    # inputdict['priordict']['GALAGE'] = {}
    # inputdict['priordict']['GALAGE']['lb_coords'] = [float(mcat['L']),float(mcat['B'])]
    # inputdict['priordict']['GALAGE']['pars'] = ({
    #     'thin': {'min':1.0,'max':14.0},
    #     'thick':{'min':6.0,'max':14.0,'mean':10.0,'sigma':2.0},
    #     'halo': {'min':8.0,'max':14.0,'mean':12.0,'sigma':2.0},
    #     })

    # inputdict['priordict']['VROT'] = {'giant':{'a':-10.0,'c':7.0,'n':1.0},'dwarf':{'a':-10.0,'c':10.0,'n':0.4}}

    inputdict['output'] = datadir + 'output/samp_run_{0}_spec_{1}_{2}.dat'.format(ind,str(mcat['GAIAEDR3_ID']),version)
    print('    ... Writing to: {}'.format(inputdict['output']))

    print('---------------')
    if 'spec' in inputdict.keys():
        print('    Median Spec Flux: ')
        print('       {0}'.format(np.median(inputdict['spec']['obs_flux']*medflux)))
        print('    Median Spec Err_Flux:')
        print('       {0}'.format(np.median(inputdict['spec']['obs_eflux']*medflux)))
        print('    Fitting Wavelengths:')
        print('       {0} -- {1}'.format(
            min(inputdict['spec']['obs_wave']),
            max(inputdict['spec']['obs_wave'])))
        print('    Min Resolution of Spectrum: {0:n}'.format(np.min(spec['WAVE']/(spec['WRESL'] * sigma_to_fwhm))))
        print('    Max Resolution of Spectrum: {0:n}'.format(np.max(spec['WAVE']/(spec['WRESL'] * sigma_to_fwhm))))        
        print('    Mean Resolution of Spectrum: {0:n}'.format(np.median(spec['WAVE']/(spec['WRESL'] * sigma_to_fwhm))))
        print('    Fitting w/ ANN: {0}'.format(inputdict['NNtype']))

    print('---------------')
    if 'phot' in inputdict.keys():
        print('    PHOT:')
        for kk in inputdict['phot'].keys():
            print('       {0} = {1} +/- {2}'.format(kk,inputdict['phot'][kk][0],inputdict['phot'][kk][1]))

    print('    PRIORS:')
    for kk in inputdict['priordict'].keys():
      if kk == 'blaze_coeff':
           pass
      elif kk in ['IMF','GAL','GALAGE','VROT','VTOT','ALPHA']:
           print('       Turned on {0} prior'.format(kk))
      else:
           try:
                for kk2 in inputdict['priordict'][kk].keys():
                     if (kk2 == 'uniform') or (kk2 == 'pv_uniform'):
                          print('       {0}: min={1} max={2}'.format(kk,inputdict['priordict'][kk][kk2][0],inputdict['priordict'][kk][kk2][1]))
                     if (kk2 == 'gaussian') or (kk2 == 'pv_gaussian'):
                          print('       {0}: N({1},{2})'.format(kk,inputdict['priordict'][kk][kk2][0],inputdict['priordict'][kk][kk2][1]))
                     if (kk2 == 'tgaussian') or (kk2 == 'pv_tgaussian'):
                          print('       {0}: N({1},{2}) [{3} - {4}]'.format(
                               kk,inputdict['priordict'][kk][kk2][2],inputdict['priordict'][kk][kk2][3],
                               inputdict['priordict'][kk][kk2][0],inputdict['priordict'][kk][kk2][1]))
                     if kk2 == 'pv_exp':
                          print('       {0}: EXP({1},{2})'.format(kk,inputdict['priordict'][kk][kk2][0],inputdict['priordict'][kk][kk2][1]))
                     if kk2 == 'pv_loguniform':
                          print('       {0}: log(uniform({1},{2}))'.format(kk,inputdict['priordict'][kk][kk2][0],inputdict['priordict'][kk][kk2][1]))
                     if kk2 == 'fixed':
                          if kk == 'Inst_R':
                            if not isinstance(inputdict['priordict'][kk][kk2],float):
                                print('       {0}: fixed'.format(kk))
                            else:
                                print('       {0}: fixed({1})'.format(kk,inputdict['priordict'][kk][kk2]))
                          else:
                            print('       {0}: fixed({1})'.format(kk,inputdict['priordict'][kk][kk2]))

           except:
                print('       {0}: {1}'.format(kk,inputdict['priordict'][kk]))

    print('--------------')


    FS = fitstar.FitMS()
    results = FS.run(inputdict=inputdict)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("starind", help="index of star being processed",type=int)
    parser.add_argument("--version", "-v", help="run version",type=str,default='VX')
    parser.add_argument("--overwrite", "-o", help="overwrite previous run", action='store_true')
    args = parser.parse_args()

    run(args.starind,version=args.version,overwrite=args.overwrite)
