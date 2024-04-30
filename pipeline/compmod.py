from minesweeper.fitutils import polycalc
from minesweeper.fitutils import airtovacuum

from astropy.table import Table
import numpy as np
import sys,os,argparse,glob
from quantiles import quantile

import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

from scipy.ndimage import zoom,gaussian_filter
from scipy.stats import gaussian_kde,scoreatpercentile
from scipy import constants
speedoflight = constants.c / 1000.0

# useful constants
# speedoflight = 2.997924e+10
speedoflight_kms = 2.997924e+5
speedoflight_nms = 2.997924e+17
lsun = 3.846e33
pc = 3.085677581467192e18  # in cm
jansky_cgs = 1e-23
# value to go from L_sun to erg/s/cm^2 at 10pc
log_rsun_cgs = np.log10(6.955) + 10.0
log_lsun_cgs = np.log10(lsun)
log4pi = np.log10(4 * np.pi)

# import socket,os
# hostname = socket.gethostname()
# if hostname[:4] == 'holy':
#     holypath   = os.environ.get('HOLYSCRATCH')
#     # specNN = '{}/conroy_lab/pacargile/ThePayne/Hecto_FAL/lowres/YSTANN_4000_7000_spec.h5'.format(holypath)
#     # contNN = '{}/conroy_lab/pacargile/ThePayne/Hecto_FAL/lowres/YSTANN_4000_7000_cont.h5'.format(holypath)
#     specNN = '{}/conroy_lab/pacargile/ThePayne/train/optfal/v256/modV0_spec_LinNet_R5K_WL445_565.h5'.format(holypath)
#     contNN = '{}/conroy_lab/pacargile/ThePayne/train/optfal/v256/modV0_cont_LinNet_R12K_WL445_565.h5'.format(holypath)
#     NNtype = 'LinNet'
#     photNN = '{}/conroy_lab/pacargile/ThePayne/SED/VARRV/'.format(holypath)
#     SBlib = '{}/conroy_lab/pacargile/CKC/ckc_R500.h5'.format(holypath)
#     MISTgrid = '{}/conroy_lab/pacargile/MIST/MIST_2.0_spot_EEPtrk_small.h5'.format(holypath)
#     datadir = '/n/holystore01/LABS/conroy_lab/Lab/SEGUE/data/'
#     outdir  = '/n/holyscratch01/conroy_lab/pacargile/SEGUE/'
# else:
#     # specNN = '/Users/pcargile/Astro/ThePayne/YSdata/lowres/YSTANN_4000_7000_spec.h5'
#     # specNN = '/Users/pcargile/Astro/ThePayne/YSdata/YSTANN.h5'
#     # contNN = '/Users/pcargile/Astro/ThePayne/YSdata/lowres/YSTANN_4000_7000_cont.h5'
#     specNN = '/Users/pcargile/Astro/ThePayne/train_grid/optfal/v256/modV0_spec_LinNet_R5K_WL445_565.h5'
#     contNN = '/Users/pcargile/Astro/ThePayne/train_grid/optfal/v256/modV0_cont_LinNet_R12K_WL445_565.h5'
#     NNtype = 'LinNet'
#     photNN = '/Users/pcargile/Astro/GITREPOS/ThePayne/data/photANN/'
#     SBlib = '/Users/pcargile/Astro/ckc/ckc_R500.h5'
#     MISTgrid = '/Users/pcargile/Astro/MIST/MIST_v2.0_spot/MIST_2.0_spot_EEPtrk_small.h5'
#     datadir = '/Users/pcargile/Astro/SEGUE/data/'
#     outdir  = '/Users/pcargile/Astro/SEGUE/'

# datadir = '/n/holyscratch01/conroy_lab/vchandra/sdss5/'
# specNN = datadir + 'ms/NN/modV0_spec_LinNet_R5K_WL445_565.h5'
# contNN = datadir + 'ms/NN/modV0_cont_LinNet_R12K_WL445_565.h5'
# photNN = datadir + 'ms/VARRV/'
# MISTgrid = datadir + 'ms/MIST_2.0_spot_EEPtrk_small.h5'
# outdir = datadir
# NNtype = 'LinNet'

specNN = '/n/home03/vchandra/software/MS_files/NN/R12K/modV0_spec_LinNet_R12K_WL445_565.h5' # CHANGE THES
contNN = '/n/home03/vchandra/software/MS_files/NN/R12K/modV0_cont_LinNet_R12K_WL445_565.h5' #'msdata/lowres/YSTANN_4000_7000_cont.h5' # FIT CONTINUUM_NORMALIZED
photNN = '/n/home03/vchandra/software/MS_files/VARRV/'
MISTgrid = '/n/home03/vchandra/software/MS_files/MIST_2.0_spot_EEPtrk_small.h5'
datadir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/'
outdir = datadir
NNtype = 'LinNet'
SBlib = '/n/home03/vchandra/software/MS_files/CKC/ckc_R500.h5'

from minesweeper import genmod
from minesweeper.fastMISTmod import GenMIST 

from getdata import getdata

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

GM = genmod.GenMod() 
GM._initspecnn(
    nnpath=specNN,
    NNtype=NNtype,
    Cnnpath=contNN)
GM._initphotnn(
    [photbands[x] for x in photbands.keys()],
    nnpath=photNN,
    )

GMIST = GenMIST(
    MISTpath=MISTgrid,
    ageweight=False)


def mksed(axSED,samples,photdata,bfdict):
    sedstr = (
        'GaiaEDR3 G = {0:.2f}'.format(photdata['GaiaEDR3_G'][0])
        )
    if 'PS_g' in photdata.keys():
        sedstr += '\n PS g = {0:.2f}'.format(photdata['PS_g'][0])
    if '2MASS_J' in photdata.keys():
        sedstr += '\n 2MASS J = {0:.2f}'.format(photdata['2MASS_J'][0])
    if 'WISE_W1' in photdata.keys():
        sedstr += '\n WISE W1 = {0:.2f}'.format(photdata['WISE_W1'][0])


    axSED.text(0.97,0.97,sedstr,
        horizontalalignment='right',verticalalignment='top', 
        transform=axSED.transAxes,fontsize=8)


    import star_basis
    import photsys
    from ccm_curve import ccm_curve

    SB = star_basis.StarBasis(
        libname=SBlib,
        use_params=['logt','logg','feh'],
        n_neighbors=1)

    WAVE_d = photsys.photsys()
    photbands_i = WAVE_d.keys()
    photbands = [x for x in photbands_i if x in photdata.keys()]
    WAVE = {pb:WAVE_d[pb][0] for pb in photbands}
    zeropts = {pb:WAVE_d[pb][2] for pb in photbands}
    fitsym = {pb:WAVE_d[pb][-2] for pb in photbands}
    fitcol = {pb:WAVE_d[pb][-1] for pb in photbands}
    filtercurves_i = photsys.filtercurves()
    filtercurves = {pb:filtercurves_i[pb] for pb in photbands}

    if bfdict['[Fe/H]'] >= 0.25:
        SEDfeh = 0.25
    elif bfdict['[Fe/H]'] <= -2.0:
        SEDfeh = -2.0
    else:
        SEDfeh = bfdict['[Fe/H]']

    if bfdict['Teff'] <= 3500.0:
        SEDTeff = 3500.0
    else:
        SEDTeff = bfdict['Teff']

    if bfdict['log(g)'] >= 5.0:
        SEDlogg = 5.0
    else:
        SEDlogg = bfdict['log(g)']

    spec_w,spec_f,_ = SB.get_star_spectrum(
        logt=np.log10(SEDTeff),logg=SEDlogg,feh=SEDfeh)

    to_cgs_i = lsun/(4.0 * np.pi * (pc*bfdict['Dist'])**2)
    nor = SB.normalize(logr=bfdict['log(R)'])*to_cgs_i
    spec_f = spec_f*nor
    spec_f = spec_f*(speedoflight/((spec_w*1E-8)**2.0))

    spec_f = np.nan_to_num(spec_f)
    spcond = spec_f > 1e-32
    spec_f = spec_f[spcond]
    spec_w = spec_w[spcond]

    extratio = ccm_curve(spec_w/10.0,bfdict['Av']/3.1)                    

    axSED.plot(spec_w/(1E+4),np.log10(spec_f/extratio),ls='-',lw=0.5,
        alpha=1.0,zorder=-1,c='C0')

    sedoutkeys = photdata.keys()
    sedpars = ([
        bfdict['Teff'],
        bfdict['log(g)'],
        bfdict['[Fe/H]'],
        bfdict['[a/Fe]'],
        bfdict['log(R)'],
        bfdict['Dist'],
        bfdict['Av'],
        ])
    sed = GM.genphot(sedpars)
    modmag = [sed[kk] for kk in sedoutkeys]

    # split out data into phot and error dict
    initphot = {kk:photdata[kk][0] for kk in sedoutkeys if kk in photbands}
    initphoterr = {kk:photdata[kk][1] for kk in sedoutkeys if kk in photbands}

    obswave   = np.array([WAVE[kk] for kk in sedoutkeys])
    fitsym    = np.array([fitsym[kk] for kk in sedoutkeys])
    fitcol    = np.array([fitcol[kk] for kk in sedoutkeys])
    fc        = [filtercurves[kk] for kk in sedoutkeys]
    obsmag    = np.array([initphot[kk] for kk in sedoutkeys if kk in photbands])
    obsmagerr = np.array([initphoterr[kk] for kk in sedoutkeys if kk in photbands])
    # modmag    = np.array([initphot[kk] for kk in sedoutkeys if kk in photbands])
    # modmag    = np.array([sedout[kk] for kk in sedoutkeys])
    obsflux_i = np.array([zeropts[kk]*10.0**(initphot[kk]/-2.5) for kk in sedoutkeys if kk in photbands])
    obsflux   = [x*(jansky_cgs)*(speedoflight/((lamb*1E-8)**2.0)) for x,lamb in zip(obsflux_i,obswave)]
    modflux_i = np.array([zeropts[kk]*10.0**(x/-2.5) for x,kk in zip(modmag,sedoutkeys)])
    modflux   = [x*(jansky_cgs)*(speedoflight/((lamb*1E-8)**2.0)) for x,lamb in zip(modflux_i,obswave)]

    # plot the observed SED and MAGS
    minobsflx = np.inf
    maxobsflx = -np.inf
    for w,f,mod,s,clr in zip(obswave,obsflux,modflux,fitsym,fitcol):
        if np.log10(f) > -30.0:
            axSED.scatter(w/1E+4,np.log10(mod),marker=s,c='C0',zorder=0,s=20)
            axSED.scatter(w/1E+4,np.log10(f),marker=s,c='k',zorder=1,s=10)
            if np.log10(f) < minobsflx:
                 minobsflx = np.log10(f)
            if np.log10(f) > maxobsflx:
                 maxobsflx = np.log10(f)

    # for w,m,me,mod,s,clr in zip(obswave,obsmag,obsmagerr,modmag,fitsym,fitcol):
    #     if np.abs(m-mod)/me > 5.0:
    #         me = np.abs(m-mod)
    #     if (m < 30) & (m > -30):
    #         axMAG.scatter(w/1E+4,mod,marker=s,c='C0',zorder=-1,s=20)
    #         axMAG.errorbar(w/1E+4,m,yerr=me,ls='',marker=',',c='k',zorder=0)
    #         axMAG.scatter(w/1E+4,m,marker=s,c='k',zorder=1,s=10)

    # plot filter curves
    # for fc_i,clr in zip(fc,fitcol):
    #     trans_i = 0.25*float(fc_i['trans'])*(0.9*float(maxobsflx)-1.1*float(minobsflx))+1.1*float(minobsflx)
    #     axSED.plot(fc_i['wave']/1E+4,trans_i,ls='-',lw=0.5,c=clr,alpha=1.0)

    axSED.set_ylim(1.1*minobsflx,0.9*maxobsflx)

    axSED.set_xlim([0.25,6.0])
    axSED.set_xscale('log')

    # axMAG.set_xlim([0.25,6.0])
    # axMAG.set_xscale('log')

    # axMAG.set_ylim(axMAG.get_ylim()[::-1])

    axSED.set_xticks([0.3,0.5,0.7,1.0,3,5])
    axSED.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    # axMAG.set_xticks([0.3,0.5,0.7,1,3,5])
    # axMAG.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    axSED.set_ylabel(r'log(F$_{\lambda}$) [erg s$^{-1}$ cm$^{-2}$]')

    axSED.set_xlabel(r'$\lambda$ [$\mu$m]')
    # axMAG.set_ylabel('mag')

    axSED.yaxis.tick_right()
    # axMAG.yaxis.tick_right()
    axSED.yaxis.set_label_position('right')
    # axMAG.yaxis.set_label_position('right')
    # axSED.set_xticklabels([])

def mkkiel(axkiel,samples):

    samples = samples.filled(np.nan)

    samples = samples[(samples['Pr'] > 1E-10)&(~np.isnan(samples['Pr']))&(np.isfinite(samples['Pr']))] # fixed bug of NaN params
    samples['Pr'] = samples['Pr'] / np.nansum(samples['Pr'])



    Teffbf = quantile(samples['Teff'],  [0.0001,0.16,0.5,0.84,0.9999],weights=samples['Pr'])
    loggbf = quantile(samples['log(g)'],[0.0001,0.16,0.5,0.84,0.9999],weights=samples['Pr'])
    # fehbf  = quantile(samples['[Fe/H]'],[0.0001,0.16,0.5,0.84,0.9999],weights=samples['Pr'])
    # afebf  = quantile(samples['[a/Fe]'],[0.0001,0.16,0.5,0.84,0.9999],weights=samples['Pr'])

    # kielstr = (
    #     'Teff   = {0:.0f} +{1:.0f}/-{2:.0f}\n'.format(Teffbf[2],Teffbf[3]-Teffbf[2],Teffbf[2]-Teffbf[1]) + 
    #     'log(g) = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(loggbf[2],loggbf[3]-loggbf[2],loggbf[2]-loggbf[1]) + 
    #     '[Fe/H] = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(fehbf[2],fehbf[3]-fehbf[2],fehbf[2]-fehbf[1]) + 
    #     '[a/Fe] = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(afebf[2],afebf[3]-afebf[2],afebf[2]-afebf[1])
    #     )
    # axkiel.text(0.01,0.0,kielstr,
    #     horizontalalignment='left',verticalalignment='bottom', 
    #     transform=axkiel.transAxes,fontsize=8)

    data = np.vstack([samples['Teff'], samples['log(g)']])
    print(np.isnan(data).any())
    print(data)
    print(np.isnan(samples['Pr']).any())
    print(samples['Pr'])
    kde = gaussian_kde(data,weights=samples['Pr'])

    # evaluate on a regular grid
    xgrid = np.linspace(Teffbf[0]-100.0,Teffbf[-1]+100.0, 100)
    ygrid = np.linspace(loggbf[0]-0.25,loggbf[-1]+0.25, 100)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))

    Z = Z.reshape(Xgrid.shape)
    Zcl = scoreatpercentile(Z,95)
    cond = Z < Zcl
    Z[cond] = np.nan

    axkiel.imshow(Z,
               origin='lower', aspect='auto',
               extent=([
                Teffbf[0]-50.0,Teffbf[-1]+50.0,
                loggbf[0]-0.1,loggbf[-1]+0.1]),
               # cmap=sns.light_palette('blue',as_cmap=True,reverse=False,n_colors=25),
               cmap='BrBG',
               alpha=0.75)


    # norm = plt.Normalize(1, 800)
    # cmap = ListedColormap(['C0', 'C1', 'C2','C3'])
    cmap = ListedColormap(['steelblue','forestgreen','plum'])
    norm = BoundaryNorm([203, 405, 606, 808], cmap.N)

    for samp_i in np.random.choice(samples,50,p=samples['Pr'],replace=False):
        Teffarr = []
        loggarr = []
        for eep_i in range(204,800,1):
            MISTpred = GMIST.getMIST(
                eep=eep_i,
                mass=samp_i['initial_Mass'],
                feh=samp_i['initial_[Fe/H]'],
                afe=samp_i['initial_[a/Fe]'],
                verbose=False)
            if MISTpred is None:
                continue
            MISTdict = ({
                kk:pp for kk,pp in zip(
                    GMIST.modpararr,MISTpred)
                })
            Teffarr.append(10.0**MISTdict['log(Teff)'])
            loggarr.append(MISTdict['log(g)'])

        points = np.array([Teffarr, loggarr]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(np.arange(204,800,1))
        lc.set_linewidth(2)
        lc.set_zorder(-1)
        lc.set_alpha(0.25)
        line = axkiel.add_collection(lc)

        # axkiel.plot(Teffarr,loggarr,ls='-',lw=2.0,c='C0',alpha=0.25,zorder=-1)

    # axkiel.set_xlim([Teffbf[-1]+50.0,Teffbf[0]-50.0])
    # axkiel.set_ylim([loggbf[-1]+0.1,loggbf[0]-0.1])
    axkiel.set_xlim([Teffbf[2]+0.1 * Teffbf[2],Teffbf[2]-0.1 * Teffbf[2]])
    axkiel.set_ylim([min([loggbf[2]+0.25 * loggbf[2],5.5]),max([loggbf[2]-0.50 * loggbf[2],-1])])
    axkiel.set_xlabel(r'$T_{eff}$ [K]')
    axkiel.set_ylabel('log(g)')

def run(index=None,GaiaID=None,version='VX',verbose=False,catalog = None,
    acat_id = None):
    # read in SEGUE data
    # SD = SegueData(survey=survey,tileID=tileID,mjd=mjd)


    if GaiaID is not None:
        data = getdata(GaiaID=GaiaID)
    if acat_id is not None:
        data = getdata(acat_id = acat_id)
    try:
        assert data['phot']
    except AssertionError:
        print('Warning: User did not pass a valid selection criteria')
        raise

    if catalog is None:
        catalog = 'solo'

    if verbose:
        print('... Phot Data Prep')

    phot = {}
    usedphot = {}
    for pp in photbands.keys():
        if ((data['phot'][pp] > 5.0) &
            (data['phot'][pp] < 90.0) & 
            (np.abs(data['phot'][pp+'_ERR']) < 90.0) & 
            ~np.isnan(data['phot'][pp]) & 
            ~np.isnan(data['phot'][pp+'_ERR'])
            ):
            filtersys = photbands[pp].split('_')[0]
            filtername = photbands[pp].split('_')[1]

            if filtersys == 'PS':
                if data['phot'][pp] <= 14.0:
                    continue
                photdata = data['phot'][pp]
                photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == '2MASS':
                if data['phot'][pp] <= 5.0:
                    continue
                photdata = data['phot'][pp]
                photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
            elif (filtersys == 'WISE') or (filtersys == 'UNWISE'):
                if data['phot'][pp] <= 8.0:
                    continue
                photdata = data['phot'][pp]
                photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
            elif filtersys == 'SDSS':
                if data['phot'][pp] <= 12.0:
                    continue
                if filtername == 'u':
                    photdata = data['phot'][pp] - 0.04
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
                elif filtername == 'i':
                    photdata = data['phot'][pp] + 0.02
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.02**2.0)))
                else:
                    photdata = data['phot'][pp]
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.02**2.0)))
            elif filtersys == 'GaiaEDR3':
                if filtername == 'G':
                    photdata = data['phot'][pp+'_CORRECTED']
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'BP':
                    photdata = data['phot'][pp]
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
                if filtername == 'RP':
                    photdata = data['phot'][pp]
                    photerr = np.sqrt((data['phot'][pp+'_ERR']**2.0)+((0.05**2.0)))
            else:
                photdata = data['phot'][pp]
                photerr = data['phot'][pp+'_ERR']

            usedphot[photbands[pp]] = pp
            phot[photbands[pp]] = [float(photdata),float(photerr)]

    if verbose:
        print('... Spec Data Prep')

    specdata = data['spec']
    spec = {}
    # cond = specdata[2] != 0.0
    spec['WAVE']   = specdata['wave']
    spec['FLUX']   = specdata['flux']
    spec['E_FLUX'] = 1.0/np.sqrt(specdata['ivar'])
    spec['WRESL']  = specdata['wresl']

    # cond = np.isfinite(spec['FLUX']) & np.isfinite(spec['E_FLUX']) & (spec['LSF'] > 0.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['LSF']    = spec['LSF'][cond]

    # create the WRESL array
    # spec['WRESL'] = spec['LSF'] #(spec['WAVE'] * spec['LSF']) / speedoflight

    # cond = (spec['WAVE'] > 3850.0) & (spec['WAVE'] < 8900.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

    # cond = (spec['WAVE'] > 4000.0) & (spec['WAVE'] < 7000.0)
    # cond = (spec['WAVE'] > 4455.0) & (spec['WAVE'] < 5645.0)
    # cond = (spec['WAVE'] > 5000.0) & (spec['WAVE'] < 5500.0)    
    # cond = (spec['WAVE'] > 4750.0) & (spec['WAVE'] < 5500.0)
    # spec['WAVE']   = spec['WAVE'][cond]
    # spec['FLUX']   = spec['FLUX'][cond]
    # spec['E_FLUX'] = spec['E_FLUX'][cond]
    # spec['WRESL']  = spec['WRESL'][cond]

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

    # spec['WAVE'] = airtovacuum(spec['WAVE'])

    samplefile = 'mage_{GAIAID}_{MJD}_{VER}_samp.dat'.format(
            GAIAID=data['phot']['GAIAEDR3_ID'],
            MJD=data['phot']['date'],
            VER=version)

    samplefile = '{OUTDIR}samples/{CATALOG}/{VER}/{SAMPLEFILE}'.format(
            OUTDIR=outdir,
            SAMPLEFILE=samplefile,
            CATALOG = catalog,
            VER=version)

    samplefile_gz = samplefile + '.gz'

    if verbose:
        print('... Read in Samples from {}'.format(samplefile))

    try:
        try:
            samp = Table.read(samplefile,format='ascii')
        except:
            samp = Table.read(samplefile_gz,format='ascii')
        if samp['delta(log(z))'][-1] > 0.01:
            print('COMPMOD: delta(log(z)) DID NOT CONVERGE!')
            return
    except:
        print('COMPMOD: Problem with reading sample file!')
        return

    samp['Pr'] = np.exp(samp['log(wt)']-samp['log(z)'][-1])
    samp = samp[samp['Pr'] > 1e-10]
    samp['Pr'] = samp['Pr'] / np.sum(samp['Pr'])
    samp.remove_column('Inst_R')
    sampval = np.lib.recfunctions.structured_to_unstructured(np.array(samp))
    print(np.isnan(sampval).any())
    samp = samp.filled(np.nan)

    maxlike = np.argmax(samp['log(lk)'])

    fitpars = ([
        x for x in samp.keys() if x not in 
        ['Iter','Agewgt','log(lk)','log(vol)','log(wt)','h','nc','log(z)','delta(log(z))','Pr']
        ])

    bf = {}
    for ff in fitpars:
        bf[ff] = samp[ff][maxlike]

    ii = 0
    for x in fitpars:
        if 'pc_' in x:
            ii += 1
    if verbose:
        print('... Calcluate BF model')
    bfpcpars = [bf['pc_{}'.format(x)] for x in range(ii)]
    # bfpc  = polycalc(bfpcpars,spec['WAVE'])
    bfpars = [bf['Teff'],bf['log(g)'],bf['[Fe/H]'],bf['[a/Fe]'],bf['Vrad'],bf['Vrot'],np.nan,spec['WRESL']]
    bfspecpars = bfpars + bfpcpars

    bfmod = GM.genspec(bfspecpars,outwave=spec['WAVE'],modpoly=True)
    bfflux = bfmod[1]*medflux

    if verbose:
        print('... Making Plots')

    fig = plt.figure(figsize=(11,8.5),constrained_layout=True)
    gs = fig.add_gridspec(6,4)
    axspec = fig.add_subplot(gs[:3,:-1])
    axres  = fig.add_subplot(gs[3,:-1])

    axSED = fig.add_subplot(gs[-2:,-2:])
    # axMAG = fig.add_subplot(gs[-1,-2:])

    axkiel = fig.add_subplot(gs[-2:,:-2])

    if verbose:
        print('... Making Kiel')

    mkkiel(axkiel,samp)

    if verbose:
        print('... Making SED')

    mksed(axSED,samp,phot,bf)

    Teffbf  = quantile(samp['Teff'],  [0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    loggbf  = quantile(samp['log(g)'],[0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    fehbf   = quantile(samp['[Fe/H]'],[0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    afebf   = quantile(samp['[a/Fe]'],[0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    vradbf  = quantile(samp['Vrad'],  [0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    vrotbf  = quantile(samp['Vrot'],  [0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])
    distbf  = quantile(samp['Dist'],  [0.0001,0.16,0.5,0.84,0.9999],weights=samp['Pr'])

    parstr = (
        'GaiaEDR3 ID = {0} \n'.format(data['phot']['GAIAEDR3_ID']) +
        'GaiaEDR3 Para = {0:.3f} +/- {1:.3f} \n'.format(data['phot']['GAIAEDR3_PARALLAX_CORRECTED'],data['phot']['GAIAEDR3_PARALLAX_ERROR']) + 
        'Teff   = {0:.0f} +{1:.0f}/-{2:.0f} K\n'.format(Teffbf[2],Teffbf[3]-Teffbf[2],Teffbf[2]-Teffbf[1]) + 
        'log(g) = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(loggbf[2],loggbf[3]-loggbf[2],loggbf[2]-loggbf[1]) + 
        '[Fe/H] = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(fehbf[2],fehbf[3]-fehbf[2],fehbf[2]-fehbf[1]) + 
        '[a/Fe] = {0:.3f} +{1:.3f}/-{2:.3f}\n'.format(afebf[2],afebf[3]-afebf[2],afebf[2]-afebf[1]) +
        'Vrad   = {0:.3f} +{1:.3f}/-{2:.3f} km/s\n'.format(vradbf[2],vradbf[3]-vradbf[2],vradbf[2]-vradbf[1]) +
        'Vstar  = {0:.3f} +{1:.3f}/-{2:.3f} km/s\n'.format(vrotbf[2],vrotbf[3]-vrotbf[2],vrotbf[2]-vrotbf[1]) +
        'Dist   = {0:.3f} +{1:.3f}/-{2:.3f} kpc\n'.format(distbf[2]/1000.0,(distbf[3]-distbf[2])/1000.0,(distbf[2]-distbf[1])/1000.0)
        )
    plt.figtext(0.75,0.35,parstr,
        horizontalalignment='left',verticalalignment='bottom', 
        fontsize=8)

    # zoomwr  = [[4000,4500],[4750.0,5000.0],[5150,5300],[5800,6000],[6450,6650],[6675,6750]]
    # zoomlab = [r'H$\delta$,H$\gamma$',r'H$\beta$','Mg b','Na D',r'H$\alpha$','Li I']

    # zoomwr  = [[4750.0,4850.0],[4825.0,5000.0],[5150,5300]]
    # zoomlab = [r'Fe-lines',r'H$\beta$','Mg b']

    zoomwr  = [[5150,5300]]
    zoomlab = ['Mg b']


    # zoomwr  = [[5000.0,5150.0],[5150,5300]]
    # zoomlab = [r'Fe-lines','Mg b']

    axline = ({x:{
        'ax':fig.add_subplot(gs[ii,-1]),'xran':xran} for x,ii,xran 
        in zip(
            zoomlab,
            range(len(zoomlab)),
            zoomwr)
        })

    axspec.plot(spec['WAVE']/10.0,spec['FLUX']*(speedoflight_nms/(spec['WAVE']**2.0)),ls='-',lw=0.5,c='k',zorder=0)
    # axspec.plot(spec['WAVE']/10.0,bfpc*medflux,ls='-',lw=0.5,c='C1',zorder=1)
    axspec.plot(spec['WAVE']/10.0,bfflux*(speedoflight_nms/(spec['WAVE']**2.0)),ls='-',lw=1.0,c='C0',zorder=2)

    axres.plot(spec['WAVE']/10.0,(bfflux-spec['FLUX'])/spec['E_FLUX'],ls='-',lw=0.5,c='C0',zorder=0)
    axres.axhline(y=0.0,lw=1.0,c='C3',alpha=0.5,ls='-',zorder=1)
    axres.axhline(y=-1.0,lw=1.0,c='C3',alpha=0.5,ls='--',zorder=1)
    axres.axhline(y=1.0,lw=1.0,c='C3',alpha=0.5,ls='--',zorder=1)

    # axspec.axvspan(585.0,595.0,color='k',alpha=0.1,zorder=-1)
    # axres.axvspan(585.0,595.0,color='k',alpha=0.1,zorder=-1)

    axspec.set_ylim(0.9*np.nanmin(bfflux*(speedoflight_nms/(spec['WAVE']**2.0))),1.1*np.nanmax(bfflux*(speedoflight_nms/(spec['WAVE']**2.0))))
    # axres.set_ylim(quantile(bfflux-spec['FLUX'],0.001)[0],quantile(bfflux-spec['FLUX'],0.999)[0])
    # axres.set_ylim(-5,5)

    axspec.set_ylabel(r'f$_{\lambda}$ [cts s$^{-1}$ cm$^{-2}$ nm$^{-1}$]')
    axspec.set_xlim(spec['WAVE'].min()/10.0,spec['WAVE'].max()/10.0)
    axspec.set_xticklabels([])

    axres.set_xlabel('Wavelength [nm]')
    axres.set_ylabel(r'$\chi$')
    axres.set_xlim(spec['WAVE'].min()/10.0,spec['WAVE'].max()/10.0)

    anno_opts = dict(xy=(0.85, 0.1), xycoords='axes fraction',va='center', ha='center')
    for x in axline.keys():
        cond = (spec['WAVE'] > axline[x]['xran'][0]) & (spec['WAVE'] < axline[x]['xran'][1])
        axline[x]['ax'].plot(spec['WAVE'][cond]/10.0,spec['FLUX'][cond],ls='-',lw=0.5,c='k',zorder=0)
        # axline[x]['ax'].plot(spec['WAVE'][cond]/10.0,bfpc[cond]*medflux,ls='-',lw=0.75,c='C1',zorder=1)
        axline[x]['ax'].plot(spec['WAVE'][cond]/10.0,bfflux[cond],ls='-',lw=0.75,c='C0',zorder=2)
        axline[x]['ax'].set_xlim([axline[x]['xran'][0]/10.0,axline[x]['xran'][1]/10.0])
        axline[x]['ax'].set_yticklabels([])
        axline[x]['ax'].annotate(x, **anno_opts)
        axline[x]['ax'].set_ylim(0.75*bfflux[cond].min(),1.25*bfflux[cond].max())


    compmodfile = '{OUTDIR}{CATALOG}/{VER}/mage_{GAIAID}__{MJD}_{VER}_compmod.png'.format(
        OUTDIR=outdir + 'plots/',
        GAIAID=data['phot']['GAIAEDR3_ID'],
        MJD=data['phot']['date'],
        VER=version,
        CATALOG = catalog)

    fig.savefig(compmodfile,dpi=200)

    plt.close(fig)

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--survey',help='SEGUE Survey ID',type=str,
    #     choices=['SEGUE','SEGUE_clusters'],default='SEGUE')
    # parser.add_argument('--tileID',help='tile ID number',type=int,default=1660)
    # parser.add_argument('--index',help='Index of star in acat',type=int,default=None)
    # parser.add_argument('--GaiaID',help='Gaia EDR3 ID of star',type=int,default=None)
    # parser.add_argument('--FiberID',help='Fiber ID of star',type=int,default=None)
    # parser.add_argument('--mjd',help='MJD of plate to run',type=int,default=None)
    # parser.add_argument("--version", "-v", help="run version",type=str,default='VX')
    # args = parser.parse_args()
    # run(
    #     survey=args.survey,
    #     tileID=args.tileID,
    #     index=args.index,
    #     GaiaID=args.GaiaID,
    #     FiberID=args.FiberID,
    #     mjd=args.mjd,
    #     version=args.version)

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
        version=args.version,)
