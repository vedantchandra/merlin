from astropy.table import Table
import numpy as np
import sys,itertools,argparse,glob
from quantiles import quantile

import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines

from scipy.stats import norm, truncnorm, expon, reciprocal,beta

from getdata import getdata

from minesweeper.advancedpriors import AdvancedPriors


# specNN = '{}/conroy_lab/pacargile/ThePayne/Hecto_FAL/lowres/YSTANN_4000_7000_spec.h5'.format(holypath)
# contNN = '{}/conroy_lab/pacargile/ThePayne/Hecto_FAL/lowres/YSTANN_4000_7000_cont.h5'.format(holypath)
# specNN = '{}/conroy_lab/pacargile/ThePayne/train/optfal/v256/modV0_spec_LinNet_R5K_WL445_565.h5'.format(holypath)
# contNN = '{}/conroy_lab/pacargile/ThePayne/train/optfal/v256/modV0_cont_LinNet_R12K_WL445_565.h5'.format(holypath)
# NNtype = 'LinNet'
# photNN = '{}/conroy_lab/pacargile/ThePayne/SED/VARRV/'.format(holypath)
# SBlib = '{}/conroy_lab/pacargile/CKC/ckc_R500.h5'.format(holypath)
# MISTgrid = '{}/conroy_lab/pacargile/MIST/MIST_2.0_spot_EEPtrk_small.h5'.format(holypath)
# datadir = '/n/holystore01/LABS/conroy_lab/Lab/SEGUE/data/'
# outdir  = '/n/holyscratch01/conroy_lab/pacargile/SEGUE/'

datadir = '/n/holyscratch01/conroy_lab/vchandra/sdss5/'
specNN = datadir + 'ms/NN/modV0_spec_LinNet_R5K_WL445_565.h5'
contNN = datadir + 'ms/NN/modV0_cont_LinNet_R12K_WL445_565.h5'
photNN = datadir + 'ms/VARRV/'
SBlib = datadir + 'ms/CKC/ckc_R500.h5'
MISTgrid = datadir + 'ms/MIST_2.0_spot_EEPtrk_small.h5'
outdir = datadir
NNtype = 'LinNet'


def run(index=None,GaiaID=None,version='VX',catalog=None,
    acat_id = None):
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

    if catalog is None:
        catalog = 'solo'

    photcat = data['phot']

    samplefile = '{OUTDIR}/{CATALOG}/{VER}/mwm_gaiaID_{GAIAID}_fieldID_{FIELDID}_mjd_{MJD}_catID_{CATID}_{VER}_samp.dat'.format(
            OUTDIR = outdir + 'samples/',
            FIELDID=data['phot']['FIELD'],
            GAIAID=data['phot']['GAIAEDR3_ID'],
            CATID=data['phot']['CATALOGID'],
            MJD=data['phot']['MJD'],
            VER=version,
            CATALOG = catalog)
    samplefile_gz = samplefile + '.gz'

    # try:
    try:
        samp = Table.read(samplefile,format='ascii')
    except:
        samp = Table.read(samplefile_gz,format='ascii')
    if samp['delta(log(z))'][-1] > 0.01:
        print('CORNER: delta(log(z)) DID NOT CONVERGE!')
        return
    # except:
    #     print('CORNER: Problem with reading sample file!')
    #     return

    samp['Pr'] = np.exp(samp['log(wt)']-samp['log(z)'][-1])
    samp = samp[samp['Pr'] > 1E-10]

    samp['Dist'] = samp['Dist']/1000.0
    if 'Para' not in samp.keys():
        samp['Para'] = 1.0/samp['Dist']
    if 'Age' not in samp.keys():
        samp['Age'] = 10.0**(samp['log(Age)']-9.0)

    pltpars = np.array([
        'Teff','log(g)','Vrad','Vrot','pc_0','pc_1','pc_2','pc_3','Dist',
        'Para','Av','EEP','Age','initial_[Fe/H]','initial_[a/Fe]','initial_Mass','log(L)'
        ])

    fitpars = ([
        x for x in samp.keys() if x not in 
        ['Iter','Agewgt','log(lk)','log(vol)','log(wt)','h','nc','log(z)','delta(log(z))','Pr']
        ])

    pardict = {}
    for ff in fitpars:
        pardict[ff] = quantile(samp[ff],[0.0001,0.01,0.16,0.5,0.84,0.99,0.9999],weights=samp['Pr'])

    gaia_parallax = [photcat['GAIAEDR3_PARALLAX_CORRECTED'],photcat['GAIAEDR3_PARALLAX_ERROR']]
    lb_coords = [photcat['L'],photcat['B']]

    parstring = (
        'Run index: {} '.format(photcat['ACAT_ID'])+
        'Field ID: {}\n'.format(photcat['FIELD'])+
        'GaiaEDR3 ID: {}\n'.format(photcat['GAIAEDR3_ID'])+
        'RA: {:n} '.format(photcat['RA'])+'Dec: {:n}\n'.format(photcat['DEC'])+
        'L: {:n} '.format(photcat['L'])+'B: {:n}\n'.format(photcat['B'])+
        'PS_g = {0:n} +/- {1:n}\n'.format(photcat['PS_G'],photcat['PS_G_ERR'])+
        'GaiaEDR3_G = {0:n} +/- {1:n}\n'.format(photcat['GAIAEDR3_G_CORRECTED'],photcat['GAIAEDR3_G_ERR'])+
        'GaiaEDR3_Parallax = {0:n} +/- {1:n}\n'.format(photcat['GAIAEDR3_PARALLAX_CORRECTED'],photcat['GAIAEDR3_PARALLAX_ERROR'])+
        'SNR = {:n}\n'.format(photcat['SN_MEDIAN_ALL'])
        )

    for ff in fitpars:
        quant = quantile(samp[ff],[0.5,0.16,0.84],weights=samp['Pr'])
        parstring = parstring + '{0} = {1:n} +{2:n}/-{3:n} \n'.format(ff,quant[0],quant[2]-quant[0],quant[0]-quant[1])

    parind = np.array(range(len(pltpars)))

    AP = AdvancedPriors()

    fig = plt.figure(figsize=(15,15))
    plt.figtext(0.6,0.45,parstring,fontsize=12)
    gs = gridspec.GridSpec(len(pltpars),len(pltpars))
    gs.update(wspace=0.05,hspace=0.05)

    nbins = 35

    for kk in itertools.product(pltpars,pltpars):
        kkind1 = parind[pltpars == kk[0]][0]
        kkind2 = parind[pltpars == kk[1]][0]
        ax = fig.add_subplot(gs[kkind1,kkind2])

        if kkind1 < kkind2:
            ax.set_axis_off()
            continue

        xarr_range = [pardict[kk[0]][0],pardict[kk[0]][-1]]
        yarr_range = [pardict[kk[1]][0],pardict[kk[1]][-1]]


        if kk[0] == kk[1]:
            n,bins,_ = ax.hist(
            samp[kk[0]],
            bins=nbins,
            histtype='step',
            linewidth=1.5,
            density=True,
            range=xarr_range,
            weights=samp['Pr'],
            )

            if (kk[0] == 'pc_0'):
                xarr = np.linspace(0.5,2.0,200)
                yarr = n.max()*np.ones(len(xarr))
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')

            pcterms = {'pc_1':[0.0,0.5],'pc_2':[0.0,0.5],'pc_3':[0.0,0.25]}
            if kk[0] in ['pc_1','pc_2','pc_3']:
                # xarr = np.linspace(xarr_range[0],xarr_range[1],200)
                xarr = np.linspace(-0.5,0.5,200)
                yarr = np.exp( -0.5*((xarr-pcterms[kk[0]][0])**2.0)/(pcterms[kk[0]][1]**2.0))
                yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')

            if (kk[0] == 'Para'):
                # xarr = np.linspace(xarr_range[0],xarr_range[1],200)
                xarr = np.logspace(-4,1,500)
                yarr = np.exp( -0.5*((xarr-gaia_parallax[0])**2.0)/(gaia_parallax[1]**2.0) )
                yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')

            if (kk[0] == 'Vrot'):
                xarr = np.logspace(0,250,500)
                a = 1.05
                b = 1.5
                loc = 0.0
                scale = 250.0
                yarr = beta.pdf(xarr,a,b,loc=loc,scale=scale)
                yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')

            if (kk[0] == 'initial_Mass'):
                xarr = np.linspace(xarr_range[0],xarr_range[1],200)
                # xarr = np.linspace(0.25,3.5,100)
                yarr = np.exp(AP.imf_lnprior(xarr))
                yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')


            if (kk[0] == 'Dist'):
                if gaia_parallax[0] > 0.0:
                    if 3.0*gaia_parallax[1] < gaia_parallax[0]:
                        maxdist = 1.0/(gaia_parallax[0]-3.0*gaia_parallax[1])
                    else:
                        maxdist = 200.0

                    mindist = 1.0/(gaia_parallax[0]+3.0*gaia_parallax[1])
                else:
                    mindist = max([1.0,1.0/(3.0*gaia_parallax[1])])
                    maxdist = 200.0

                # xarr = np.linspace(self.SAMPLEStab['Dist'].min(),self.SAMPLEStab['Dist'].max(),200)
                xarr = np.linspace(mindist,maxdist,200)
                yarr = np.exp(AP.gal_lnprior(xarr,coords=lb_coords))
    
                if isinstance(yarr,float):
                    yarr = 0.0*np.ones_like(xarr)
                else:
                    yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')

            if (kk[0] == 'Age'):
                # xarr = np.linspace(self.SAMPLEStab['Dist'].min(),self.SAMPLEStab['Dist'].max(),200)
                xarr = np.linspace(xarr_range[0],xarr_range[1],200)

                lnp_dist,comp = AP.gal_lnprior(pardict['Dist'][3],coords=lb_coords,return_components=True)
                # Compute component membership probabilities.
                logp_thin  = comp['number_density'][0]
                logp_thick = comp['number_density'][1]
                logp_halo  = comp['number_density'][2]

                lnprior_thin = logp_thin - lnp_dist
                lnprior_thick = logp_thick - lnp_dist
                lnprior_halo = logp_halo - lnp_dist

                galmodpars = ({
                    'thin': {'min':1.0,'max':14.0},
                    'thick':{'min':6.0,'max':14.0,'mean':10.0,'sigma':2.0},
                    'halo': {'min':8.0,'max':14.0,'mean':12.0,'sigma':2.0},
                    })

                yarr = np.exp(AP.age_lnprior(xarr,
                    lnp_thin=lnprior_thin,
                    lnp_thick=lnprior_thick,
                    lnp_halo=lnprior_halo,**galmodpars))

                yarr = n.max()*(yarr-yarr.min())/(yarr.max()-yarr.min())
                ax.plot(xarr,yarr,ls='-',lw=1.0,c='green')


            ax.set_xlim(xarr_range[0],xarr_range[1])
            ylimtmp = ax.get_ylim()
            ax.set_ylim(ylimtmp[0],1.25*ylimtmp[1])

            ax.set_yticks([])
            if kk[0] != pltpars[-1]:
                 ax.set_xticks([])
            else:
                 ax.set_xlabel(kk[0])

        else:
            ax.hist2d(
                samp[kk[1]],
                samp[kk[0]],
                bins=nbins,
                cmap='Blues',
                range=[yarr_range,xarr_range],
                weights=samp['Pr'],
                )

            ax.set_xlim(yarr_range[0],yarr_range[1])
            ax.set_ylim(xarr_range[0],xarr_range[1])

        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        [l.set_rotation(45) for l in ax.get_xticklabels()]
        [l.set_fontsize(6) for l in ax.get_xticklabels()]
        [l.set_fontsize(6) for l in ax.get_yticklabels()]

        labelcol = 'k'

        try:
            isfc = ax.get_subplotspec().is_first_col()
            islc = ax.get_subplotspec().is_last_col()
            isfr = ax.get_subplotspec().is_first_row()
            islr = ax.get_subplotspec().is_last_row()
        except AttributeError:
            isfc = ax.is_first_col()
            islc = ax.is_last_col()
            isfr = ax.is_first_row()
            islr = ax.is_last_row()

        if not isfc:
            ax.set_yticks([])
        elif isfc & isfr:
            ax.set_yticks([])
        elif kk[0] == pltpars[0]:
            pass
        else:
            if 'initial' in kk[0]:
                 ax.set_ylabel(kk[0].split('_')[1]+r'$_{i}$')
            elif kk[0] == '[a/Fe]':
                 ax.set_ylabel('['+r'$\alpha$'+'/Fe]')
            elif kk[0] == 'Vrad':
                 ax.set_ylabel(r'V$_{rad}$')
            elif kk[0] == 'Vrot':
                 ax.set_ylabel(r'V$_{\bigstar}$')
            elif kk[0] == 'Para':
                 ax.set_ylabel(r'$\pi$')
            elif 'pc_' in kk[0]:
                ax.set_ylabel(r'pc$_{0}$'.format(kk[0].split('_')[1]))
            else:
                 ax.set_ylabel(kk[0])

        if not islr:
            ax.set_xticks([])
        else:
            if 'initial' in kk[1]:
                 ax.set_xlabel(kk[1].split('_')[1]+r'$_{i}$')
            elif kk[1] == '[a/Fe]':
                 ax.set_xlabel('['+r'$\alpha$'+'/Fe]')
            elif kk[1] == 'Teff':
                 ax.set_xlabel(r'T$_{eff}$'+'\n[K]')
            elif kk[1] == 'Dist':
                 ax.set_xlabel('Dist.'+'\n[kpc]')
            elif kk[1] == 'Age':
                 ax.set_xlabel('Age'+'\n[Gyr]')
            elif kk[1] == 'Vrad':
                 ax.set_xlabel(r'V$_{rad}$'+'\n[km/s]')
            elif kk[1] == 'Vrot':
                 ax.set_xlabel(r'V$_{\bigstar}$'+'\n[km/s]')
            elif kk[1] == 'Para':
                 ax.set_xlabel(r'$\pi$'+'\n["]')
            elif 'pc_' in kk[1]:
                ax.set_xlabel(r'pc$_{0}$'.format(kk[1].split('_')[1]))
            else:
                 ax.set_xlabel(kk[1])

    fig.align_labels()

    cornerfile = '{OUTDIR}{CATALOG}/{VER}/mwm_gaiaID_{GAIAID}_fieldID_{FIELDID}_mjd_{MJD}_catID_{CATID}_{VER}_corner.png'.format(
        OUTDIR=outdir + 'plots/',
        FIELDID=data['phot']['FIELD'],
        GAIAID=data['phot']['GAIAEDR3_ID'],
        CATID=data['phot']['CATALOGID'],
        MJD=data['phot']['MJD'],
        VER=version,
        CATALOG=catalog)


    fig.savefig(cornerfile,dpi=200)

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