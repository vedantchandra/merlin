import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
from tqdm import tqdm

import numpy as np

import astropy.units as u
from astropy.table import Table

import argparse

outdir  = '/n/holyscratch01/conroy_lab/vchandra/mage/'

class mkrcat(object):
    """docstring for mkrcat"""
    def __init__(self, catalog, version, runtype, outfile):
        super(mkrcat, self).__init__()
        self.catalog = catalog
        self.version  = version
        self.runtype  = runtype
        
        self.outfile = outfile
        
        # add additional columns with nan's
        self.newcols = ([
            "EEP",
            "EEP_lerr",
            "EEP_uerr",
            "EEP_err",
            "initial_[Fe/H]",
            "initial_[Fe/H]_lerr",
            "initial_[Fe/H]_uerr",
            "initial_[Fe/H]_err",
            "initial_[a/Fe]",
            "initial_[a/Fe]_lerr",
            "initial_[a/Fe]_uerr",
            "initial_[a/Fe]_err",
            "initial_Mass",
            "initial_Mass_lerr",
            "initial_Mass_uerr",
            "initial_Mass_err",
            "pc_0",
            "pc_0_lerr",
            "pc_0_uerr",
            "pc_0_err",
            "pc_1",
            "pc_1_lerr",
            "pc_1_uerr",
            "pc_1_err",
            "pc_2",
            "pc_2_lerr",
            "pc_2_uerr",
            "pc_2_err",
            "pc_3",
            "pc_3_lerr",
            "pc_3_uerr",
            "pc_3_err",
            "Inst_R",
            "Inst_R_lerr",
            "Inst_R_uerr",
            "Inst_R_err",
            "Teff",
            "Teff_lerr",
            "Teff_uerr",
            "Teff_err",
            "log(g)",
            "log(g)_lerr",
            "log(g)_uerr",
            "log(g)_err",
            "log(R)",
            "log(R)_lerr",
            "log(R)_uerr",
            "log(R)_err",
            "[Fe/H]",
            "[Fe/H]_lerr",
            "[Fe/H]_uerr",
            "[Fe/H]_err",
            "[a/Fe]",
            "[a/Fe]_lerr",
            "[a/Fe]_uerr",
            "[a/Fe]_err",
            "Vrad",
            "Vrad_lerr",
            "Vrad_uerr",
            "Vrad_err",
            "Vrot",
            "Vrot_lerr",
            "Vrot_uerr",
            "Vrot_err",
            "Dist",
            "Dist_lerr",
            "Dist_uerr",
            "Dist_err",
            "Av",
            "Av_lerr",
            "Av_uerr",
            "Av_err",
            "log(Age)",
            "log(Age)_lerr",
            "log(Age)_uerr",
            "log(Age)_err",
            "Mass",
            "Mass_lerr",
            "Mass_uerr",
            "Mass_err",
            "log(L)",
            "log(L)_lerr",
            "log(L)_uerr",
            "log(L)_err",
            "Para",
            "Para_lerr",
            "Para_uerr",
            "Para_err",
            "Age",
            "Age_lerr",
            "Age_uerr",
            "Age_err",
            "lnZ",
            "lnL",
            "lnP",
            "chisq_spec",
            "nspecpix",
            "chisq_phot",
            "nbands",
            "R_gal",
            "R_gal_err",
            "X_gal",
            "X_gal_err",
            "Y_gal",
            "Y_gal_err",
            "Z_gal",
            "Z_gal_err",
            "Vx_gal",
            "Vx_gal_err",
            "Vy_gal",
            "Vy_gal_err",
            "Vz_gal",
            "Vz_gal_err",
            "Vr_gal",
            "Vr_gal_err",
            "Vphi_gal",
            "Vphi_gal_err",
            "Vtheta_gal",
            "Vtheta_gal_err",
            "V_tan",
            "V_tan_err",
            "V_gsr",
            "V_gsr_err",
            "Lx",
            "Lx_err",
            "Ly",
            "Ly_err",
            "Lz",
            "Lz_err",
            "Ltot",
            "Ltot_err",
            "E_kin_pot1",
            "E_kin_pot1_err",
            "E_pot_pot1",
            "E_pot_pot1_err",
            "E_tot_pot1",
            "E_tot_pot1_err",
            "circLz_pot1",
            "circLz_pot1_err",
            "circLtot_pot1",
            "circLtot_pot1_err",
            ])

    def filltab(self,rcat_i):
        mjd     = rcat_i['date']
        gaiaID  = rcat_i['GAIAEDR3_ID']

        parfile = '{OUTDIR}pars/{CATALOG}/{VER}/mage_{GAIAID}_{MJD}_{VER}.pars'.format(
            OUTDIR=outdir,
            GAIAID=gaiaID,
            MJD=mjd,
            VER=self.version,
            RUNTYPE=self.runtype,
            CATALOG=self.catalog)

        try:
            pars = Table.read(parfile,format='ascii')
        except:
            raise
            return

        #print(pars)

        print('... working on : {0}'.format(parfile))

        for nc in self.newcols:
            rcat_i[nc] = pars[nc][0]

        return rcat_i

        # cond = (
        #     (rcat_i['Vrot'] > 200.0) | 
        #     (rcat_i['SN_MEDIAN_ALL'] < 10.0) | 
        #     (rcat_i['initial_[Fe/H]'] < -4.0) |
        #     (rcat_i['Teff'] > 7000.0) |
        #     (rcat_i['chisq_spec']/rcat_i['nspecpix'] > 3.0) |
        #     (rcat_i['nbands'] < 5.0) |
        #     (rcat_i['GAIAEDR3_RUWE'] > 1.5)
        # )

        # if cond:
        #     rcat_i['FLAG'] = 1
        # else:
        #     rcat_i['FLAG'] = 0


    def calcSgr(self,outdict):
        # create Sgr coordinates
        #   outdict['Sgr_l'] = np.nan 
        #   outdict['Sgr_b'] = np.nan 
        #   outdict['Sgr_FLAG'] = np.int16(0)

        ra_deg = outdict['GAIAEDR3_RA']*u.degree
        dec_deg = outdict['GAIAEDR3_DEC']*u.degree

        arg1 = (-0.93595354*np.cos(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian)) 
            -0.31910658*np.sin(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian)) 
            + 0.14886895*np.sin(dec_deg.to(u.radian)))

        arg2 = (0.21215555*np.cos(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian))
            -0.84846291*np.sin(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian)) 
            -0.48487186*np.sin(dec_deg.to(u.radian)))

        outdict['Sgr_l'] = np.mod(360.0*u.degree-np.arctan2(arg1,arg2).to(u.degree),360.0*u.degree)

        outdict['Sgr_b'] = -1.0*u.degree*(
            np.arcsin(0.28103559*np.cos(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian))
                -0.42223415*np.sin(ra_deg.to(u.radian))*np.cos(dec_deg.to(u.radian))
                +0.86182209*np.sin(dec_deg.to(u.radian))).to(u.degree))

        outdict['Sgr_l'] = outdict['Sgr_l'].value
        outdict['Sgr_b'] = outdict['Sgr_b'].value

        outdict['Sgr_FLAG'] = np.zeros(len(outdict['Sgr_l']),dtype=np.int16)

          # add Sgr flag
        cond = (
               (outdict['logg'] < 3.5) & 
               ((outdict['Ly']/(1E+4)) < (-0.3*(outdict['Lz']/(1E+4))-0.25)) & 
               ((outdict['Ly']/(1E+4)) < ( 1.2*(outdict['Lz']/(1E+4))+0.18))
               )
        
        outdict['Sgr_FLAG'][cond] = np.int16(1)

        #   if isinstance(outdict['Sgr_FLAG'],str):
        #        outdict['Sgr_FLAG'] = np.int16(-1)

        return outdict
      
    def run(self):
        # copy acat to rcat
        acat = Table.read(outdir + 'catalogs/%s_acat.fits' % self.catalog,format='fits')
        rcat = acat.copy()
        del acat

        # rcat['FLAG'].name          = 'SEGUE_FLAG'
        # rcat['G_MAG'].name         = 'SEGUE_G_MAG'
        # rcat['TEFF_ADOP'].name     = 'SEGUE_TEFF'
        # rcat['TEFF_ADOP_UNC'].name = 'SEGUE_TEFF_ERR'
        # rcat['LOGG_ADOP'].name     = 'SEGUE_LOGG'
        # rcat['LOGG_ADOP_UNC'].name = 'SEGUE_LOGG_ERR'
        # rcat['FEH_ADOP'].name      = 'SEGUE_FEH'
        # rcat['FEH_ADOP_UNC'].name  = 'SEGUE_FEH_ERR'
        # rcat['AFE'].name           = 'SEGUE_AFE'
        # rcat['AFE_UNC'].name       = 'SEGUE_AFE_ERR'
        # rcat['DIST_DWARF'].name    = 'SEGUE_DIST_DWARF'
        # rcat['DIST_TO'].name       = 'SEGUE_DIST_TO'
        # rcat['DIST_GIANT'].name    = 'SEGUE_DIST_GIANT'
        # rcat['DIST_AGB'].name      = 'SEGUE_DIST_AGB'
        # rcat['DIST_FHB'].name      = 'SEGUE_DIST_FHB'
        # rcat['DIST_AP'].name       = 'SEGUE_DIST_AP'
        # rcat['RV_FLAG'].name       = 'SEGUE_RV_FLAG'
        # rcat['RV_ADOP'].name       = 'SEGUE_RV'
        # rcat['RV_ADOP_UNC'].name   = 'SEGUE_RV_ERR'
        # rcat['mcat_matsep'].name   = 'SEGUE_H3_sep'

        for nc in self.newcols:
            if nc == 'nbands':
                rcat[nc] = np.zeros(len(rcat),dtype=int)
            elif nc == 'nspecpix':
                rcat[nc] = np.zeros(len(rcat),dtype=int)            
            else:
                rcat[nc] = np.nan * len(rcat)

        # make FLAG column
        rcat['FLAG'] = -1 * np.ones(len(rcat),dtype=int)

        print('making rcat...')
        # pass each set of PLATEID, FIBERID, MJD to search for pars file
        for rcat_i in tqdm(rcat):
            rcat_i = self.filltab(rcat_i)
        #list(map(self.filltab,rcat))
        #print(rcat[0])

        # translate columns names to be IDL readable
        for nc in rcat.keys():
            if '[Fe/H]' in nc:
                rcat[nc].name = nc.replace('[Fe/H]','FeH')
            if '[a/Fe]' in nc:
                rcat[nc].name = nc.replace('[a/Fe]','aFe')
            if 'log(g)' in nc:
                rcat[nc].name = nc.replace('log(g)','logg')
            if 'log(R)' in nc:
                rcat[nc].name = nc.replace('log(R)','logR')
            if 'log(L)' in nc:
                rcat[nc].name = nc.replace('log(L)','logL')
            if 'log(Age)' in nc:
                rcat[nc].name = nc.replace('log(Age)','logAge')
        for nc in rcat.keys():
            if 'initial_' in nc:
                rcat[nc].name = nc.replace('initial_','init_')

        # del rcat['RCHI2']
        # del rcat['DOF']
        # del rcat['SEGUE_G_MAG']
        # del rcat['SEGUE_DIST_DWARF']
        # del rcat['SEGUE_DIST_TO']
        # del rcat['SEGUE_DIST_GIANT']
        # del rcat['SEGUE_DIST_AGB']
        # del rcat['SEGUE_DIST_FHB']
        # del rcat['SEGUE_DIST_AP']
        # del rcat['MP_FLAG']

        # find duplicates and flag lower-SNR duplicates

        # make dup column
        rcat['dup'] = np.zeros(len(rcat),dtype=int)

        # find all dups
        # dup = np.unique(rcat['GAIAEDR3_ID'],return_inverse=True,return_counts=True)
        # dupGID = dup[0][dup[-1] > 1]

        # # for dups, figure out max SNR, set dup = 1 otherwise dup =2
        # for GID in dupGID:
        #     dupind = np.argwhere(rcat['GAIAEDR3_ID'] == GID).flatten()
        #     rcat['dup'][dupind] = 2.0

        #     maxsnr = np.argmax([rcat['SNR'][x] for x in dupind]).flatten()
        #     dupind_max = dupind[maxsnr]        
        #     rcat['dup'][dupind_max] = 1.0

        # for all dup = 2, set FLAG = 1
        # cond = rcat['dup'] > 1
        # rcat['FLAG'][cond] = 1
        
        rcat = self.calcSgr(rcat)
        
        # write rcat out
        if self.outfile is None:            
            outfile = "{OUTDIR}catalogs/{CATALOG}_rcat_{VERSION}_{RUNTYPE}.fits".format(
                OUTDIR=outdir,
                VERSION=self.version,
                RUNTYPE=self.runtype,
                CATALOG=self.catalog)
        else:
            outfile = self.outfile
        rcat.write(outfile,overwrite=True)

        print('success, rcat is made!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", "-a", help="which acat to use",type=str)
    parser.add_argument("--version", "-v", help="run version",type=str,default='V0.0')
    parser.add_argument("--runtype", "-r", help="runtype of results",type=str,default='MSG')
    parser.add_argument("--outfile", "-o", help="name of outfile",type=str,default=None)
    args = parser.parse_args()

    MK = mkrcat(args.catalog,args.version,args.runtype,args.outfile)
    MK.run()
