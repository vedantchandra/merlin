import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

import gala.potential as gp
import gala.dynamics as gd
from gala.units import DimensionlessUnitSystem

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)
warnings.filterwarnings("ignore",category=UserWarning)

import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from potentials import defaultMW, ToyPot

import sys,os,glob,shutil,gzip,subprocess,time,ast
from datetime import datetime
import numpy as np
from scipy.stats import truncnorm
from scipy.optimize import minimize
import h5py
import argparse
import functools,operator

class phaseafy(object):
     """docstring for phaseafy"""
     def __init__(self,):
          super(phaseafy, self).__init__()

          self.potentials = {'pot1':defaultMW(verbose=False)}
          self.galcen_distance = 8.122 * u.kpc
          self.galcen_v_sun = (12.9, 245.6, 7.78) * u.km / u.s
          self.z_sun = 20.8 * u.pc

          # circ setup from Rohan
          def pot(R):
               return gp.MilkyWayPotential().energy([R, 0, 0]*u.kpc).value[0]
          pot_vec = np.vectorize(pot)
          def Lcirc(Etot,R):
               return -R*((2*(Etot - pot_vec(R)))**0.5) 
          def maxLcirc(Etot):
               optfunc = functools.partial(Lcirc,Etot)
               res = minimize(optfunc, np.array([0.1]), method='BFGS')
               return np.abs(res.fun)
          maxLcirc_vec = np.vectorize(maxLcirc)
          self.maxLcirc_arr = maxLcirc_vec(np.linspace(-0.175, 0, 1000))
          

     def calcphase(self,outtab, nsamples=50000,verbose=False):
          if verbose:
               print('   ... Creating Phase Space Parameters')

          # nsamples == number of subsamples, if n = 1 then just use the value in outtab

          pararr = ['R_gal','X_gal','Y_gal','Z_gal',
               'Vx_gal','Vy_gal','Vz_gal',
               'Vr_gal','Vphi_gal','Vtheta_gal','V_tan','V_gsr',
               'Lx','Ly','Lz','Ltot']
          for POTN in self.potentials.keys():
               pararr.append('E_kin_{0}'.format(POTN))
               pararr.append('E_pot_{0}'.format(POTN))
               pararr.append('E_tot_{0}'.format(POTN))
               pararr.append('circLz_{0}'.format(POTN))
               pararr.append('circLtot_{0}'.format(POTN))

          for pp in pararr:
               outtab[pp] = np.nan
               outtab[pp+'_err'] = np.nan

          if verbose:
               print('   ... Using {} Samples'.format(nsamples))

          if nsamples == 1:
               RAarr    = [outtab['GAIAEDR3_RA']]
               Decarr   = [outtab['GAIAEDR3_DEC']]
               PMRAarr  = [outtab['GAIAEDR3_PMRA']]
               PMDecarr = [outtab['GAIAEDR3_PMDEC']]
               Distarr  = [outtab['Dist']]
               Vradarr  = [outtab['Vrad']]
          else:
               # build input maxtrix
               # mu = [RA,DEC,PMRA,PMDEC]
               # cov = [
               #    [V_ra,      V_dec_ra,   V_pmra_ra,    V_pmdec_ra  ],
               #    [V_ra_dec,  V_dec,      V_pmra_dec,   V_pmdec_dec ],
               #    [V_ra_pmra, V_dec_pmra, V_pmra,       V_pmdec_pmra],
               #    [V_ra_pmdec,V_dec_pmdec,V_pmra_pmdec, V_pmdec     ],
               #    ]

               mastodeg = 1.0 / 3600000.0

               V_ra         = (outtab['GAIAEDR3_RA_ERROR']*mastodeg)**2.0
               V_dec        = (outtab['GAIAEDR3_DEC_ERROR']*mastodeg)**2.0
               V_pmra       = outtab['GAIAEDR3_PMRA_ERROR']**2.0
               V_pmdec      = outtab['GAIAEDR3_PMDEC_ERROR']**2.0
               V_ra_dec     = (outtab['GAIAEDR3_RA_ERROR'] *mastodeg) * (outtab['GAIAEDR3_DEC_ERROR']*mastodeg) * outtab['GAIAEDR3_RA_DEC_CORR']
               V_ra_pmra    = (outtab['GAIAEDR3_RA_ERROR'] *mastodeg) * outtab['GAIAEDR3_PMRA_ERROR']  * outtab['GAIAEDR3_RA_PMRA_CORR']
               V_dec_pmra   = (outtab['GAIAEDR3_DEC_ERROR']*mastodeg) * outtab['GAIAEDR3_PMRA_ERROR']  * outtab['GAIAEDR3_DEC_PMRA_CORR']
               V_ra_pmdec   = (outtab['GAIAEDR3_RA_ERROR'] *mastodeg) * outtab['GAIAEDR3_PMDEC_ERROR'] * outtab['GAIAEDR3_RA_PMDEC_CORR']
               V_dec_pmdec  = (outtab['GAIAEDR3_DEC_ERROR']*mastodeg) * outtab['GAIAEDR3_PMDEC_ERROR'] * outtab['GAIAEDR3_DEC_PMDEC_CORR']
               V_pmra_pmdec = outtab['GAIAEDR3_PMRA_ERROR'] * outtab['GAIAEDR3_PMDEC_ERROR'] * outtab['GAIAEDR3_PMRA_PMDEC_CORR']

               mu = [outtab['GAIAEDR3_RA'],outtab['GAIAEDR3_DEC'],outtab['GAIAEDR3_PMRA'],outtab['GAIAEDR3_PMDEC']]
               cov = ([
                    [V_ra,      V_ra_dec,   V_ra_pmra,    V_ra_pmdec  ],
                    [V_ra_dec,  V_dec,      V_dec_pmra,   V_dec_pmdec ],
                    [V_ra_pmra, V_dec_pmra, V_pmra,       V_pmra_pmdec],
                    [V_ra_pmdec,V_dec_pmdec,V_pmra_pmdec, V_pmdec     ],
                    ])
               try:
                    astsamples = np.random.multivariate_normal(mu,cov,nsamples) 
               except:
                    # print(outtab['GAIADR2_PMRA'],outtab['GAIADR2_PMDEC'],outtab['GAIAEDR3_PMRA'],outtab['GAIAEDR3_PMDEC'])
                    # print('mu:',mu)
                    # print('cov',cov)
                    # raise
                    if verbose:
                         print('   ... Issue with ND norm')
                    return outtab

               RAarr    = astsamples[...,0]
               Decarr   = astsamples[...,1]
               PMRAarr  = astsamples[...,2]
               PMDecarr = astsamples[...,3]

               # RAarr    = np.random.normal(loc=outtab['GAIAEDR3_RA'],   scale=outtab['GAIAEDR3_RA_ERROR'],   size=nsamples)
               # Decarr   = np.random.normal(loc=outtab['GAIAEDR3_DEC'],  scale=outtab['GAIAEDR3_DEC_ERROR'],  size=nsamples)
               # PMRAarr  = np.random.normal(loc=outtab['GAIAEDR3_PMRA'], scale=outtab['GAIAEDR3_PMRA_ERROR'], size=nsamples)
               # PMDecarr = np.random.normal(loc=outtab['GAIAEDR3_PMDEC'],scale=outtab['GAIAEDR3_PMDEC_ERROR'],size=nsamples)

               Vradarr  = np.random.normal(loc=outtab['Vrad'],         scale=outtab['Vrad_err'],size=nsamples)

               distmean, diststd = outtab['Dist'],np.nanmean([outtab['Dist_uerr'],outtab['Dist_lerr']])

               if np.isfinite(distmean) and np.isfinite(diststd):
                    pass
               else:
                    # print(distmean,diststd,outtab['dist_adpt'],outtab['dist_adpt_uerr'],outtab['dist_adpt_lerr'],outtab['dist_adpt_err'])
                    if verbose:
                         print('   ... Issue with distance mean/std')
                    return outtab
               try:
                    a, b = (0.0 - distmean) / diststd, (200.0 - distmean) / diststd
                    Distarr = truncnorm.rvs(a,b,loc=distmean,scale=diststd,size=nsamples)
               except:
                    if verbose:
                         print('   ... Issue with drawing distances')
                    return outtab


          # now calculate parameters using gala
          dist = Distarr*u.kpc
          vrad = Vradarr*u.km/u.s
          ra = RAarr*u.deg
          dec = Decarr*u.deg
          pmra = PMRAarr*u.mas/u.yr
          pmdec = PMDecarr*u.mas/u.yr

          ceq = coord.ICRS(
               ra=ra, dec=dec, distance=dist, 
               pm_ra_cosdec=pmra, pm_dec=pmdec,
               radial_velocity=vrad)

          cgal = ceq.transform_to(coord.Galactocentric(
               galcen_distance=self.galcen_distance,
               galcen_v_sun=self.galcen_v_sun,
               z_sun=self.z_sun,
               ))
          w0 = gd.PhaseSpacePosition(cgal.cartesian)

          x = np.array([w0.pos.x.value, w0.pos.y.value, w0.pos.z.value]) * w0.pos.x.unit
          v = np.array([w0.vel.d_x.value, w0.vel.d_y.value, w0.vel.d_z.value]) * w0.vel.d_x.unit
          L = np.cross(x.value, v.value, axis=0) * w0.pos.x.unit * w0.vel.d_x.unit
          Ltot = np.linalg.norm(L.value, axis=0) * L.unit
          Lx = L[0]
          Ly = L[1]
          Lz = L[2]

          energydict = {}
          for POTN in self.potentials.keys():
               # Ek = 0.5*(np.linalg.norm(v.value, axis=0) * v.unit)**2
               Ek   = w0.kinetic_energy().to(u.km**2*u.s**-2)
               Epot = w0.potential_energy(self.potentials[POTN].pot).to(u.km**2*u.s**-2)
               Etot = w0.energy(self.potentials[POTN].pot).to(u.km**2*u.s**-2)
               energydict[POTN] = {'Ek':Ek,'Epot':Epot,'Etot':Etot}

          Vtan = 4.74 * np.sqrt( (PMRAarr**2.0) + (PMDecarr**2.0)) * dist

          w0 = gd.PhaseSpacePosition(cgal.sphericalcoslat)

          # calculate RV_GSR
          ceq1 = coord.ICRS(
               ra=ra, dec=dec,
               radial_velocity=vrad)

          # v_sun = coord.galactocentric_frame_defaults.get()['galcen_v_sun'].to_cartesian()
          v_sun = cgal.galcen_v_sun.to_cartesian()

          cgal1 = ceq1.transform_to(coord.Galactic)
          cart_data = cgal1.data.to_cartesian()
          unit_vector = cart_data / cart_data.norm()
          v_proj = v_sun.dot(unit_vector)

          V_gsr = ceq1.radial_velocity + v_proj

          R_gal      = np.sqrt((x[0]**2.0)+(x[1]**2.0)+(x[2]**2.0)).value
          X_gal      = x[0].value
          Y_gal      = x[1].value
          Z_gal      = x[2].value
          Vx_gal     = v[0].value
          Vy_gal     = v[1].value
          Vz_gal     = v[2].value
          V_tan      = Vtan.value
          V_gsr      = V_gsr.value
          Vr_gal     = w0.radial_velocity.value
          Vtheta_gal = (w0.distance * w0.pm_lat).to(u.km/u.s, u.dimensionless_angles()).value
          Vphi_gal   = (w0.distance * w0.pm_lon_coslat).to(u.km/u.s, u.dimensionless_angles()).value
          Lx         = Lx.value
          Ly         = Ly.value
          Lz         = Lz.value
          Ltot       = Ltot.value

          outtab['R_gal']      = np.nanmean(R_gal     )
          outtab['X_gal']      = np.nanmean(X_gal     )
          outtab['Y_gal']      = np.nanmean(Y_gal     )
          outtab['Z_gal']      = np.nanmean(Z_gal     )
          outtab['Vx_gal']     = np.nanmean(Vx_gal    )
          outtab['Vy_gal']     = np.nanmean(Vy_gal    )
          outtab['Vz_gal']     = np.nanmean(Vz_gal    )
          outtab['V_tan']      = np.nanmean(V_tan     )
          outtab['V_gsr']      = np.nanmean(V_gsr     )
          outtab['Vr_gal']     = np.nanmean(Vr_gal    )
          outtab['Vtheta_gal'] = np.nanmean(Vtheta_gal)
          outtab['Vphi_gal']   = np.nanmean(Vphi_gal  )
          outtab['Lx']         = np.nanmean(Lx        )
          outtab['Ly']         = np.nanmean(Ly        )
          outtab['Lz']         = np.nanmean(Lz        )
          outtab['Ltot']       = np.nanmean(Ltot      )

          for POTN in self.potentials.keys():
               outtab['E_kin_{0}'.format(POTN)] = np.nanmean(energydict[POTN]['Ek'].value )
               outtab['E_pot_{0}'.format(POTN)] = np.nanmean(energydict[POTN]['Epot'].value )
               outtab['E_tot_{0}'.format(POTN)] = np.nanmean(energydict[POTN]['Etot'].value )

          outtab['R_gal_err']      = np.nanstd(R_gal     )
          outtab['X_gal_err']      = np.nanstd(X_gal     )
          outtab['Y_gal_err']      = np.nanstd(Y_gal     )
          outtab['Z_gal_err']      = np.nanstd(Z_gal     )
          outtab['Vx_gal_err']     = np.nanstd(Vx_gal    )
          outtab['Vy_gal_err']     = np.nanstd(Vy_gal    )
          outtab['Vz_gal_err']     = np.nanstd(Vz_gal    )
          outtab['V_tan_err']      = np.nanstd(V_tan     )
          outtab['V_gsr_err']      = np.nanstd(V_gsr     )
          outtab['Vr_gal_err']     = np.nanstd(Vr_gal    )
          outtab['Vtheta_gal_err'] = np.nanstd(Vtheta_gal)
          outtab['Vphi_gal_err']   = np.nanstd(Vphi_gal  )
          outtab['Lx_err']         = np.nanstd(Lx        )
          outtab['Ly_err']         = np.nanstd(Ly        )
          outtab['Lz_err']         = np.nanstd(Lz        )
          outtab['Ltot_err']       = np.nanstd(Ltot      )
          for POTN in self.potentials.keys():
               outtab['E_kin_{0}_err'.format(POTN)] = np.nanstd(energydict[POTN]['Ek'].value )
               outtab['E_pot_{0}_err'.format(POTN)] = np.nanstd(energydict[POTN]['Epot'].value )
               outtab['E_tot_{0}_err'.format(POTN)] = np.nanstd(energydict[POTN]['Etot'].value )

          # do circ calc
          for POTN in self.potentials.keys():
               Lmax = np.interp(energydict[POTN]['Etot'].value/(1E+6),np.linspace(-0.175,0,1000),self.maxLcirc_arr)
               outtab['circLz_{0}'.format(POTN)]   = np.nanmean(np.abs(Lz / Lmax)/1000.0)
               outtab['circLtot_{0}'.format(POTN)] = np.nanmean(np.abs(Ltot / Lmax)/1000.0)
               outtab['circLz_{0}_err'.format(POTN)]   = np.nanstd(np.abs(Lz / Lmax)/1000.0)
               outtab['circLtot_{0}_err'.format(POTN)] = np.nanstd(np.abs(Ltot / Lmax)/1000.0)

          return outtab
