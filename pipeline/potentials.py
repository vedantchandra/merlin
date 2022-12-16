import numpy as np
import time
import os,sys,time,signal
import glob

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)
warnings.filterwarnings("ignore",category=UserWarning)

from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import astropy.units as u 

import gala.potential as gp
import gala.dynamics as gd
import gala.integrate as gi
from collections import OrderedDict


# stuff for timing action angle calculation
class TimeoutError(Exception):
    pass
def _sig_alarm(sig, tb):
    raise TimeoutError

class defaultMW(object):
    """docstring for DefaultMW"""
    def __init__(self,**kwargs):
        super(defaultMW, self).__init__()
        self.verbose = kwargs.get('verbose',False)
        self.pot=gp.MilkyWayPotential()

        self.grid = [(5*u.Myr, 5*u.Gyr), (100*u.Myr, 100*u.Gyr),]# (1*u.Myr, 25*u.Gyr),] 
                # (10*u.Myr, 250*u.Gyr), (25*u.Myr, 2500*u.Gyr), (10*u.Myr, 5000*u.Gyr)]

        # self.grid = [(100*u.Myr, 1*u.Gyr),]


    def compute_orbit(self,w0, dt=10*u.Myr, t1=0*u.Gyr, t2=250*u.Gyr, Integrator=gi.DOPRI853Integrator):
        with warnings.catch_warnings(record=True):
            try:
                orbit = self.pot.integrate_orbit(w0, dt=dt, t1=t1, t2=t2, Integrator=Integrator)
                iso = gd.fit_isochrone(orbit)
                orbit = orbit[orbit.energy(iso)<0]
                return orbit, iso
            except:
                return None, None

    #compare 80% of an orbit against the rest to see if actions have converged within 1%
    def check_convergence(self, orbit, iso, N_max=16):
        # with warnings.catch_warnings(record=True):
        if orbit is None:
            if self.verbose:
                print('Orbit is None object')
                sys.stdout.flush()
            return (
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]))

        try:
            if np.max(orbit.energy())<0:
                len_orbit = orbit.shape[0]

                res0 = (gd.find_actions(orbit[:int(-0.2*len_orbit)], N_max=N_max, toy_potential=iso))
                res1 = (gd.find_actions(orbit, N_max=N_max, toy_potential=iso))
                rel_err = np.abs((res1['actions'] - res0['actions'])/(res0['actions']))
                return rel_err, res1['actions'].to(u.kpc * u.km / u.s), res1['angles'], res1['freqs']

            else:
                if self.verbose:
                    print('max(E) > 0')
                    sys.stdout.flush()
                return (
                    np.array([np.nan, np.nan, np.nan]),
                    np.array([np.nan, np.nan, np.nan]), 
                    np.array([np.nan, np.nan, np.nan]),
                    np.array([np.nan, np.nan, np.nan]))
        except:
            if self.verbose:
                print('Checking Convergence Failed')
                print(orbit)
                print(orbit.energy())
                sys.stdout.flush()
            return (
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]))
                              
            
    #Compute actions across a grid of (timestep, total integration time) and stop as soon as actions converge within <1%
    #Starts from quickly computed gridpoints before moving on to finer integrations
    def iterative_actions(self, w0, runactions):
        
        grid_flag = 0
        break_flag = 0
        max_err_grid = []
        res_grid = []
        
        if runactions:
            for grid_pt in self.grid:
                orbit,iso = self.compute_orbit(w0, dt=grid_pt[0], t2=grid_pt[1])
                rel_err, actions, angles, freqs = self.check_convergence(orbit,iso)
                res_grid.append((rel_err, actions, angles, grid_flag, freqs, orbit))
                
                max_err_grid.append(max(rel_err[0], rel_err[2]))
                
                if max(rel_err[0], rel_err[2])<1e-2:
                    break_flag = 1
                    break
                grid_flag = grid_flag + 1
                
            if break_flag!=1: #When <1% isn't reached, report values for the gridpoint that produces the least deviation as per "check_convergence". These values aren't safe to use.
                rel_err, actions, angles, grid_flag, freqs, orbit = res_grid[np.argmin(max_err_grid)]
        else:
            if self.verbose:
                print('runactions turned off')
            rel_err, actions, angles, freqs = (
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]),
                np.array([np.nan, np.nan, np.nan]))
            orbit = None
            grid_flag = -1
            
        default_orbit,_ = self.compute_orbit(w0, dt=1*u.Myr, t2=10*u.Gyr)
            
        return rel_err, actions, angles, grid_flag, freqs, orbit, default_orbit

    def get_ecc(self, w0, dt=1*u.Myr, t2=100*u.Gyr):
        try:
            orbit = self.compute_orbit(w0, dt=dt, t2=t2)
            ecc = orbit.eccentricity()
            apo = orbit.apocenter()
        except:
            ecc = -1
            apo = -1*u.kpc
        return ecc, apo

    #Wrapper function to call potential computations on a size-1 astropy Table ("cat")
    def potentialize(self, w, runactions=True):
        rel_err, actions, angles, grid_flag, freqs, orbit, default_orbit = self.iterative_actions(w, runactions)
        return {'rel_err':rel_err,'actions':actions,'angles':angles,'freqs':freqs,'grid_flag':np.int16(grid_flag),'orbit':orbit, 'default_orbit':default_orbit}



class ToyPot(object):
    """docstring for DefaultMW"""
    def __init__(self,**kwargs):
        super(ToyPot, self).__init__()
        self.pot=gp.MilkyWayPotential()

    def compute_orbit(self,w0, dt=5*u.Myr, t1=5*u.Myr, t2=100*u.Gyr,Integrator=gi.DOPRI853Integrator):
        try:
            orbit = gp.Hamiltonian(self.pot).integrate_orbit(w0, dt=dt, t1=t1,t2=t2, Integrator=Integrator)
        except (IndexError,ValueError,TypeError,RuntimeError):
            return None, None

        try:
            toy_potential = gd.fit_isochrone(orbit)
            orbit = orbit[orbit.energy(toy_potential)<0]
            return orbit, toy_potential
        except (IndexError,ValueError,TypeError,RuntimeError):            
            return orbit, None

                             
    def iterative_actions(self, w0, N_max=8):
        orbit,toypot = self.compute_orbit(w0, dt=1.0*u.Myr, t1=0.0*u.Myr,t2=10*u.Gyr)
        # if toypot is None:
        #     # print("Could not fit orbit")
        return (
            np.array([np.nan, np.nan, np.nan]),
            np.array([np.nan, np.nan, np.nan]),
            np.array([np.nan, np.nan, np.nan]),
            np.array([np.nan, np.nan, np.nan]),
            orbit, orbit)

        # signal.signal(signal.SIGALRM, _sig_alarm)
        # try:
        #     signal.alarm(1)
        #     # try:
        #     #     # res0 = gd.find_actions(orbit[:int(-0.2*orbit.shape[0])], N_max=N_max, toy_potential=toypot)
        #     #     res1 = gd.find_actions(orbit, N_max=N_max, toy_potential=toypot)       
        #     #     # rel_err = np.abs((res1['actions'] - res0['actions'])/(res0['actions']))
        #     #     rel_err = np.array([0.0,0.0,0.0])
        #     #     return rel_err, res1['actions'].to(u.kpc * u.km / u.s), res1['angles'],res1['freqs'],orbit
        #     # except (IndexError,ValueError,TypeError,RuntimeError) as err:
        #     #     # print("Error: {0}".format(err))
        #     return (
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         orbit)
        # except TimeoutError:
        #     # print('Timed out')
        #     return (
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         np.array([np.nan, np.nan, np.nan]),
        #         orbit)

    #Wrapper function to call potential computations on a size-1 astropy Table ("cat")
    def potentialize(self, w):
        rel_err, actions, angles, freqs, orbit, default_orbit = self.iterative_actions(w)            
        return {'rel_err':rel_err,'actions':actions,'angles':angles,'freqs':freqs,'grid_flag':-1,'orbit':orbit, 'default_orbit':orbit}

