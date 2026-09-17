from astropy.table import Table
import sys,argparse
import os
import glob
import numpy as np

datadir = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/'
pipedir = '/n/home03/vchandra/outerhalo/08_mage/pipeline/'


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--catalog',help='catalog to use as input',type=str,default=None)
    parser.add_argument('--version',help='Version of run',type=str,default='V0.0')
    parser.add_argument('--npoints',help='number of points',type=int,default=500)
    parser.add_argument('--skipfit',help='skip MS fit',type=int,default=0)
    parser.add_argument('--overwrite',help='overwrite already-fitted stars',
                                type=int,default=0)
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--no-submit', dest='submit', action='store_false')
    parser.add_argument('--arr1', help='start array idx',type=int,default=0)
    parser.add_argument('--arr2', help='end array idx',type=int,default=None)
    parser.add_argument('--sel', help = 'fit single selection', type = str, default = None)
    parser.set_defaults(submit=True)
    args = parser.parse_args()

    catalog = args.catalog
    version = args.version 
    npoints = args.npoints
    skipfit = args.skipfit
    submit = args.submit
    arr1 = args.arr1 
    arr2 = args.arr2
    selection = args.sel
    overwrite = bool(int(args.overwrite))

    if overwrite:
        print('overwriting past fits!')

    if catalog is None:
        print('please pass a valid catalog to use!')
        raise

    acatpath = datadir + 'catalogs/' + catalog + '_acat.fits'
    print('reading %s...' % acatpath)

    tab = Table.read(acatpath)
    nstar = len(tab)

    if arr2 is None:
        arr2 = nstar - 1

    try:
        os.makedirs(datadir + 'logs/' + catalog + '/' + version + '/')
    except:
        print('logs dir exists for this catalog...')

    if overwrite:
        print('submitting %i stars to MINESWEEPER...' % (arr2 - arr1 + 1))
    else:
        print('skipping stars from version %s that have already been fitted...' % version)
        sampledir = datadir + 'pars/%s/%s/' % (catalog, version)
        samplefiles = glob.glob(sampledir + '*.pars')
        print('there are %i sample files' % len(samplefiles))
        # a fit is identified by (Gaia ID, observation date): the same star observed on
        # two nights needs two fits. Filenames are mage_<GaiaID>_<date>_<VER>.pars, where
        # <date> itself contains underscores (e.g. b2026_02_07).
        done = set()
        for file in samplefiles:
            stem = os.path.basename(file)[len('mage_'):-len('_%s.pars' % version)]
            gid, date = stem.split('_', 1)
            done.add((gid, date))

        print('%i (GaiaID, date) pairs already fitted' % len(done))

        isdone = np.array([(str(g), str(d)) in done for g, d in zip(tab['GAIAEDR3_ID'], tab['date'])])
        tofit_idx = tab['ACAT_ID'][~isdone]
        tofit_idx = [str(x) for x in tofit_idx]
        print('ACAT_IDs to fit:', tofit_idx)

        print('there are %i stars that need to be fit' % len(tofit_idx))

        tofit_idx_list = ','.join(tofit_idx)

    with open(pipedir + 'slurm/runms_template.txt', 'r') as f:
        script = f.readlines()

    outscript = []
    for row in script:
        row = row.replace('{CATALOG}', catalog)
        row = row.replace('{VERSION}', version)
        row = row.replace('{ARRAY}', '0-%i' % (nstar - 1))
        row = row.replace('{NPOINTS}', '%i' % (npoints))
        row = row.replace('{SKIPFIT}', '%i' % skipfit)
        
        if overwrite:
            row = row.replace('{ARR1}-{ARR2}', '%i-%i' % (arr1, arr2))
        else:
            row = row.replace('{ARR1}-{ARR2}', '%s' % tofit_idx_list)

        outscript.append(row)

    outpath = pipedir + 'slurm/04_runms_%s_%s.sh' % (catalog, version)

    with open(outpath, 'w') as f:
        for line in outscript:
            f.write(line)

    if submit:
        os.system('sbatch %s' % outpath)
