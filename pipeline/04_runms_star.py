from astropy.table import Table
import sys,argparse
import os



datadir = '/n/holyscratch01/conroy_lab/vchandra/mage/'
outdir = datadir

def run(catalog=None,ind=None,version='V0.0', npoints = 250, skipfit = 0):

    if skipfit == 0:
        import runMS


    mattab = Table.read(datadir + 'catalogs/' + catalog + '_acat.fits', format='fits')
    mattab_i = mattab[ind]
    gaiaID  = mattab_i['GAIAEDR3_ID']
    #field=mattab_i['FIELD']
    mjd=mattab_i['date']
    acat_id = mattab_i['ACAT_ID']

    try:
        os.makedirs(datadir + 'plots/' + catalog + '/' + version + '/')
    except:
        print('plot dir exists for this catalog...')

    try:
        os.makedirs(datadir + 'samples/' + catalog + '/' + version + '/')
    except:
        #raise
        print('samples dir exists for this catalog...')

    try:
        os.makedirs(datadir + 'pars/' + catalog + '/' + version + '/')
    except:
        print('pars dir exists for this catalog...')



    samplefile = '{OUTDIR}{CATALOG}/{VER}/mage_{GAIAID}_{MJD}_{VER}_samp.dat'.format(
            GAIAID=gaiaID,
            MJD=mjd,
            VER=version,
            OUTDIR=outdir,
            CATALOG=catalog)

    print('sample file is %s' %samplefile)

    if os.path.exists(samplefile) and skipfit == 0:
        print('sample file exists, deleting...')
        os.remove(samplefile)

    # sys.stdout.flush()
    # runMS.run(GaiaID=gaiaID,version=version,catalog = catalog, npoints = npoints,
    #     field = field, mjd = mjd)
    # sys.stdout.flush()
    # compmod.run(GaiaID=gaiaID,version=version,catalog = catalog,
    #     field = field, mjd = mjd)
    # sys.stdout.flush()
    # corner.run(GaiaID=gaiaID,version=version,catalog = catalog,
    #     field = field, mjd = mjd)
    # sys.stdout.flush()
    # calcpars.run(GaiaID=gaiaID,version=version,catalog = catalog,
    #     field = field, mjd = mjd)
    # sys.stdout.flush()

    if skipfit == 0:
        print('running MS...')
        sys.stdout.flush()
        runMS.run(acat_id = acat_id,version=version,catalog = catalog, npoints = npoints)
    elif skipfit == 1:
        print('skipping MS fit...')
    else:
        print('invalid skipfit value! must be 0 or 1')
        raise
    sys.stdout.flush()
    print('running calcpars...')

    import calcpars
    calcpars.run(acat_id = acat_id,version=version,catalog = catalog)
    sys.stdout.flush()
    print('running compmod...')
    import compmod
    compmod.run(acat_id = acat_id,version=version,catalog = catalog)
    sys.stdout.flush()
    print('running corner...')
    import corner
    corner.run(acat_id = acat_id,version=version,catalog = catalog)
    sys.stdout.flush()
    print('finished 04_runms_star!')



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    # parser.add_argument('--survey',help='SEGUE Survey ID',type=str,
    #     choices=['SEGUE','SEGUE_clusters'],default='SEGUE')
    parser.add_argument('--catalog',help='catalog to use as input',type=str,default=None)
    parser.add_argument('--ind',help='Index of star',type=int,default=None)
    # parser.add_argument('--batch',help='number of 10k batches',type=int,default=0)
    parser.add_argument('--version',help='Version of run',type=str,default='V0.0')
    parser.add_argument('--npoints', help = "number of points", type = int, default = 250)
    parser.add_argument('--skipfit', help = 'skip MS fit', type = int, default = 0)
    args = parser.parse_args()

    run(catalog=args.catalog,ind=args.ind,version=args.version,npoints=args.npoints,
        skipfit = args.skipfit)
