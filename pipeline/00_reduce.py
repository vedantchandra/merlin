from astropy.table import Table
import sys,argparse
import os
import glob

datadir = '/n/holyscratch01/conroy_lab/vchandra/mage/'
pipedir = '/n/home03/vchandra/outerhalo/08_mage/pipeline/'


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir',help='reduction directory',type=str,default=None)
    parser.add_argument('--version',help='Version of run',type=str,default='0')
    parser.add_argument('--dryrun',help='whether to only generate pypeit files without reducing',type=int,default=0)
    parser.add_argument('--skipred',help='whether to skip main reduction step',type=int,default=0)
    args = parser.parse_args()

    with open(pipedir + 'slurm/reduce_template.txt', 'r') as f:
        script = f.readlines()

    if args.dryrun == 0:
        dryrun_str = ''
    elif args.dryrun == 1:
        dryrun_str = '--dryrun'

    if args.skipred == 0:
        skipred_str = '--skipred=False'
    elif args.skipred == 1:
        skipred_str = '--skipred=True --restart=False'


    if args.dir == 'all':
        folders = glob.glob(datadir + 'data/202*')
        dirs = [os.path.basename(x) for x in folders]
    elif args.dir == 'bonaca':
        folders = glob.glob(datadir + 'data/b202*')
        dirs = [os.path.basename(x) for x in folders]
    else:
        dirs = [args.dir]

    print('reducing %i nights of data!' % len(dirs))

    for directory in dirs:

        print('making slurm script for %s' % directory)

        outscript = []
        for row in script:
            row = row.replace('{DIR}', directory)
            row = row.replace('{VERSION}', args.version)
            row = row.replace('{DRYRUN}', dryrun_str)
            row = row.replace('{SKIPRED}', skipred_str)
            outscript.append(row)

        outpath = pipedir + 'slurm/04_reduce_%s_%s.sh' % (directory, args.version)

        with open(outpath, 'w') as f:
            for line in outscript:
                f.write(line)

        print('submitting %s to slurm!' % directory)
        os.system('sbatch %s' % outpath)
