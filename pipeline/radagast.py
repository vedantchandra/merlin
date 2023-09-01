###############################################################################################
# RADAGAST
# an end-to-end pipeline to reduce data from the Magellan Echellete spectrograph
# on the 6.5m Baade telescope at Las Campanas Observatory
# written by Vedant Chandra (vedant.chandra@cfa.harvard.edu)

# example usage:
# python radagast.py --dir=/Users/vedantchandra/magedata/night1/ --version=1 --restart=True --skipred=False
# add --dryrun flag to generate obslog and pypeit file for inspection without running full script

# CHANGELOG:
# 11/09/2022: initial version
# 12/12/2022: disable spatial flexure, increase RMS threshold to 0.4 for wavecal
# 			  change input file handling to match pypeit v1.11.1 syntax
#			  skipred now automatically disables restart, to keep master files
# 12/12/2022: change xe-flash to trace only
#			  turn off global_sky_std 
# 12/20/2022: add hip108327 to standards. raise error now if flux standard missing from lib
#             set edge thresh to 3, now only use xe-flash for tracing (not red flat)
#             set findobj snr_thresh 10 -> 3

###############################################################################################

# make a library of your flux standards here
# the code will find the relevant flux standard within the directory

flux_standards = ['hip77', 'hip67523', 'hip104326', 'hip108327', 'hip17946', 'hip51633',
		  		 'hip18271'] # these are for Ana's program. none are A stars. 

###############################################################################################
# IMPORTS
###############################################################################################

import matplotlib
matplotlib.use('agg')
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from scipy.signal import medfilt
from numpy.polynomial.polynomial import polyval,polyfit
import sys
import shutil
import os
from astropy.io import ascii
import scipy
import argparse
import time
import gc
start = time.time()

#plt.style.use('vedant')
class bcolors:
	HEADER =  '\033[91m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

###############################################################################################
# TAKE COMMAND-LINE INPUTS
###############################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--dir', '-d', help='reduction directory with a raw/ subdirectory', type=str, default='./')
parser.add_argument('--version','-v', help='version of reduction', type=str, default='0')
parser.add_argument('--restart', '-r', help='whether to clear previous cals', type=str, default='True')
parser.add_argument('--skipred', '-s', help='skip main reduction, only flux and make QA plots. if so, set --restart=False too!', type=str, default='False')
parser.add_argument('--dryrun', dest='dryrun', action='store_true')
# parser.add_argument('--no-feature', dest='feature', action='store_false')
parser.set_defaults(dryrun=False)
args = parser.parse_args()

###############################################################################################
# CONTROL PANEL
###############################################################################################

print(args.dir)

workdir = args.dir
version = args.version
dryrun = args.dryrun

if workdir[-1] != '/':
	print('adding / to work directory!!')
	workdir += '/'

print('                      ')
print(bcolors.HEADER + '##################################' + bcolors.ENDC)
print(bcolors.HEADER + '             RADAGAST             ' + bcolors.ENDC)
print(bcolors.HEADER + '          MAGE REDUCTION          ' + bcolors.ENDC)
print(bcolors.HEADER + '  vedant.chandra@cfa.harvard.edu  ' + bcolors.ENDC)
print(bcolors.HEADER + '##################################' + bcolors.ENDC)
print('                      ')
print(bcolors.HEADER + 'work directory: %s' % workdir + bcolors.ENDC)
print(bcolors.HEADER + 'version: %s' % version + bcolors.ENDC)
print('                      ')

if args.restart.lower() == 'true':
	restart = True
elif args.restart.lower() == 'false':
	restart = False
else:
	print('invalid restart boolean! deleting previous cals...')
	restart = True


if args.skipred.lower() == 'true':
	skipred = True
	restart = False
	print('skipping main reduction, keeping old master files...')
elif args.skipred.lower() == 'false':
	skipred = False
	print('performing full main reduction step...')
else:
	print('invalid skip reduction boolean! redoing reduction...')
	skipred = True


gunzip = True # set to True if any raw files need to be gunzipped first
spectral_flexure = 'skip' # can be 'skip' or 'boxcar'
spatial_flexure = 'False' # can be true or false

###############################################################################################
# MAKE DIRECTORIES 
###############################################################################################

rawdir = workdir + 'raw/' # assume this exists with raw files in it
reddir = workdir + 'reduced_v%s/' % version

try:
        print('making reduction directory: %s' % reddir)
        os.mkdir(reddir)
except Exception as e:
        
	print('reduction directory already exists...')
	print(e)
        #raise
	
	if restart:
		print('deleting old reduction directory...')
		shutil.rmtree(reddir)
		os.mkdir(reddir)
		
print('there are %i raw frames' % len(glob.glob(rawdir + '*')))

if gunzip:
	os.system('gunzip %s*.fits.gz' % rawdir)

os.chdir(reddir)

print('now working in this directory: %s' % os.getcwd())

###############################################################################################
# PAIR ARCS TO SCIENCE EXPOSURES 
###############################################################################################

print('pairing arcs to science exposures...')
os.system('pypeit_obslog -r %s -d %s -o -f obslog.txt magellan_mage' % (rawdir, rawdir))

log = ascii.read(rawdir + 'obslog.txt', format = 'fixed_width')
logcol = list(log.columns)

log['target'] = [target.lower() for target in log['target']]

print('targets are: ')
print(list(log['target']))

# ANY SPECIAL-CASE FILES THAT NEED PROCESSING

for row in log:
	if row['filename'] == 'mage1067.fits' and row['target'] == 'None':
		row['target'] = 'thar' # rename bad filename from nov3

log.remove_column('calib')

sci = np.array(['j' in name for name in log['target']])
std = np.array(['hip' in name.lower() for name in log['target']])
arc = np.array(['thar' in name.lower() for name in log['target']])

print(std)

n_cal = np.sum(sci) + np.sum(std)
dummy_cal = 0

log['arcfile'] = '                 '
log['calib'] = '                 '

calib_idx = 1

for row in log:
	
	
	# find nearest arc for each science frame
	
	if 'j' in row['target'] or 'hip' in row['target']:
		arc_dists = log['mjd'] - row['mjd'] #np.sqrt(((row['ra'] - log['ra']) * np.cos(np.radians(log['dec'])))**2 + (row['dec'] - log['dec'])**2)
		arc_dists[arc_dists < 0] = 999
		nearest_arc = np.argmin(arc_dists[arc])
		row['calib'] = calib_idx
		row['arcfile'] = log[arc][nearest_arc]['filename']
		
		calib_idx += 1
		
	else: 
		row['calib'] = 'all'
		
		
	## reset file types
	
	if 'flat' in row['target'].lower():
		row['frametype'] = 'illumflat,pixelflat'
	elif 'flash' in row['target'].lower():
		row['frametype'] = 'trace'
	elif 'j' in row['target'].lower():
		row['frametype'] = 'science'
	elif 'hip' in row['target'].lower():
		row['frametype'] = 'science'
	elif 'thar' in row['target'].lower():
		row['frametype'] = 'arc,tilt'
		
	row['comb_id'] = -1

# combine science frames with the same target name

uniq_targets = list(np.unique(log['target']))
print('unique targets are:')
print(uniq_targets)

coadd_ctr = 1
for ctr,targ in enumerate(uniq_targets):

	if 'j' in targ or 'hip' in targ:
		sel = log['target'] == targ

		log['comb_id'][sel] = coadd_ctr
		log['calib'][sel] = log['calib'][sel][0]

		coadd_ctr += 1

# Set calib strings


log['arcfile_str'] = [file.strip() for file in log['arcfile']]
for row in log:
	if 'thar' in row['target'].lower():
		scisel = log['arcfile_str'] == row['filename'].strip()
		
		if np.sum(scisel) == 0:
			row['calib'] = str(dummy_cal)
			continue
		calibs = np.unique(list(log[scisel]['calib'])) # np unique here prevents dupes
		cal_str = ",".join(str(x) for x in calibs)
		
		row['calib'] = cal_str

# find which star is the flux standard

telluric = None
for name in log['target']:
	if name.lower() in flux_standards:
		telluric = name
		break # use the first telluric in the list
		
if telluric is None:
	print('there is no flux standard from the library in this folder!')
	raise
else:
	print('using %s as the flux standard...' % telluric)


# prevent 'all' in error by setting explicit ints for all calibs. nope this is not a problem

# calib_ints = list(log['calib'])
# calib_ints = ','.join(calib_ints)
# calib_ints = calib_ints.replace('all,', '')
# print(calib_ints)

# calib_ints = np.unique(np.array(calib_ints.split(',')))
# print(calib_ints)

# allcalib_str = ','.join(calib_ints)

# for row in log:
# 	if row['calib'] == 'all':
# 		row['calib'] = allcalib_str

# write log

ascii.write(log[logcol], format = 'fixed_width', output = rawdir + 'obslog_edited.txt', overwrite = True)

###############################################################################################
# GENERATE PYPEIT FILE
###############################################################################################

print('generating pypeit setup file...')

pypfile = reddir + 'magellan_mage_A/magellan_mage_A.pypeit'

try:
	print('deleting pypeit file if it already exists...')
	os.remove(pypfile)
except:
	print('no pypeit file exists, making a new one...')

os.system('pypeit_setup -b -r %s -s magellan_mage -c A' % rawdir)

os.chdir(reddir + 'magellan_mage_A/')

pypfile = reddir + 'magellan_mage_A/magellan_mage_A.pypeit'

outdir = reddir + 'magellan_mage_A/'

# Make edits to pypeit file

with open(pypfile) as f:
	lines = f.readlines()

newlines = [];

for line in lines:
	
	newlines.append(line)
	
	if 'spectrograph' in line:
					
		newlines.append('[scienceframe]\n')
		newlines.append('  [[process]]\n')
		newlines.append('    spat_flexure_correct=%s\n' % spatial_flexure)
		
		
		newlines.append('[flexure]\n')
		newlines.append('  spec_method = %s\n' % spectral_flexure) #FLEXURE CORRECTION
			
		newlines.append('[baseprocess]\n')
		newlines.append('  use_biasimage=False\n')


		newlines.append('[calibrations]\n')
		newlines.append('  [[wavelengths]]\n')
		newlines.append('    rms_threshold=0.4\n')
		newlines.append('  [[slitedges]]\n')
		newlines.append('    edge_thresh=1\n')

		newlines.append('[reduce]\n')
		newlines.append('  [[findobj]]\n')
		newlines.append('    maxnumber_sci=1\n')
		newlines.append('    snr_thresh=3\n')
		newlines.append('  [[skysub]]\n')
		newlines.append('    global_sky_std=False\n')
		
	if 'path' in line:
		break

os.remove(pypfile)

# Write new pypeit file

with open(pypfile, 'w') as f:
	for line in newlines:
		f.write("%s" % line)


with open(rawdir + 'obslog_edited.txt') as f:
	obslog = f.readlines()[3:]

with open(pypfile, 'a') as f:
	for line in obslog:
		f.write(line)
		
	f.write('data end')


# Check pypeit file

print('pypeit file has been edited and generated:')

with open(pypfile, 'r') as f:
	lines = f.readlines()
	
for line in lines:
	print(line)


if dryrun:
	print(bcolors.HEADER + 'dry run ended, terminating script!' + bcolors.ENDC)
	sys.exit()

###############################################################################################
# RUN MAIN REDUCTION
###############################################################################################

if skipred:
	print(bcolors.HEADER + 'skipping main reduction...' + bcolors.ENDC)
else:
	print(bcolors.HEADER + 'running main reduction! buckle up and grab some coffee...' + bcolors.ENDC)
	os.system('run_pypeit %s -o' % pypfile)
	print(bcolors.HEADER + 'main reduction completed!' + bcolors.ENDC)

gc.collect()

###############################################################################################
# MAKE PREVIEW PLOTS
###############################################################################################

outfiles = glob.glob(outdir + 'Science/spec1d*.fits')
print('there are %i reduced science exposures...' % len(outfiles))
try:
	os.mkdir(outdir + 'plots/')
except:
	print('output plot directory exists, deleting it and remaking...')
	shutil.rmtree(outdir + 'plots/')
	os.mkdir(outdir + 'plots/')


print('making order preview plots...')

for outfile in outfiles:

	name = outfile.split('_')[-3].split('-')[1]

	f,axs = plt.subplots(4, 3, figsize = (35, 10), sharey = False)

	with fits.open(outfile) as f:

		nspec = f[0].header['NSPEC']

		if nspec < 12:
			print('there are only %i/12 orders for %s!!!' % (nspec, name))

		order_list = np.arange(nspec)

		for order in order_list:

			plt.sca(axs.ravel()[order])

			try:
				wl, fl, ivar = f[order+1].data['OPT_WAVE'], f[order+1].data['OPT_COUNTS'], f[order+1].data['OPT_COUNTS_IVAR']
			except Exception as e:
				print('order preview failed for %s order %i' % (name, order))
				print(e)
				continue
			sig = 1 / np.sqrt(ivar)
			snr = np.nanmedian(fl * np.sqrt(ivar))

			#cont = scipy.signal.medfilt(fl, 251)
			plt.plot(wl, fl, color = 'k')

			cl = np.isfinite(fl)

			y1 = np.nanquantile(fl, 0.01)
			y2 = np.nanquantile(fl, 0.95)
			plt.ylim(0.01 * y1, 1.5 * y2)

			plt.text(0.95, 0.9, 'S/N = %.1f' % snr, ha = 'right', va = 'top', bbox = dict(color = 'w', boxstyle = 'round'), transform = plt.gca().transAxes)
			
	plt.suptitle(name, y = 0.95)

	plt.tight_layout()

	plt.savefig(outdir + 'plots/preview_%s.png' % name, dpi = 200)
	plt.close()


###############################################################################################
# FLUX CALIBRATION
###############################################################################################

print(bcolors.HEADER + 'starting flux calibration...' + bcolors.ENDC)

scidir = outdir + 'Science/'

os.chdir(scidir)

mksens_list = ['[sensfunc]', # ASSUMES an A0 type telluric standard. CHANGE THE STAR_MAG according to your star
'    algorithm = UVIS',
'    star_mag = 7.27',
'    star_type = A0',
'    polyorder = 5',
'    [[UVIS]]',
'        telluric=True',
'        nresln=20',
'        telluric_correct=True',
'        trans_thresh=0.9',
'        extinct_correct = False']

with open(scidir + 'make_sens.sens', 'w') as f:
	for line in mksens_list:
		f.write('%s\n' % line)


tellfiles = glob.glob(scidir + 'spec1d*%s*.fits' % telluric)
print('there are %i telluric files, picking the first one...' % len(tellfiles))
tellfile = tellfiles[0] # assume you only use 1 telluric exposure
print('using %s as the telluric' % tellfile)

print('making sensitivity function...')

status = os.system('pypeit_sensfunc %s -o sensfunc.fits -s make_sens.sens' % tellfile) # if you want debug plots, add --debug here

print(status)

if status != 0:
	print('error running pypeit_sensfunc!!')
	raise

scifiles = glob.glob(scidir + 'spec1d*.fits')

fluxfile_list = [];

fluxfile_list.append('flux read')
fluxfile_list.append('filename | sensfile')
for ii,scifile in enumerate(scifiles):
	line = scifile
	fluxfile_list.append(line + ' | sensfunc.fits')
fluxfile_list.append('flux end')

with open(scidir + 'fluxfile.txt', 'w') as f:
	for line in fluxfile_list:
		f.write('%s\n' % line)


print(bcolors.HEADER + 'flux calibrating science exposures...' + bcolors.ENDC)

os.system('pypeit_flux_calib fluxfile.txt')

targets = [];

for scifile in scifiles:

	f = fits.open(scifile)

	target = f[0].header['TARGET']

	f.close()

	targets.append(target)
	
targets = list(np.unique(targets))

fluxed_tellfile = glob.glob(scidir + 'spec1d*%s*.fits' % telluric)[0]

print('making fluxed standard star plot...')

plt.figure(figsize = (15, 5))

with fits.open(fluxed_tellfile) as f:
	for order in range(len(f) - 2):
		wl, fl = f[order+1].data['OPT_WAVE'], f[order+1].data['OPT_FLAM']
		
		if order == 0:
			medfl = np.median(fl)
			
		else:
			medfl = np.median([medfl, np.median(fl)])
		plt.plot(wl, fl, label = f[order + 1].header['NAME'])
		#print(f[order].header)
		
plt.title('fluxed standard: %s' % telluric)
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('flux')
plt.ylim(0.1 * medfl, 15 * medfl)
plt.legend(fontsize = 12)
plt.savefig(outdir + 'plots/fluxed_standard_%s.png' % telluric)


print(bcolors.HEADER + 'collating and stitching science spectra...' + bcolors.ENDC)

coadd_dir = scidir + 'coadd/'

try:
	os.mkdir(coadd_dir)
except:
	print("dir exists!")
	
try:
	os.mkdir(scidir + 'red_plots/')
except:
	print("dir exists!")


for target in targets:
	
	print('stiching %s...' % target)

	targetfiles = glob.glob(scidir + 'spec1d*%s*.fits' % target)

	coadd_list = [];

	coadd_list.append('[coadd1d]')
	coadd_list.append('  coaddfile=coadd/%s_coadd.fits' % target)
	coadd_list.append('  sensfuncfile = \'sensfunc.fits\'')
	coadd_list.append('  wave_method = velocity')
	#coadd_list.append('  spec_samp_fact = 1')


	coadd_list.append('  coadd1d read')
	coadd_list.append('    filename | obj_id')

	for scifile in targetfiles:

		txtfile = scifile[:-5] + '.txt'

		tab = ascii.read(txtfile, names = ['adsf', 'order', 'name', 'spat', 'frac', 'box', 'fwhm', 's2n', 'wv', 'blah'])

		for obj in tab[1:]:
			line = '    ' + scifile + ' | ' + obj['name'].strip()   
			coadd_list.append(line)
			break

		#print(line)

	coadd_list.append('  coadd1d end')

	coadd_list

	with open(scidir + 'make_coadd_%s.txt' % target, 'w') as f:
		for line in coadd_list:
			f.write('%s\n' % line)

	os.system('pypeit_coadd_1dspec make_coadd_%s.txt' % (target))


# Make coadd preview plots

print('making coadd preview plots...')

for target in targets:
	coaddfile = scidir + 'coadd/' + target + '_coadd.fits'

	try:
		f = fits.open(coaddfile)
	except:
		print('failed to open coadd file for %s' % coaddfile)
		continue
	
	f[0].header

	plt.figure(figsize = (75, 5))
	plt.subplot(121)
	plt.plot(f[1].data['wave'], f[1].data['flux'], color = 'k', lw = 0.5)
	plt.ylabel('flux')
	plt.xlim(3500, 9000)
	
	medfl = np.nanmedian(f[1].data['flux'])
	
	if target == telluric:
		height = 5
	else:
		height = 2
	
	plt.ylim(- 0.4 * medfl, height * medfl)
	plt.xlabel('wavelength [$\mathrm{\AA}$]')
	
	
	snr = np.nanmedian(f[1].data['flux']*np.sqrt(f[1].data['ivar']))
	
	plt.title('%s (S/N = %.1f)' % (target, snr))
	plt.savefig(outdir + 'plots/fluxed_%s.png' % target)
	plt.close()
	f.close()

print('finished making coadd previews...')

###############################################################################################
# END
###############################################################################################

end = time.time()
time_min = (end - start) / 60

print('                      ')
print(bcolors.HEADER + '##################################' + bcolors.ENDC)
print(bcolors.HEADER + '        REDUCTION COMPLETE        ' + bcolors.ENDC)
print(bcolors.HEADER + '        ELAPSED: %.1f MINS        ' % time_min + bcolors.ENDC)
print(bcolors.HEADER + '##################################' + bcolors.ENDC)
print('                      ')
