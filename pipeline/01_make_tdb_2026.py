"""
Build catalogs/tdb/targetdb_2026a.fits for the two Bonaca-program MagE nights
b2026_02_07 and b2026_04_06.

Sources of truth, in order:
  1. the FITS OBJECT header of every science/standard frame (defines `name`)
  2. the per-night observing plans (name -> approximate RA/Dec, G, selection)
       b2026_02_07/Observing plan for MagE on Feb 7 2026 - Plan Overview.csv
       b2026_04_06/OC_MagE_xue_mage20260406_cat_apr6_v1.txt
  3. a Gaia DR3 cone search around the plan position, filtered on G, which
     supplies the precise coordinates and source_id written to the tdb
  4. existing targetdb_*.fits rows, used only to validate (3)

Every anomaly is printed and written to targetdb_2026a_issues.txt.
Run in the `outerhalo` env (needs astroquery + network).
"""
import re, glob, sys, csv
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astropy.coordinates import SkyCoord
import astropy.units as u

D   = '/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/'
TDB = D + 'catalogs/tdb/'
OUT = TDB + 'targetdb_2026a.fits'
ISS = TDB + 'targetdb_2026a_issues.txt'

NIGHTS = {
    'b2026_02_07': D + 'data/b2026_02_07/Observing plan for MagE on Feb 7 2026 - Plan Overview.csv',
    'b2026_04_06': D + 'data/b2026_04_06/OC_MagE_xue_mage20260406_cat_apr6_v1.txt',
}
RAD_SCI, RAD_STD = 15.0, 20.0     # cone radius [arcsec]; plan RA is only good to 1 s in Feb
DG_SCI,  DG_STD  = 0.4, 1.5       # |G_gaia - G_plan| tolerance (HIP stars are near saturation)

issues = []
def flag(msg):
    issues.append(msg); print('  !! ' + msg)

# ------------------------------------------------------------------ plans
def parse_feb(path):
    rows = []
    with open(path) as f:
        for r in csv.reader(f):
            if len(r) < 6 or not re.match(r'^(j\d{4}[pm]\d{4}b?|hip\d+)$', r[1].strip().lower()):
                continue
            rows.append(dict(name=r[1].strip().lower(), ra=r[2].strip(), dec=r[3].strip(),
                             gmag=float(r[4]), selection=r[5].strip().lower(), night='b2026_02_07'))
    return rows

def parse_apr(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        body, _, comment = line.partition('#')
        f = body.split()
        m = re.search(r'G\s*=\s*([\d.]+)\s*,\s*(\w+)', comment)
        rows.append(dict(name=f[1].lower(), ra=f[2], dec=f[3], gmag=float(m.group(1)),
                         selection=m.group(2).lower(), night='b2026_04_06'))
    return rows

plan = parse_feb(NIGHTS['b2026_02_07']) + parse_apr(NIGHTS['b2026_04_06'])
plan = Table(plan)
c = SkyCoord(plan['ra'], plan['dec'], unit=(u.hourangle, u.deg))
plan['ra_plan'], plan['dec_plan'] = c.ra.deg, c.dec.deg
print('plan entries: %i (Feb %i, Apr %i)' % (len(plan), np.sum(plan['night']=='b2026_02_07'), np.sum(plan['night']=='b2026_04_06')))

# same name in both plans -> must agree
for nm in set(plan['name']):
    s = plan[plan['name'] == nm]
    if len(s) > 1:
        sep = SkyCoord(s['ra_plan'], s['dec_plan'], unit='deg')
        d = sep[0].separation(sep[1]).arcsec
        if d > 10 or len(set(s['selection'])) > 1:   # Feb plan RA is rounded to 1 s (~7\" max)
            flag('%s appears in both plans with sep=%.1f" and selections %s' % (nm, d, list(s['selection'])))

# ------------------------------------------------------------------ observed frames
obs = []
for night in NIGHTS:
    for fn in sorted(glob.glob(D + 'data/%s/raw/mage*.fits' % night)):
        h = fits.getheader(fn)
        nm = str(h['OBJECT']).strip().lower()
        if not (nm.startswith('j') or nm.startswith('hip') or nm.startswith('ltt')):
            continue
        hc = SkyCoord(h['RA'], h['DEC'], unit=(u.hourangle, u.deg))
        obs.append(dict(name=nm, night=night, file=fn.split('/')[-1], ra_hdr=hc.ra.deg, dec_hdr=hc.dec.deg))
obs = Table(obs)
print('observed science/standard frames: %i' % len(obs))

# one row per (name, night); keep first header coordinate
targets = unique(obs, keys=['name', 'night'])
print('unique observed targets: %i' % len(targets))

# ------------------------------------------------------------------ join to plan
targets['ra_plan'] = np.nan; targets['dec_plan'] = np.nan; targets['gmag'] = np.nan
targets['selection'] = np.array(['unknown'] * len(targets), dtype='U16')
for r in targets:
    p = plan[(plan['name'] == r['name']) & (plan['night'] == r['night'])]
    if len(p) == 0:
        p = plan[plan['name'] == r['name']]
        if len(p): flag('%s (%s) not in that night\'s plan; using the other night\'s entry' % (r['name'], r['night']))
    if len(p) == 0:
        flag('%s (%s) is NOT in any plan; falling back to FITS header coordinates' % (r['name'], r['night']))
        r['ra_plan'], r['dec_plan'] = r['ra_hdr'], r['dec_hdr']
        continue
    r['ra_plan'], r['dec_plan'], r['gmag'], r['selection'] = p['ra_plan'][0], p['dec_plan'][0], p['gmag'][0], p['selection'][0]
    # header pointing vs plan: catches mislabeled frames
    d = SkyCoord(r['ra_hdr'], r['dec_hdr'], unit='deg').separation(SkyCoord(p['ra_plan'][0], p['dec_plan'][0], unit='deg')).arcsec
    if d > 30:
        near = plan[SkyCoord(plan['ra_plan'], plan['dec_plan'], unit='deg').separation(SkyCoord(r['ra_hdr'], r['dec_hdr'], unit='deg')).arcsec < 30]
        flag('%s (%s, %s): header pointing is %.0f" from its plan position; nearest plan entries within 30": %s'
             % (r['name'], r['night'], r['file'], d, list(near['name'])))
        if len(near) == 1:
            flag('   -> assuming the telescope was on %s; using that plan entry (name kept as %s to match the FITS header)' % (near['name'][0], r['name']))
            r['ra_plan'], r['dec_plan'], r['gmag'], r['selection'] = near['ra_plan'][0], near['dec_plan'][0], near['gmag'][0], near['selection'][0]

# ------------------------------------------------------------------ Gaia DR3
from astroquery.gaia import Gaia
targets['rad'] = np.where([n.startswith('hip') for n in targets['name']], RAD_STD, RAD_SCI)
# batched cone searches with constant centres (uses the archive's spatial index; a joined
# upload-table crossmatch against gaia_source times out on the sync endpoint)
print('querying Gaia DR3 for %i targets in batches...' % len(targets))
cands = []
B = 15
for i in range(0, len(targets), B):
    chunk = targets[i:i+B]
    circles = ' OR '.join("1=CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.7f, %.7f, %.6f))"
                          % (r['ra_plan'], r['dec_plan'], 25./3600.) for r in chunk)
    adql = ("SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag AS g_gaia, g.pmra, g.pmdec, g.parallax "
            "FROM gaiadr3.gaia_source AS g WHERE " + circles)
    res = Gaia.launch_job(adql).get_results()
    print('  batch %i-%i: %i sources' % (i, i+len(chunk)-1, len(res)))
    cands.append(res)
cand = vstack(cands)
for c in list(cand.colnames): cand.rename_column(c, c.lower())   # archive may return upper-case names
cand = unique(cand, keys='source_id')
cc_all = SkyCoord(cand['ra'], cand['dec'], unit='deg')
print('  %i unique candidate sources' % len(cand))

targets['source_id'] = np.int64(999999); targets['ra'] = np.nan; targets['dec'] = np.nan
targets['g_gaia'] = np.nan; targets['sep_plan'] = np.nan; targets['ncand'] = 0
for r in targets:
    cc = cand.copy()
    cc['sep'] = SkyCoord(r['ra_plan'], r['dec_plan'], unit='deg').separation(cc_all).arcsec
    cc = cc[cc['sep'] < 25.]
    dg = DG_STD if r['name'].startswith('hip') else DG_SCI
    ok = cc[(cc['sep'] < r['rad']) & (np.abs(cc['g_gaia'] - r['gmag']) < dg)] if np.isfinite(r['gmag']) else cc[cc['sep'] < r['rad']]
    r['ncand'] = len(ok)
    if len(ok) == 0:
        flag('%s (%s): NO Gaia DR3 source within %.0f" with |dG|<%.1f (G_plan=%.1f); %i sources in cone: %s'
             % (r['name'], r['night'], r['rad'], dg, r['gmag'], len(cc), [(round(float(s), 1), round(float(g), 1)) for s, g in zip(cc['sep'], cc['g_gaia'])]))
        r['ra'], r['dec'] = r['ra_plan'], r['dec_plan']
        continue
    if len(ok) > 1:
        flag('%s (%s): %i Gaia sources pass; taking nearest. (sep",G): %s'
             % (r['name'], r['night'], len(ok), [(round(float(s), 1), round(float(g), 2)) for s, g in zip(ok['sep'], ok['g_gaia'])]))
    b = ok[np.argmin(ok['sep'])]
    r['source_id'], r['ra'], r['dec'], r['g_gaia'], r['sep_plan'] = b['source_id'], b['ra'], b['dec'], b['g_gaia'], b['sep']

# ------------------------------------------------------------------ validate against existing tdbs
old = vstack([Table.read(f) for f in sorted(glob.glob(TDB + 'targetdb_20*.fits')) + [TDB + 'targetdb_bonaca.fits'] if 'targetdb_2026' not in f])
old['name'] = [str(x).strip() for x in old['name']]
old = unique(old, keys='name')
print('\nvalidation against existing target DBs:')
nval = 0
for r in targets:
    o = old[old['name'] == r['name']]
    if len(o) == 0: continue
    nval += 1
    d = SkyCoord(r['ra'], r['dec'], unit='deg').separation(SkyCoord(o['ra'][0], o['dec'][0], unit='deg')).arcsec
    same_id = (o['source_id'][0] == r['source_id']) or o['source_id'][0] == 999999
    print('  %-12s old sep=%.2f" source_id %s sel old=%s new=%s' % (r['name'], d, 'same' if same_id else 'DIFFERENT', o['selection'][0], r['selection']))
    if d > 1.0 or not same_id:
        flag('%s: disagrees with existing tdb (sep=%.2f", id same=%s)' % (r['name'], d, same_id))
    if str(o['selection'][0]).strip() != r['selection']:
        flag('%s: selection differs from existing tdb (%s vs %s); keeping plan value %s' % (r['name'], o['selection'][0], r['selection'], r['selection']))
print('  %i targets validated' % nval)

# header-pointing offset statistics (sanity on the header-vs-plan checks above)
hs = SkyCoord(targets['ra_hdr'], targets['dec_hdr'], unit='deg').separation(SkyCoord(targets['ra'], targets['dec'], unit='deg')).arcsec
print('\nFITS header pointing vs Gaia position: median %.1f", max %.1f"' % (np.nanmedian(hs), np.nanmax(hs)))
for r, s in zip(targets, hs):
    if s > 15: flag('%s (%s): header pointing %.0f" from adopted Gaia position' % (r['name'], r['night'], s))

# ------------------------------------------------------------------ write
# one row per name (the pipeline joins on name); if a star was observed both nights keep the first
out = unique(targets, keys='name')
tdb = Table()
tdb['source_id'] = np.array(out['source_id'], dtype=np.int64)
tdb['ra'] = np.array(out['ra'], dtype=float)
tdb['dec'] = np.array(out['dec'], dtype=float)
tdb['name'] = np.array(out['name'], dtype='S32')
tdb['selection'] = np.array(out['selection'], dtype='S11')
tdb.write(OUT, overwrite=True)
print('\nwrote %s with %i rows' % (OUT, len(tdb)))
print('selections:', dict(zip(*np.unique(tdb['selection'], return_counts=True))))
print('no Gaia match (placeholder source_id): %i' % np.sum(tdb['source_id'] == 999999))

# full diagnostic table alongside
targets[['name', 'night', 'file', 'selection', 'gmag', 'g_gaia', 'source_id', 'ra', 'dec', 'ra_plan', 'dec_plan', 'sep_plan', 'ncand']].write(TDB + 'targetdb_2026a_diagnostics.csv', overwrite=True)
with open(ISS, 'w') as f:
    f.write('\n'.join(issues) + '\n')
print('%i issues written to %s' % (len(issues), ISS))
