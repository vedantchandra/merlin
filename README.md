<p align="center">
  <img src="merlin.png" alt="logo" width="400" align="center">
</p>

# merlin

End-to-end reduction and analysis pipeline for Magellan/MagE echellette spectra
of outer-halo stars. Raw frames are reduced with PypeIt (v1.15.0), matched to
all-sky photometry, and fit star-by-star with MINESweeper to produce an "rcat"
of stellar parameters, distances, and Galactic kinematics.

Everything runs on the Harvard FAS Cannon cluster via SLURM. This README
describes the pipeline as it exists in September 2026.

## Overview

The pipeline has two halves:

1. **Reduction** (`radagast.py`): one SLURM job per observing night turns raw
   MagE frames into flux-calibrated, order-stitched 1D coadds using PypeIt.
2. **Analysis** (`01`–`06`): coadds are collated into a survey catalog
   ("spall"), cross-matched to Gaia EDR3 / PS1 / 2MASS / unWISE / SDSS
   photometry ("acat"), fit one-star-per-job with MINESweeper, and the fit
   results are collated into a results catalog ("rcat").

```
raw frames (per night)
   │  00_reduce.py  →  radagast.py   [pypeit2 env, 8 cores/night]
   ▼
data/<night>/reduced_v0/magellan_mage_A/Science/coadd/<target>_coadd.fits
   │  01_make_spall.py
   ▼
catalogs/spall.fits            (one row per coadd, joined to target DBs)
   │  02_xmatch_gall.sh → 02_xmatch_gall.py   [SLURM array, 27 chunks]
   ▼
catalogs/xgall/mage_xgall_<i>.fits
   │  03_merge_acat.py
   ▼
catalogs/mage_acat.fits        (spall + photometry, ACAT_ID assigned)
   │  05_runms_cat.py → 04_runms_star.py   [SLURM array, one star per task]
   │       runMS.py → calcpars.py → compmod.py → corner.py
   ▼
samples/mage/<VER>/*_samp.dat, pars/mage/<VER>/*.pars, plots/mage/<VER>/
   │  06_mkrcat.py
   ▼
catalogs/mage_rcat_<VER>_MSG.fits
```

## Data layout

All data live on lab storage, not in this repo:

```
/n/holystore01/LABS/conroy_lab/Lab/vchandra/mage/
├── data/
│   ├── <YYYY_MM_DD>/raw/          raw MagE frames for one night (survey nights)
│   ├── <YYYY_MM_DD>/reduced_v0/   PypeIt output for that night
│   ├── b<YYYY_MM_DD>/             same, for Ana Bonaca's program nights
│   ├── reduced/v0/                all coadds copied here as <night>_<target>.fits
│   └── *.zip                      raw data not yet unpacked/reduced
├── catalogs/
│   ├── tdb/targetdb_*.fits        target databases per semester (2022b–2024a, bonaca)
│   ├── spall.fits                 collated observation log
│   ├── xgall/                     per-chunk photometry crossmatch outputs
│   ├── mage_acat.fits             analysis catalog (input to MINESweeper)
│   └── mage_rcat_<VER>_MSG.fits   results catalogs
├── samples/<cat>/<VER>/           dynesty posterior samples per star
├── pars/<cat>/<VER>/              one-row summary per star (.pars)
├── plots/<cat>/<VER>/             best-fit spectrum / SED / corner plots
└── logs/                          SLURM stdout/stderr
```

The repo itself lives at `/n/home03/vchandra/outerhalo/08_mage/` and all
scripts hard-code both that path and the data path above.

The `99_store_to_scratch.sh` / `99_scratch_to_store.sh` scripts are relics of
when the data lived on `holyscratch01`. Everything has since been moved to
`holystore01` (commits of April 2024), and the SLURM templates now point there.

## Environments

Two conda environments are used:

| Env       | Used by                                | Notes                                   |
|-----------|----------------------------------------|-----------------------------------------|
| `pypeit2` | `radagast.py` (reduction)              | PypeIt v1.15.0; see `environment.yml`   |
| `outerhalo` | everything else (MINESweeper, catalogs) | needs `minesweeper`, `gala`, astropy   |

MINESweeper model files are expected under `~/software/MS_files/`:
the R=12k LinNet spectral and continuum ANNs (`NN/R12K/modV0_*_LinNet_R12K_WL445_565.h5`),
the photometric ANN directory (`VARRV/`), the MIST grid
(`MIST_2.0_spot_EEPtrk_small.h5`), and the CKC library used by `calcpars.py`.

## Pipeline stages

### 00 — Reduce a night (`00_reduce.py` → `radagast.py`)

```bash
cd pipeline
python 00_reduce.py --dir=2024_04_23 --version=0        # one night
python 00_reduce.py --dir=all --version=0               # every data/202* night
python 00_reduce.py --dir=bonaca --version=0            # every data/b202* night
python 00_reduce.py --dir=<night> --dryrun=1            # only build obslog + .pypeit file
python 00_reduce.py --dir=<night> --skipred=1           # keep cals, redo flux/coadd only
```

`00_reduce.py` fills `slurm/reduce_template.txt` and submits one 8-core job
per night, which runs `radagast.py`. Radagast does the following for a night:

1. Gunzips raw frames, builds an obslog with `pypeit_obslog`, and un-comments
   any frames PypeIt mistyped.
2. Classifies frames by target name: `j*` = science, `hip*`/`ltt*` = standard,
   `thar`/`arc` = arc, `flat`/`flash` = trace+flats. Only the last exposure of
   each standard is kept.
3. Pairs every science/standard frame with the nearest **following** ThAr arc
   and assigns calib groups. Frames with the same target name get the same
   `comb_id` so PypeIt coadds them.
4. Writes `obslog_edited.txt` and a `.pypeit` file with custom parameters
   (no bias, no flexure correction, `edge_thresh=1`, `snr_thresh=3`,
   `maxnumber_sci=1`, `max_mask_frac=1.0`), then runs `run_pypeit`.
5. Picks the flux standard from the `flux_standards` list at the top of the
   file. **If you observe a new standard star, add it there** or the script
   raises. The sensfunc block assumes an A0 star with `star_mag=7.27`.
6. "Fudge-fixes" order 6 of the standard's spec1d (copies boxcar extraction
   into the optimal columns) so `pypeit_sensfunc` doesn't fail on a missing
   order.
7. Runs `pypeit_sensfunc`, `pypeit_flux_calib`, and `pypeit_coadd_1dspec` per
   target (velocity-grid wave method), then writes preview PNGs to
   `reduced_v0/magellan_mage_A/plots/`.

Nights that fail are usually fixed by hand: adding a standard to the library,
renaming a mis-labelled frame in the special-case block, or the
`nb/99_fix_files.ipynb` notebook. `nb/99_check_reductions.ipynb` summarises
which nights/targets have coadds.

### 01 — Build a target database for new nights (`01_make_tdb_2026.py`)

Before `01_make_spall.py` can keep a coadd, its target name must exist in one
of the `catalogs/tdb/targetdb_*.fits` files (columns `source_id, ra, dec,
name, selection`). `01_make_tdb_2026.py` builds `targetdb_2026a.fits` for the
two 2026 Bonaca nights and is the template for future nights:

1. Takes the `name` of every science/standard frame from the FITS `OBJECT`
   header (lower-cased), so the join in make_spall is guaranteed to hit.
2. Looks each name up in the night's observing plan for an approximate
   position, G magnitude, and selection tag. Plan coordinates are **not**
   written to the tdb: the February plan lists RA to 1 s (~11"), far coarser
   than the 3" photometry match.
3. Runs batched Gaia DR3 cone searches (25" around the plan position,
   constant centres so the archive's index is used) and adopts the source
   whose G matches the plan within 0.4 mag (1.5 mag for HIP standards),
   nearest first. This gives the precise `ra, dec` and `source_id`.
4. Cross-checks against the FITS header pointing (catches mislabeled frames)
   and against any existing tdb rows (catches method errors), and prints every
   anomaly to `targetdb_2026a_issues.txt` with a full per-target table in
   `targetdb_2026a_diagnostics.csv`.

Run it in the `outerhalo` env; it needs network access to the Gaia archive.
Known 2026 quirks it handles: the frame labelled `hip21020` on 2026-02-07 was
actually pointed at hip21024 (kept under the header name with hip21024's
coordinates), and the plan's `j1200m2755`/`j1200m2755b` pair 11" apart
resolves to the unsuffixed star.

### 01 — Collate coadds (`01_make_spall.py`)

Copies every `*_coadd.fits` and preview PNG into `data/reduced/v0/` and
`plots/v0/`, named `<night>_<target>`. Builds `spall.fits` with one row per
coadd (name, date, header keys prefixed `mage_`), then left-joins to the union
of the target DBs (`tdb_` columns: RA, Dec, selection, etc.). Rows with no
target-DB match are dropped and printed. Use `--no-transfer` to skip the rsync
step and only rebuild the table.

Note the special case `j2035m2245 → j2035m2445` for a mis-named header.
The list of target DBs is hard-coded in this script; add new ones there.

### 02 — Photometry crossmatch (`02_xmatch_gall.sh`)

```bash
sbatch 02_xmatch_gall.sh
```

A 27-task SLURM array; each task loads one chunk of Charlie Conroy's Gaia EDR3
"gall2" all-band photometry catalog and sky-matches spall to it (3 arcsec, or
20 arcsec for `rvs`/`tell` targets). Output: `catalogs/xgall/mage_xgall_<i>.fits`.

### 03 — Merge into the acat (`03_merge_acat.py`)

Stacks the xgall chunks, drops a few unWISE columns, sorts by MJD, and assigns
`ACAT_ID` (a running integer). Output: `catalogs/mage_acat.fits`.
`ACAT_ID` is the index used by the SLURM array in the next step, so
**re-running this script renumbers stars**.

### 04/05 — MINESweeper fits (`05_runms_cat.py` → `04_runms_star.py`)

```bash
python 05_runms_cat.py --catalog=mage --version=V0.08 --npoints=500
python 05_runms_cat.py --catalog=mage --version=V0.08 --overwrite=1 --arr1=0 --arr2=100
python 05_runms_cat.py --catalog=mage --version=V0.08 --skipfit=1     # redo postprocessing only
python 05_runms_cat.py --catalog=mage --version=V0.08 --no-submit     # just write the .sh
```

`05_runms_cat.py` fills `slurm/runms_template.txt` and writes
`slurm/04_runms_<catalog>_<VER>.sh`. By default it looks in
`pars/<catalog>/<VER>/` and only submits stars without a `.pars` file; with
`--overwrite=1` it submits the index range `arr1–arr2`. Each array task runs
`04_runms_star.py --ind=$SLURM_ARRAY_TASK_ID`, which for one acat row:

- `runMS.py` — builds the MINESweeper input: photometry with per-survey floors
  and cuts, spectrum from `getdata.py`, priors (uniform EEP/mass/[Fe/H]/[a/Fe],
  Gaia parallax, SFD-based Av, truncated-Gaussian age, beta Vrot, ±500 km/s
  Vrad), and runs dynesty (`rwalk`, `npoints` live points, 25 walks). Writes
  `samples/<cat>/<VER>/mage_<GaiaID>_<date>_<VER>_samp.dat`.
- `calcpars.py` — summarises the posterior (median, 16/84 percentiles), computes
  best-fit spectrum/SED chi-square and S/N, and runs `phaseafy.py` to get
  Galactocentric positions, velocities, angular momenta, and energies in the
  default `potentials.py` Milky Way model (via gala). Writes the `.pars` file.
- `compmod.py` — best-fit spectrum + SED comparison plot.
- `corner.py` — corner plot with priors overlaid.

`getdata.py` is what defines the fitted spectrum: it reads the coadd from
`spall['specfile']`, keeps 4800–5500 Å with good mask/ivar, rejects 10σ
outliers against a 5-pixel median filter, normalises by the median flux, and
sets the instrumental resolution from the linear `control/res_sigma_p.txt` fit
(measured in `nb/01_measure_res.ipynb`) times a 0.5 fudge factor.
`control/redux.txt` holds the reduction version string (`v0`).

The `--sel` flag on `05_runms_cat.py` is parsed but not yet used.

### 06 — Build the rcat (`06_mkrcat.py`)

```bash
python 06_mkrcat.py --catalog=mage --version=V0.08
```

Copies the acat, adds the ~150 MINESweeper/kinematic columns, fills them from
each star's `.pars` file (stars without one stay NaN), renames columns to be
IDL-safe (`[Fe/H]→FeH`, `log(g)→logg`, `initial_→init_`), computes Sgr
coordinates and a Sgr flag, and writes `catalogs/<cat>_rcat_<VER>_MSG.fits`.
A `FLAG` column exists but is currently all −1; the duplicate-flagging logic
is commented out.

## Helper modules

| File            | Purpose                                                        |
|-----------------|----------------------------------------------------------------|
| `getdata.py`    | Load one star's photometry row + cleaned spectrum from the acat |
| `phaseafy.py`   | Posterior samples → Galactic phase-space quantities (gala)      |
| `potentials.py` | Milky Way potential definitions used by `phaseafy`              |
| `photsys.py`    | Photometric system / filter-curve lookup                        |
| `star_basis.py` | Stellar spectral basis interpolation (used by `compmod`)        |
| `ccm_curve.py`  | CCM extinction curve                                            |
| `quantiles.py`  | Weighted quantiles                                              |
| `99_runstar.py` | Older standalone MINESweeper runner (pre-acat); superseded      |

## Notebooks (`pipeline/nb/`)

| Notebook                   | What it does                                                        |
|----------------------------|---------------------------------------------------------------------|
| `00_survey_overview`       | Survey progress / sky maps / comparisons to H3 (figures in `nb/fig/`) |
| `01_measure_res`           | Measure MagE resolution vs wavelength from arcs → `control/res_sigma_p.txt` |
| `02_make_acats`            | Inspect/subset the acat (e.g. the `rvs` RV-standard sub-catalog)     |
| `03_preview_cals`          | Look at raw frames / calibrations for one night before reducing      |
| `04_preview_chains`        | Inspect dynesty chains for a version                                 |
| `05_h3cal`                 | Compare MagE fits of H3 stars to the H3 rcat (`h3cal` catalog)       |
| `06_preview_fits`          | Plot best-fit models against spectra                                 |
| `99_check_reductions`      | Per-night audit of what reduced and what failed                      |
| `99_collate_coadds`, `99_dev_getspec`, `99_fix_files`, `99_newpsc` | scratch / one-off fixes |

## SLURM conventions

- Templates live in `pipeline/slurm/*_template.txt`; generated scripts are
  written next to them and committed as a record of what was run.
- Reduction: 8 cores, 4.5 GB/core, 5 h, partitions `conroy_priority,shared,itc_cluster`.
- MINESweeper: 1 core, 10 GB, 15 h, partitions `conroy,shared,itc_cluster,sapphire`,
  Intel nodes only.
- Logs go to `<datadir>/logs/reduce/` and `<datadir>/logs/<catalog>/<VER>/`.

## Current state (September 2026)

- All 30 survey nights (2022_08_03 → 2024_04_23) and the four 2024 Bonaca
  nights are reduced at `v0`; there are 782 coadd files in `data/reduced/v0/`.
  Five older `b2022_*` Bonaca nights have raw data but no reduction.
- MINESweeper versions `V0.07` and `V0.08` have `.pars` outputs on disk;
  rcats exist for `V0.0`, `V0.01`, `V0.03`, `V0.07`, and `V0.08`. `V0.08` is
  the latest.
- The `ut2403*`/`ut2404*`/`ut2405*` zips in `data/` are the source archives
  for the four reduced `b2024_*` nights and can be ignored.
- Two new Bonaca-program nights were unpacked on 2026-09-16 and pass a
  radagast dry run but are **not yet reduced**: `b2026_02_07` (128 frames,
  standard hip70050) and `b2026_04_06` (137 frames, standard hip44395). The
  source zips (`ut260207_08-*.zip`, `20260406-*.zip`) remain in `data/`.
  Observing logs and target lists from the archives sit at the night level
  next to `raw/`. Three unnamed startup ThAr frames from the February night
  were moved to `b2026_02_07/raw_unused/`.
- `catalogs/tdb/targetdb_2026a.fits` (68 targets, built by
  `01_make_tdb_2026.py` from the observing plans + Gaia DR3) covers both 2026
  nights and is loaded by `01_make_spall.py`. Selections: hvs, jet, ngc1851,
  ngc5904, rvs, tell.
- A separate `rvs` catalog (RV standards) and `h3cal` catalog (H3 overlap
  stars) have been fit for calibration.

## Known issues / to do

- **Duplicate night `b2024_03_31`.** The 2024-03-31 data were ingested twice:
  once as survey night `2024_03_31` (reduced April 2024) and again as Bonaca
  night `b2024_03_31` (reduced August 2025, from `ut240330_31.zip`, with two
  extra ThAr frames). The raw frames are byte-identical and all 19 targets
  appear twice in spall/acat/rcat. Only the `2024_03_31` copies have been
  fitted. Fix: exclude `b2024_03_31` in `01_make_spall.py`, then rerun steps
  01 → 03 and rebuild the rcat.
- HIP standards within ~5° of the Galactic plane never get photometry
  because gall2 has no coverage there (hip38789, hip41926, hip42922,
  hip43656, hip44395, hip84267, hip84881). They are dropped at the acat
  stage. Harmless for science targets.
- `radagast.py` only fudge-fixes missing orders on the flux standard, so other
  standards with a missing order fail the 1D coadd (e.g. hip17819 and
  hip22865 on b2026_02_07).

## Adding a new night

1. Unpack raw frames to `<datadir>/data/<YYYY_MM_DD>/raw/`.
2. Check the standard star is in `flux_standards` in `radagast.py`.
3. `python 00_reduce.py --dir=<night> --dryrun=1`, inspect the obslog and
   `.pypeit` file, then run for real.
4. Add the night's targets to a target DB under `catalogs/tdb/` and load it in
   `01_make_spall.py` if it is a new semester.
5. Re-run steps 01 → 03, then `05_runms_cat.py` (it will only fit new stars),
   then `06_mkrcat.py`.

## Other directories

- `science/01_process_ngc19.py` — converts N-body tracer particle files for the
  NGC 19 simulation into observed-frame coordinates (unrelated to MagE
  reduction; lives here for convenience).
- `old/` — the 2022 pre-pipeline scripts and their SLURM logs, kept for
  reference.
