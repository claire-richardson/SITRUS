## Purpose of script:
`part2.parse_data.py` establishes the rest of the file system and parses the raw dataset into phase-wise subdirectories. It then calculates 1D raypaths for each of the residuals in the raw dataset using the TauP toolkit (Crotwell et al., 1999; https://github.com/crotwell/TauP). If indicated in `mod_input.py`, the following are applied to all residuals:
    - ellipticity corrections (EllipticiPy; Russell et al., 2022)
    - crustal corrections (CRUST 1.0; Laske et al., 2013)
    - any other desired special corrections however these must be configured by the user


## Output/saved files:
1. establishes remainder of necessary file system for storing parsed data (`/phases` and subdirectories for each phase).
2. makes and saves all master data files for each individual phase in `phases/{phase}/phase_data`. These will include added columns with populated values for ellipticity corrections (`ELLIP_CORR`) and residuals with ellipticity corrections applied (`ELLIP_DT`).
3. makes and saves all TauP raypath files in `phases/{phase}/raypath_files`.
4. outputs progress to `phases/{phase}/phase_data/{phase}_pt2_make_raypaths_log.txt`
5. outputs any issues with raypath making to `phases/{phase}/phase_data/{phase}_pt2_pathfile_bugs.txt`.


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. define necessary variables in `mod_input.py`, _including a list of all phases in the dataset (`mod_input.dataset`) with compatible naming convention_.
3. make sure that the columns in the raw dataset have the correct headers (see `dataset` in next section).
4. make sure that all earthquake depths (`EQ_DEP` in `mod_input.dataset`) are ≥0 (i.e., no depths are listed as being above the surface of the reference model). This may require preprocessing steps, such as applying crustal corrections, prior to using this package.


## Input variables to define in `mod_input.py`: 
1. `dataset`: file name of the raw dataset of travel time residuals. Must be a CSV file, and must contain at least the following information for each residual: 1) phase; 2) event latitude; 3) event longitude; 4) station latitude; 5) station longitude; 6) event depth; 7) travel time residual. These must be columns, and their headers must match the following: `'PHASE'`, `'EQ_LAT'`, `'EQ_LON'`, `'STA_LAT'`, `'STA_LON'`, `'EQ_DEP'`, `'DT'` in order to be properly interpreted.  If there are more columns, they will be ignored unless specified in `raw_headers_to_keep`--no need to delete manually. For example, if any of the other columns are special weights to be added or you wish to retain it for another reason, that can be indicated in `raw_headers_to_keep`. Phase names must not include integers; indicate multiples as longer strings rather than with numbers. For example, `S3` should be instead listed as `SSS`. Use `m` at the end of the entire phase name to indicate a major arc path (e.g., `ScSScSScSm`) and/or `s` at the beginning of the entire phase name to indicate a depth phase (e.g., `sS`). In the case where there are multiple datasets with the same phase, use a suffix separated by an underscore from the phase name to indicate the dataset (e.g., `Sdiff_A`, `Sdiff_B`).
2. `data_wave_type`: the polarity of the waves from which the travel time residuals were measured (i.e., `'S'` or `'P'`.).
3. `raw_headers_to_keep`: the headers of any columns you wish to retain from the raw dataset. these can include, for example, special weights that have already been calculated.
4. `reference_model`: 1D reference model you wish to use for computing raypaths and ellipticity corrections (default: PREM; currently no other models are supported).
5. `all_phases`: a list of all phases in the raw dataset, OR a list of the subset of phases in `dataset` you would like to make raypaths for (it is recommended to do all phases in `dataset`, as the user can define later which phases to include in an update).


## Necessary compute resources:
1. HPC cores: as many as the number of phases in `all_phases` (should generally match the number of phases in `dataset`)
2. Time allocation: Generally a few hours, depending on the total number of data and phases.