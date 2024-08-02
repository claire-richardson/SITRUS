## Purpose of script:
`part2.data_processing.py` is written to perform three main tasks related to data preprocessing. these are:
1. **Parse data and make original raypaths:** establishes the rest of the file system and parses the raw dataset into phase-wise subdirectories. Calculates 1D raypaths for each of the residuals in the raw dataset using the TauP toolkit (Crotwell et al., 1999; https://github.com/crotwell/TauP). If indicated in `mod_input.py`, the following are applied to all residuals:
    - ellipticity corrections (EllipticiPy; Russell et al., 2022)
    - crustal corrections (CRUST 1.0; Laske et al., 2013)
    - any other desired special corrections however these must be configured by the user

2. **Resample raypaths:** takes all of the newly made raw TauP pathfiles and performs three primary tasks with each:
    - all mesh boundaries are found wherever the raypath crosses them. these boundaries are added as a point to the raypath file. this ensures that only raypath segments that explicitly sample a given block are included during any model updates using this mesh and dataset.
    - all raypaths are resampled to have relatively evenly spaced points while retaining the boundary points. this step generally reduces the overall number of points on any given path and is performed to remove excess, redundant points, ultimately reducing overall runtime of the model update.
    - crustal corrections for CRUST1.0 (Laske et al., 2013) are calculated from the resampled paths and added to the `phases/{phase}/phase_data/{phase}_master_data.csv` files.
  
3. **Calculate coverage:** calculates the total raypath coverage for each individual phase as well as the entire dataset. For each 3D block in the mesh, it counts the total overall number of paths that sample it, the total number of azimuthal sectors covered, and the total number of paths in each azimuthal sector.


## Output/saved files:
1. establishes remainder of necessary file system for storing parsed data (`/phases` and subdirectories for each phase).
2. makes and saves all master data files for each individual phase in `phases/{phase}/phase_data`. These will include added columns with populated values for ellipticity corrections (`ELLIP_CORR`) and residuals with ellipticity corrections applied (`ELLIP_DT`). Similar columns for crustal corrections and residuals will be added (`CRUST_1.0_CORR`; `CRUST_1.0_DT`). A final column with residuals with both corrections applied will be titled `CRUST_1.0_ELLIP_DT`)
3. makes and saves all TauP raypath files in `phases/{phase}/raypath_files`.
4. outputs progress to `phases/{phase}/phase_data/{phase}_pt2_log_make_raypaths.txt`, `phases/{phase}/phase_data/{phase}_pt2_log_boundary_finding.txt`, and `phases/{phase}/phase_data/{phase}_pt2_log_coverage.txt`.
5. outputs any issues with raypath making to `phases/{phase}/phase_data/{phase}_pt2_bugs_make_raypaths.txt`.
6. outputs the file names of any paths that could not be resampled to `phases/{phase}/phase_data/{phase}_pt2_bugs_boundary_finding.txt`.
7. makes and saves all resampled raypath files in `phases/{phase}/resampled_path_files`, replacing previous original TauP pathfiles.
8. updates `phases/{phase}/phase_data/{phase}_master_data.csv` with crustal correction information
10. makes and saves individual coverage files for each phase as `coverage/{phase}/{phase}_coverage.csv`
11. makes and saves a total coverage file for the whole dataset as `coverage/total_coverage/{mod_input.dataset}.csv`


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
6. `discontinuities`: lateral discontinuities or other desired depth values to include in the resampling routine. points will be found at these depths along all raypaths.
7. `shell_bounds`: list of boundaries that define constant thickness depth shells. **this should exactly equal the list used to generate `shell_boundaries.csv` in `part1`**.
8. `target_path_length`: the target raypath segment length to resample the points along the raypath to (i.e., points along the raypath should be approximately `target_path_length` km apart from each other).
9. `target_path_length_tolerance`: the tolerance for the target raypath segment length resampling routine (i.e., points along the raypath should be approximately `target_path_length` km, +/- `target_path_length_tolerance` apart from each other).
10. `total_radius`: total radius of the Earth according to the reference model (default: PREM).
11. `azimuthal_sectors`: the total number of equal area azimuthal sectors, defined from 0-180 degrees, you would like to calculate coverage for and subsequently use for azimuthal weighting in the model update.


## Necessary compute resources:
1. HPC cores: as many as the number of phases in `all_phases` (should generally match the number of phases in `dataset`)
2. Time allocation: Generally about a day, depending on the total number of data and phases as well as the processing speed of the system.