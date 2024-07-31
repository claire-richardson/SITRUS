## Purpose of script:
`part3.resample.py` takes all of the raw TauP pathfiles made in `part2` and performs three primary tasks with each:
1. all mesh boundaries are found wherever the raypath crosses them. these boundaries are added as a point to the raypath file. this ensures that only raypath segments that explicitly sample a given block are included during any model updates using this mesh and dataset.
2. all raypaths are resampled to have relatively evenly spaced points while retaining the boundary points. this step generally reduces the overall number of points on any given path and is performed to remove excess, redundant points, ultimately reducing overall runtime of the model update.
3. crustal corrections for CRUST1.0 (Laske et al., 2013) are calculated from the resampled paths and added to the `phases/{phase}/phase_data/{phase}_master_data.csv` files.


## Output/saved files:
1. makes and saves all resampled raypath files in `phases/{phase}/raypath_files`, replacing previous original TauP pathfiles.
2. updates `phases/{phase}/phase_data/{phase}_master_data.csv` with crustal correction information
3. outputs progress to `phases/{phase}/phase_data/{phase}_pt3_boundary_finding_log.txt`
4. outputs the file names of any paths that could not be resampled to `f'./phases/{phase}/phase_data/{phase}_pt3_resample_bugs.txt'`


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. make sure that original path files for all of the residuals in the raw dataset have been made using `part2.parse_data_make_raypaths.py`.
3. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `all_phases`: a list of all phases in the raw dataset, OR a list of the subset of phases in `dataset` you would like to resample raypaths for (it is recommended to do all phases in `dataset`, specifically those that have already been processed in `part2`, as the user can define later which phases to include in an update).
2. `data_wave_type`: the polarity of the waves from which the travel time residuals were measured (i.e., `'S'` or `'P'`.).
3. `discontinuities`: lateral discontinuities or other desired depth values to include in the resampling routine. points will be found at these depths along all raypaths.
4. `shell_bounds`: list of boundaries that define constant thickness depth shells. **this should exactly equal the list used to generate `shell_boundaries.csv` in `part1`**.
5. `target_path_length`: the target raypath segment length to resample the points along the raypath to (i.e., points along the raypath should be approximately `target_path_length` km apart from each other).
8. `target_path_length_tolerance`: the tolerance for the target raypath segment length resampling routine (i.e., points along the raypath should be approximately `target_path_length` km, +/- `target_path_length_tolerance` apart from each other).
9. `total_radius`: total radius of the Earth according to the reference model (default: PREM).


## Necessary compute resources:
1. HPC cores: as many as the number of phases in `all_phases` (should generally match the number of phases in `dataset`)
2. Time allocation: Generally a few hours, depending on the total number of data and phases.