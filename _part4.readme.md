## Purpose of script:
`part4.data_coverage.py` calculates the total raypath coverage for each individual phase as well as the entire dataset. For each 3D block in the mesh, it counts the total overall number of paths that sample it, the total number of azimuthal sectors covered, and the total number of paths in each azimuthal sector.


## Output/saved files:
1. makes and saves individual coverage files for each phase as `phases/{phase}/phase_data/{phase}_coverage.csv`
2. makes and saves a total coverage file for the whole dataset as `phases/total_coverage.csv`
3. outputs progress to `phases/{phase}/phase_data/{phase}_pt4_coverage_log.txt`


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. make sure that original path files for all of the residuals in the raw dataset have been made using `part2.parse_data_make_raypaths.py`.
3. make sure that resampled path files for all of the residuals in the raw dataset have been made using `part3.resample.py`.
4. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `all_phases`: a list of all phases in the raw dataset.
2. `azimuthal_sectors`: the total number of equal area azimuthal sectors, defined from 0-180 degrees, you would like to calculate coverage for and subsequently use for azimuthal weighting in the model update.


## Necessary compute resources:
1. HPC cores: as many as the number of phases in `all_phases`.
2. Time allocation: Generally less than an hour, depending on the total number of data and phases.