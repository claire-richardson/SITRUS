## Purpose of script:
`part2.parse_data.py` establishes the rest of the file system and parses the raw dataset into phase-wise subdirectories. It then calculates 1D raypaths for each of the residuals in the raw dataset using the TauP toolkit (Crotwell et al., 1999; https://github.com/crotwell/TauP). If indicated in `mod_input.py`, the following are applied to all residuals:
    - ellipticity corrections (EllipticiPy; Russell et al., 2022)
    - crustal corrections (CRUST 1.0; Laske et al., 2013)
    - any other desired special corrections however these must be configured by the user


## Output/saved files:
1. 

## Before running:
1. run `part1.grid_definition.py`
2. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`:
1. `dataset`: raw dataset of travel time residuals. Must be a CSV file, and must contain at least the following information for each residual: 1) phase; 2) event latitude; 3) event longitude; 4) station latitude; 5) station longitude; 6) event depth; 7) travel time residual. These must be columns, and their headers must match the following: 'PHASE', 'EQ_LAT', 'EQ_LON', 'STA_LAT', 'STA_LON', 'EQ_DEP', 'DT' in order to be properly interpreted.  If there are more columns, they will just be ignored--no need to delete manually. If any of the other columns are special weights to be added, that can be indicated in a separate `mod_input.py` variable.
2. 

## Necessary compute resources:
1. HPC cores: as many as the number of phases in `dataset`
2. Time allocation: Generally a few hours, depending on the total number of data and corrections to be applied.