## Purpose of script:
`part5.model_conversion.py` converts netCDF formatted files to CSV files. This utility is not strictly necessary for the SITRUS workflow if the input tomography model files are already in .csv format.

## Output/saved files:
1. converts any `{model}.nc` file to `{model}.csv` file, stored in `models/{wave_type}/{model}.csv`


## Before running:
1. upload a `{model}.nc` file to the appropriate `models/{wave_type}` directory, depending on type of waves the model uses.
2. if desired, rename your file to a name you would like to use later for the `input_model` variable in `mod_input`.


## Input variables to define in `mod_input.py`: 
1. `input_model`: the name of the model file before the `.nc` extension.
2. `data_wave_type`: the polarity of the waves from which the model was made (i.e., `'S'` or `'P'`.).


## Necessary compute resources:
1. HPC cores: 1
2. Time allocation: Generally between several seconds and ~30 minutes, depending on the size of the model file.