## Purpose of script:
`part6.model_grid_registration.py` registers the input model to the 3D mesh made in `part1` via linearly interpolating the eight nearest model values to each 3D block node. It converts the velocity values in the input model to percent perturbation relative to the reference model (default: PREM). Finally, it computes a radial RMS profile for the perturbation values in each depth shell.


## Output/saved files:
1. `models/{data_wave_type}/{model}_update/{model}_grid_registered.csv`: the newly interpolated model, ready to be used for the model update.
2. `models/{data_wave_type}/{model}_update/{model}_RMS_{perturb}.csv`: the RMS(perturb) profile of the model.


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. make sure that original path files for all of the residuals in the raw dataset have been made using `part2.parse_data_make_raypaths.py`.
3. make sure that resampled path files for all of the residuals in the raw dataset have been made using `part3.resample.py`.
4. make sure that your input model is in a CSV format and is saved in the correct directory (`models/{data_wave_type}/{model}_update/{model}.csv`)
5. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `input_model`: the name of the input model file before the `.csv` extension
2. `data_wave_type`: the polarity of the waves from which the model was made (i.e., `'S'` or `'P'`.)
3. `lat_header`: the latitude header in the input `{model}.csv` file
4. `lon_header`: the longitude header in the input `{model}.csv` file
5. `depth_header`: the depth header in the input `{model}.csv` file
6. `property_header`: the velocity value header in the input `{model}.csv` file
7. `delimiter`: the field delimiter in the input `{model}.csv` file (typically `|`)
8. `voigt`: indication of whether or not to compute a voigt average (typically done for anisotropic models that give velocities for both horizontal and vertical components). if voigt averaging is desired, the user must indicate the headers of the values to average.


## Necessary compute resources:
1. HPC cores: 1
2. Time allocation: Generally about an hour