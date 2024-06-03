## Purpose of script:
`part6.model_grid_registration.py` registers the input model to the 3D mesh made in `part1` by taking a weighted average of the original values closest to the 3D block node. It converts the velocity values in the input model to percent perturbation relative to the reference model (default: PREM). Finally, it computes a radial RMS profile for the perturbation values in each depth shell.


## Output/saved files:
1. `models/{data_wave_type}/{model}_update/{model}_grid_registered.csv`: the newly interpolated model, ready to be used for the model update.
2. `models/{data_wave_type}/{model}_update/{model}_RMS_{perturb}.csv`: the RMS(perturb) profile of the model.


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. make sure that your input models are in a CSV format and saved in the correct directory (`models/{data_wave_type}/{model}_update/{model}.csv`)
3. make sure that the velocity value in your input models is labled with either `vsh` or `vsv`.
4. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `all_models_to_process`: a list of the names of the input model files before their `.csv` extensions
2. `data_wave_type`: the polarity of the waves from which the models were made (i.e., `'S'` or `'P'`.)
3. `lat_header`: the latitude header in the input `{model}.csv` files
4. `lon_header`: the longitude header in the input `{model}.csv` files
5. `depth_header`: the depth header in the input `{model}.csv` files
6. `property_header`: the velocity value header in the input `{model}.csv` files (either `vsh` or `vsv`
7. `delimiter`: the field delimiter in the input `{model}.csv` file (default `,`; will be the delimiter if `part5` was used for the model conversion)
8. `voigt`: indication of whether or not to compute a voigt average (typically done for anisotropic models that give velocities for both horizontal and vertical components). if voigt averaging is desired, the user must indicate the headers of the values to average.


## Necessary compute resources:
1. HPC cores: as many models as you plan to process (i.e., the length of the `mod_input.all_models_to_process` list.
2. Time allocation: Generally about an hour