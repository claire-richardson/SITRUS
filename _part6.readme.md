## Purpose of script:
`part6.model_grid_registration.py` registers the input model to the 3D mesh made in `part1` by taking a weighted average of the original values closest to the 3D block node. If indicated in `mod_input`, it converts the velocity values in the input model to percent perturbation relative to the reference model (default: PREM). Finally, it computes a radial RMS profile for the perturbation values in each depth shell.


## Output/saved files:
1. `models/{data_wave_type}/{model}_update/{model}_grid_registered.csv`: the newly interpolated model, ready to be used for the model update.
2. `models/{data_wave_type}/{model}_update/{model}_RMS_{perturb}.csv`: the RMS(perturb) profile of the model.


## Before running:
1. make sure that a mesh has been made using `part1.grid_definition.py`.
2. make sure that your input model is in a CSV format and saved in the correct directory (`models/{data_wave_type}/{model}_update/{model}.csv`)
3. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `input_model`: the name of the input model file before its `.csv` extension
2. `data_wave_type`: the polarity of the waves from which the model was made (i.e., `'S'` or `'P'`.)
3. `lat_header`: the latitude header in the input `{model}.csv` file
4. `lon_header`: the longitude header in the input `{model}.csv` file
5. `depth_header`: the depth header in the input `{model}.csv` file
6. `property_header`: the model value header in the input `{model}.csv` files
7. `voigt`: indication of whether or not to compute a voigt average (typically done for anisotropic models that give velocities for both horizontal and vertical components). if voigt averaging is desired, the user must indicate the headers of the values to average.
8. `convert_vel_to_perturb`: boolean value indicating whether or not the model values need to be converted from absolute velocities to velocity perturbations relative to PREM. The model update only handles perturbations, so if the model values are given in velocities, this conversion must be performed.


## Necessary compute resources:
1. HPC cores: as many depth shells as there are in your model space/`shell_dimensions.csv`
2. Time allocation: Can vary between an hour and a day, depending on the number of nodes in the input model.
3. Memory allocation: Can require up to tens of GB, if the input model file is very large.