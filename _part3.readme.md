## Purpose of script:
`part3.model_processing` performs two main tasks:
1. converts netCDF formatted files to CSV files. This utility is not strictly necessary for the SITRUS workflow if the input tomography model files are already in .csv format (indicated by `mod_input.convert_nc_to_csv`. **IMPORTANT: model values _must_ be given as velocities or percent perturbations relative to PREM. percent perturbations relative to a different reference model are not supported**
2. registers the input model to the 3D mesh made in `part1` by taking the average of the original values closest to each 3D block node. If indicated in `mod_input` (i.e., for the common case that model values are given in absolute velocity), it converts the velocity values in the input model to percent perturbation relative to PREM. It computes a radial RMS profile for the perturbation values in each depth shell. Finally, it makes plottable files of the model.


## Output/saved files:
1. creates a new directory for the model: `models/{data_wave_type}/{model}_update/`
2. converts the `{model}.nc` file to `{model}.csv` file, stored in `models/{data_wave_type}/{model}_update/{model}.csv`
3. `models/{data_wave_type}/{model}_update/{model}_grid_registered.csv`: the newly interpolated model, ready to be used for the model update.
4. `models/{data_wave_type}/{model}_update/{model}_RMS_{perturb}.csv`: the RMS(perturb) profile of the model.
5. `models/{data_wave_type}/{model}_update/original_model_plot_files/{model}_shell_{shell}_*.csv`: files that can be plotted with Cartopy in the included jupyter notebook.


## Before running:
1. check if your model file includes model values given as velocities (km / s) rather than perturbations. if there are no velocity values, or the perturbations are relative to a reference model that is not PREM, the model cannot be used.
2. upload a `{model}.nc` file to the appropriate `models/{wave_type}` directory, depending on type of waves the model uses.
3. if desired, rename your file to a name you would like to use later for the `input_model` variable in `mod_input`.
4. find the headers of the model file. if converting from an `.nc` file, the following commands can be used to access the headers from an interactive Python shell:
    1. `>>> import netCDF4`
    2. `>>> netCDF4.Dataset(f'./models/S/{model}.nc').variables.keys()`


## Input variables to define in `mod_input.py`: 
1. `input_model`: the name of the model file before the `.nc` extension.
2. `data_wave_type`: the polarity of the waves from which the model was made (i.e., `'S'` or `'P'`.).
3. `lat_header`: the latitude header in the input `{model}` file
4. `lon_header`: the longitude header in the input `{model}` file
5. `depth_header`: the depth header in the input `{model}` file
6. `property_header`: the model value header in the input `{model}` file
7. `voigt`: indication of whether or not to compute a voigt average (typically done for anisotropic models that give velocities for both horizontal and vertical components). if voigt averaging is desired, the user must indicate the headers of the values to average.
8. `convert_vel_to_perturb`: boolean value indicating whether or not the model values need to be converted from absolute velocities to velocity perturbations relative to PREM. The model update only handles perturbations, so if the model values are given in velocities, this conversion must be performed.
9. `convert_nc_to_csv`: boolean value for whether or not to convert a `{model}.nc` file to a `{model}.csv` file.
10. `zero_shells`: any depth shells to be registered as having model values of zero (i.e., for the case where crustal corrections have been applied and only the mantle will be updated. this is the intended case for SITRUS.)


## Necessary compute resources:
1. HPC cores: as many depth shells as there are in your model space/`shell_dimensions.csv`
2. Time allocation: Generally between several seconds and up to a day, depending on the size of the model file and the number of nodes in the input model.
3. Memory allocation: Can require up to tens of GB, if the input model file is very large.
