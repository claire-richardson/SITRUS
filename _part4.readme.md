## Purpose of script:
`part4.model_update` actually performs the iterative model update.


## Output/saved files:
1. if the model has not been updated yet, `models/{data_wave_type}/{model}_update/original_model_plot_files` is created and files for each depth shell are made that are plottable with the included SITRUS visualization tools are made.
2. a model update log is either created or appended to, with information about the update and its input parameters.
3. **in the main `models/{data_wave_type}/{model}_update/{job_id}_update` directory:**
    1. the final updated solution model: `{model}_final_updated_model_{job_id}.csv`
    2. the RMS profile of the final updated solution model: `{model}_final_updated_rms_{job_id}.csv`
    3. a file with predicted travel times, residuals, and misfit of the starting model: `{model}_starting_residuals_{job_id}.csv`
    4. a file with predicted travel times, residuals, and misfit of the final model: `{model}_final_residuals_{job_id}.csv`
    5. a file with total number of residuals used and misfit mean, standard deviation, and variance for the starting model and each subsequent iteration: `{model}_variance_reduction_{job_id}.csv`
    6. n subdirectories for each main update layer/assemblage of depth shells, starting with shell _x_ and ending with shell _y_, included in the update: `/layer_{n}_shell_{x}_to_{y}`

4. **in each main update layer subdirectory:**
    1. n subdirectories for each iteration in the layer `iteration_{m}`
    2. n subdirectories `smoothing_block_information/shell_{n}` for each depth shell in the model space, each of which contains files for each model element in that shell. These files contain all of the information for each raypath segment used to compute that element's final perturbation. _These files are only saved for the final model in the layer_.
    3. the final updated layer model: `{model}_updated_model.csv`
    4. the final updated layer RMS profile: `{model}_updated_rms.csv`
    5. the final smoothing radii for each 3D block: `smoothing_radii.csv`
  
5. **in each main update layer/iteration subdirectory**:
    1. a subdirectory with plottable files for the iteration's solution model `/plot_files`
    2. a file with predicted travel times, residuals, and misfit for the iteration's solution model: `{model}_starting_residuals_{job_id}.csv`
    3. a file with runtime information for the iteration: `{model}_update_info.txt`
    4. the final updated iteration model: `{model}_updated_model_itr_{i}.csv`
    5. the final updated iteration RMS profile: `{model}_updated_rms_itr_{i}.csv`


## Before running:
1. Successfully run all previous scripts.
5. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`: 
1. `input_model`: the global, whole mantle tomography model you would like to update. the name should match the string before the `.csv` extension in the `/models/{data_wave_type}` directory.
2. `HPC_cores`: the number of HPC cores you will be requesting. this should be the same as the number of depth shells in your model space/defined in `shell_dimensions.csv`
3. `data_wave_type`: the polarity of the waves from which the input travel time residuals were measured and the input model was made (i.e., `'S'` or `'P'`.).
4. `azimuthal_sectors`: the total number of azimuthal sectors to consider azimuthal sampling (should be the same value as was used in `part4`.
5. `starting_RMS_model_to_use`: the model whose RMS profile you would like to use for RMS weighting. this should generally be the same as `input_model`.
6. `reference_lat`: the latitudinal extent of the resolution at which to make your plot files (should generally match the value used for generating `block_dimensions.csv`
7. `reference_lon`: the longitudinal extent of the resolution at which to make your plot files (should generally match the value used for generating `block_dimensions.csv`
8. `residual_header`: for an update with _n_ layers, a list of _n_ strings indicating the header of the desired residual to use when updating the layer.
9. `layer_top_shells`: for an update with _n_ layers, a list of _n_ integers indicating the Shell IDs of the top most depth shells for each main update layer (i.e., for a layer starting at the top of the mantle, this value would be the Shell ID that corresponds to that depth shell in the model space)
10. `layer_base_shells`: for an update with _n_ layers, a list of _n_ integers indicating the Shell IDs of the bottom most depth shells for each main update layer (i.e., for a layer ending at the base of the mantle, this value would be the Shell ID that corresponds to the last depth shell in the model space)
11. `remove_residual_mean`: for an update with _n_ layers, a list of _n_ booleans indicating whether or not to remove the mean of the residuals after correcting for the starting tomography model. This acts as a proxy for earthquake relocation in the forward modeling workflow.
12. `freeze_previous_layers`: a boolean that indicates whether or not to freeze or continue to update a previous layer, during the update of the next layer.
13. `update_phases`: for an update with _n_ layers, a list of _n_ lists of the phases to use to update that layer.
14. `type_of_phase_subselection`: for an update with _n_ layers, a list of _n_ strings, either `'proportion'` or `'number'`, to indicate how to index the amount of data desired for each of the individual phases
15. `subselection_of_phase_data_to_use`: for an update with _n_ layers, a list of _n_ lists of the subselection value (values from 0 to 1 if `type` is `'proportion'`; values from 0 to the total number of a given phase's residuals if `type` is `number`).
16. `iteration_to_stop_RMS_weighting`: for an update with _n_ layers, a list of _n_ integers which indicate the iteration in that layer at which to stop applying the RMS weighting (i.e., this would be `[2]` for an update with one layer, where the user would like to apply RMS weighting for only the first iteration).
17. `residual_limits`: for an update with _n_ layers, a list of _n_ lists. If limits on residual values are desired, a given iterations list will look like; `[True, lower_lim, upper_lim]`. If limits are not desired, the list is simply `[False]`.
18. `cutoff_type`: for an update with _n_ layers, a list of _n_ strings, either `'total iterations'` or `'reduction'`, indicating the type of cutoff criteria to use.
19. `cutoff`: for an update with _n_ layers, a list of _n_ integers (for `total_iterations`) or floats (for `reduction` to indicate when layer should stop iterating.
20. `smoothing_radii`: for an update with _n_ layers, a list of _n_ lists with floats that indicate the smoothing radii to test for that layer.
21. `total_required_paths`: for an update with _n_ layers, a list of _n_ integers indicating the total number of required paths for a smoothing radius to be accepted.
22. `total_required_azimuths`: for an update with _n_ layers, a list of _n_ integers indicating the total number of required represented azimuthal sectors for a smoothing radius to be accepted.
23. `gaussian_cutoff_weight`: for an update with _n_ layers, a list of _n_ floats indicating the desired weight value to be assigned at the smoothing radius of each block, and used to apply weights to paths that fall within that smoothing radius based on their distance from the center.
24. `azimuthal_weighting`: for an update with _n_ layers, a list of _n_ booleans indicating whether to apply the weight when computing the final perturbation for each model element.
26. `special_weights`: for an update with _n_ layers, a list of _n_ lists with headers indicating any special weights added to the `{phase}_master_data.csv` files to be applied during the layer update.
27. `dataset_description`: for an update with _n_ layers, a list of _n_ strings with a short description of the dataset used to be added to the model update log. **do not use commas**
28. `reference_model`: 1D reference model to use for computing velocity perturbations (default: PREM; currently no other models are supported).



## Necessary compute resources:
1. HPC cores: as many depth shells as there are in your model space/`shell_dimensions.csv`
2. Time allocation: Between several hours and up to a week or more, depending on the input parameters chosen.