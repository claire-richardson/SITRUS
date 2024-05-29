import multiprocessing as mp
import mod_refmodels
import pandas as pd
import numpy as np
import mod_pandas
import mod_input
import mod_geo
import random
import shutil
import math
import glob
import time
import sys
import csv
import os

cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places
sectors = mod_input.azimuthal_sectors
sector_extent = 180. / mod_input.azimuthal_sectors
sector_ids = list(range(1, sectors + 1))
model = mod_input.input_model
job_id = str(sys.argv[1])
cores_requested = mod_input.HPC_cores
if mod_input.data_wave_type == 'S':
    out_property_header = 'dVs_%'
    RMS_header = 'RMS_dVs'
elif mod_input.data_wave_type == 'P':
    out_property_header = 'dVp_%'
    RMS_header = 'RMS_dVp'

# define data frames of the input model (df_model) and the RMS property values (df_rms)
# df_model will update with each iteration
df_original_model = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/{model}_grid_registered.csv')
df_model = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/{model}_grid_registered.csv')
df_rms = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.starting_RMS_model_to_use}_update/{mod_input.starting_RMS_model_to_use}_{RMS_header}.csv')

# compute azimuthal sector extents and define a dataframe with those attributes and sector ids (df_sectors)
df_sectors = pd.DataFrame(columns = ['ID', 'MIN', 'MAX'])
df_sectors['ID'] = list(range(1, mod_input.azimuthal_sectors + 1))
df_sectors['MIN'] = (df_sectors['ID'] - 1) * sector_extent
df_sectors['MAX'] = df_sectors['ID'] * sector_extent

# make a function to output the sector ID of an input azimuth.
def azimuthal_sector(az):
    if az == 360.:
        az = 0.
    if az >= 180.:
        az -= 180.
    df_az_slice = df_sectors.loc[(df_sectors['MIN'] <= az) & (df_sectors['MAX'] > az)]
    return int(df_az_slice['ID'])

# Find the standard deviation of the Gaussian function that will be used for weighting, and then define the Gaussian as a function of radius:
def gaussian(total_R, x, cutoff_weight):
    gaussian_numerator = ((total_R - mod_input.peak_center) ** 2.)
    gaussian_denominator = np.log(cutoff_weight / mod_input.peak_value)
    gaussian_standard_deviation = np.sqrt((gaussian_numerator / gaussian_denominator) * -0.5)
    exp_n = ((x - mod_input.peak_center) ** 2.)
    exp_d = (2 * gaussian_standard_deviation ** 2.)
    g = mod_input.peak_value * np.exp(-(exp_n / exp_d))
    return g

def rms(property_list):
    n = len(property_list)
    x = 0.
    for p in property_list:
        x += (p**2.)
    return np.sqrt(x / n)

# define dataframes with block and shell ID information, as well as near neighbors.
df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
all_shells = df_shells['SHELL#']
all_blocks = df_blocks['BLOCK#']

# make plot ready files for model visualization
try:
    plot_dir = f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/original_model_plot_files'
    os.mkdir(plot_dir)
    for shell_to_plot in all_shells:
        og_perturbs_file = f'{plot_dir}/{model}_shell_{shell_to_plot}_original_perturbs_plot_ready_{mod_input.reference_lat}deg_lat_by_{mod_input.reference_lon}deg_lon.csv'
        df_og_shell_data = df_original_model.loc[df_original_model['SHELL#'] == shell_to_plot]
        df_original = pd.DataFrame(data = {'LON': lons})

        # make empty lists that will be filled with numpy arrays for each latitude band, for each model.\
        og_dvs = []
        # loop through all latitudes in the model space
        for lat in np.arange(mod_input.start_lat, mod_input.final_lat, mod_input.reference_lat):
            # make empty lists to fill in the perturbations for the current latitude band
            lat_og_dvs = []
            # loop through all longitudes in the model space
            for lon in np.arange(mod_input.start_lon, mod_input.final_lon, mod_input.reference_lon):
                # find the block that the current lat/lon pair falls in
                block = mod_pandas.find_block_id(lat, lon)
                original_perturb = float(df_og_shell_data.loc[df_og_shell_data['BLOCK#'] == block][out_property_header])
                lat_og_dvs.append(original_perturb)
            og_dvs.append(np.array(lat_og_dvs))
            df_original[f'{lat}'] = lat_og_dvs
        df_original.to_csv(og_perturbs_file, index = False)
        og_dvs = np.array(og_dvs)
except:
    pass

# initialize an empty dataframe for variance reduction
df_variance = pd.DataFrame(columns = ['layer', 'iteration', 'total_paths', 'misfit_mean', 'misfit_std', 'misfit_var'])


# define functions for multiprocessing/parallelization
def backmap_residual(lists_of_paths, list_idx):
    df_idx = 0
    block_info_cols = ['SHELL#', 'BLOCK#', 'PHASE', 'PATH_ID', 'PATH_LENGTH_KM', out_property_header, 'AZIMUTH', 'SECTOR']
    for special_weight in layer_special_weights:
        block_info_cols.append(special_weight)
    df_block_info = pd.DataFrame(columns = block_info_cols)
    df_residuals = pd.DataFrame(columns = ['PHASE', 'PATH_ID', 'REF_pred_time', 'ORIG_DT', 'STARTING_pred_time', 'STARTING_DT', 'MISFIT'])

    for list_of_paths in lists_of_paths:
        phase = str(list_of_paths[0].split('/')[2])
        phase_name = str(phase.split('_')[0])
        df_data = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_master_data.csv')
        
        for path in list_of_paths:
            backmap = True
            tot_time = 0.
            df_p = pd.read_csv(path)
            path_id = int(path.split('_')[-3])
            # check if there are any segments in df_p that have a length of zero and drop these from the dataframe
            df_p = df_p.loc[df_p['SEG_DIST_KM'] != 0.].reset_index(drop = True)
            path_max_depth = df_p['START_DEPTH'].max()

            # check if the current path is above the layer stripping depth. continue only if it is.
            if path_max_depth <= layer_stripping_depth and path_max_depth > layer_depth_min:

                # find the reference model (e.g., prem) total travel time for the path. this will just be the total time in df_p.
                ref_model_time = float(df_p['TIME'].iloc[-1])

                # grab the dataset residual for the current path
                dataset_residual = float(df_data.loc[df_data['PATH_ID'] == path_id][layer_residual_header])

                if layer_residual_limits[0] == True:
                    if layer_residual_limits[1] <= dataset_residual <= layer_residual_limits[2]:
                        backmap = True
                    else:
                        backmap = False
                
                elif layer_residual_limits[0] == False:
                    backmap = True
                    
                if backmap == True:
                    ## LOOP 1: compute the travel time of the current path through the starting model (df_model).
                    for segment in range(len(df_p)):
                        # for each path segment in the current working path, find its attributes
                        dist_km = float(df_p['SEG_DIST_KM'].iloc[segment])
                        shell_no = int(df_p['SEG_SHELL#'].iloc[segment])
                        block_no = int(df_p['SEG_BLOCK#'].iloc[segment])
                        mid_depth = float(df_p['MID_DEPTH'].iloc[segment])
                        prem_val = float(df_p[f'SEG_V{mod_input.data_wave_type}'].iloc[segment])

                        # df_model_slice is the single line in the model dataframe for the element that the current path segment is in.
                        # from df_model_slice, find the velocity perturbation of that element (model_element_perturb)
                        model_element_perturb = float(df_model.loc[(df_model['SHELL#'] == shell_no) & (df_model['BLOCK#'] == block_no)][out_property_header])

                        # find the velocity value for the path segment through the starting model.
                        seg_v = prem_val + ((model_element_perturb / 100.) * prem_val)
                        df_p.loc[segment, f'SEG_V{mod_input.data_wave_type}'] = seg_v

                        # use that velocity value and the length of the path (in km) to find the travel time for just this path segment
                        seg_time = dist_km / seg_v
                        df_p.loc[segment, 'SEG_TIME'] = seg_time

                        # add that path segment travel time to the cumulative `tot_time` variable
                        tot_time += seg_time
                        df_p.loc[segment, 'TIME'] = tot_time

                    # update df_p, which is now the path through the model
                    df_p['TIME_PERCENT'] = (df_p['TIME'] / df_p['TIME'].iloc[-1]) * 100.

                    # step 1 of flow chart:
                    # total starting model time - total reference model time (residual between starting model and reference model)
                    predicted_residual = tot_time - ref_model_time

                    # step 2 of flow chart:
                    # difference between Hongyu's residual and the predicted residual from the starting model
                    unexplained_difference = dataset_residual - predicted_residual

                    # path id, PREM, total time through the starting model, predicted resid for starting model, residual from our dataset, difference between the two
                    df_residuals.loc[len(df_residuals)] = [phase, path_id, ref_model_time, dataset_residual, tot_time, predicted_residual, unexplained_difference]

                    # step 3 of flow chart (backmap unexplained difference into paths through starting model):
                    # Check where we are in the layer stripping situation
                    # first, if we DONT want to freeze the previous layer: 
                    if freeze_previous_layer == False or freeze_previous_layer == None:
                        df_p_update = df_p.copy()
                        df_layer_rms = df_rms.copy()
                        # x is a multiplyer for mapping the unexplained difference above into velocity perturbations along the raypath in the next loop
                        x = tot_time / (tot_time + unexplained_difference)
                        df_p_update['norm_rms'] = 0.

                    elif freeze_previous_layer == True:
                        df_p_update = df_p.loc[df_p['SEG_SHELL#'] >= layer_top_shell].copy().reset_index(drop = True)
                        df_layer_rms = df_rms.loc[df_rms['SHELL#'] >= layer_top_shell].copy().reset_index(drop = True)
                        
                        if df_layer_rms[f'NORM_{RMS_header}'].all() == 0.:
                            df_layer_rms[f'NORM_{RMS_header}'] = 1.
                        else:
                            df_layer_rms[f'NORM_{RMS_header}'] = df_layer_rms[f'NORM_{RMS_header}'] / df_layer_rms[f'NORM_{RMS_header}'].max()
                        
                        if math.isnan(df_p_update['SEG_TIME'].sum()) == True or df_p_update['SEG_TIME'].sum() == 0.:
                            with open(f'{iteration_directory_path}/{list_idx}_nan_bug.txt', 'a') as fout:
                                fout.write(f"BUG WITH PROCESS {list_idx}! total time sum: {df_p_update['SEG_TIME'].sum()}\n")
                            df_p_update.to_csv(f'{iteration_directory_path}/idx_{list_idx}_{phase}_{path_id}_bugged_path.csv', index = False)

                        x = df_p_update['SEG_TIME'].sum() / (df_p_update['SEG_TIME'].sum() + unexplained_difference)
                            
                        df_p_update['norm_rms'] = 0.

                    for segment in range(len(df_p_update)):
                        segment_shell_no = df_p_update['SEG_SHELL#'].iloc[segment]
                        norm_rms_update = float(df_layer_rms.loc[df_layer_rms['SHELL#'] == segment_shell_no][f'NORM_{RMS_header}'])
                        df_p_update.loc[segment, 'norm_rms'] = norm_rms_update                    


                    # v_update ensures that velocity perturbations are updated by equal partitioning of velocity percent, rather than time percent
                    df_p_update['v_update'] = df_p_update[f'SEG_V{mod_input.data_wave_type}'] * x

                    # compute the corresponding travel time for the segment
                    df_p_update['t_update'] = df_p_update['SEG_DIST_KM'] / df_p_update['v_update']

                    # Find the difference between the updated time and the original model + reference model time:
                    df_p_update['delta_t'] = df_p_update['t_update'] - df_p_update['SEG_TIME']

                    # Weight delta_t according to the RMS values:
                    df_p_update['delta_t_rms'] = df_p_update['delta_t'] * df_p_update['norm_rms']

                    total_weighted_delta_t = df_p_update['delta_t_rms'].sum()

                    if total_weighted_delta_t == 0.:
                        df_p_update['final_time_portions'] = 0.
                    else:
                        df_p_update['final_time_portions'] = (df_p_update['delta_t_rms'] / total_weighted_delta_t) * unexplained_difference
                    
                    df_p_update['SEG_TIME'] = df_p_update['SEG_TIME'] + df_p_update['final_time_portions']
                    df_p_update[f'SEG_V{mod_input.data_wave_type}'] = df_p_update['SEG_DIST_KM'] / df_p_update['SEG_TIME']

                    t = 0.
                    for val in range(len(df_p_update['SEG_TIME'])):
                        t += df_p_update['SEG_TIME'].iloc[val]
                        df_p_update.loc[val, 'TIME'] = t

                    df_p_update['TIME_PERCENT'] = (df_p_update['TIME'] / df_p_update['TIME'].iloc[-1]) * 100.


                    # Isolate the unique block and shell IDS in each path and condense to representative block values
                    df_path_elements = df_p_update.groupby(['SEG_SHELL#','SEG_BLOCK#']).size().reset_index().rename(columns = {0: 'TOTAL_PATHS'})

                    # loop through each unique combination
                    for pair in range(len(df_path_elements)):
                        shell_no = df_path_elements['SEG_SHELL#'].iloc[pair]
                        block_no = df_path_elements['SEG_BLOCK#'].iloc[pair]

                        # find all of the segments in the current block
                        df_slice = df_p_update.loc[(df_p_update['SEG_SHELL#'] == shell_no) & (df_p_update['SEG_BLOCK#'] == block_no)]

                        updated_velocity = df_slice[f'SEG_V{mod_input.data_wave_type}'].iloc[0]
                        mid_depth = df_slice['MID_DEPTH'].iloc[0]
                        prem_velocity = mod_refmodels.prem_vel(mod_input.data_wave_type, mid_depth)
                        percent_perturbation = (100. * (updated_velocity / prem_velocity)) - 100.

                        # if we're in the first iteration, then add all new data. if not, just update the velocity perturbation column.
                        total_length_km = df_slice['SEG_DIST_KM'].sum()
                        azimuth = df_slice['AZIMUTH'].mean()
                        az_sector = azimuthal_sector(azimuth)

                        # add all of that information to the last row of the block file
                        new_block_info = [shell_no, block_no, phase, path_id, total_length_km, percent_perturbation, azimuth, az_sector]

                        for special_weight in layer_special_weights:
                            new_block_info.append(float(df_data.loc[df_data['PATH_ID'] == path_id][special_weight]))
                        df_block_info.loc[len(df_block_info)] = new_block_info

                    df_p_update.to_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.backmapped_paths_directory}_{job_id}/{phase_name}_{path_id}_backmapped_segments.csv', index = False)

                    if len(df_block_info) > 10000 or path == lists_of_paths[-1][-1]:
                        df_block_info.to_csv(f'{layer_directory}/tmp_block_info_idx_{list_idx}/df_index_{df_idx}.csv', index = False)
                        df_idx += 1
                        df_block_info = df_block_info.drop(df_block_info.index)

    df_residuals.to_csv(f'{iteration_directory_path}/tmp_residuals_{list_idx}.csv', index = False)
    
    # here, i need to merge these medium sized files i just made into shell specific files. {layer_directory}/tmp_block_info_idx_{list_idx}/shell_{shell#}.csv
    for d_id in range(df_idx):
        df_block_info = pd.read_csv(f'{layer_directory}/tmp_block_info_idx_{list_idx}/df_index_{d_id}.csv')
        for shell_no in shells_to_update:
            df_block_info_shell = df_block_info.loc[df_block_info['SHELL#'] == shell_no].drop(columns = ['SHELL#'])
            try:
                df_shell = pd.read_csv(f'{layer_directory}/tmp_block_info_idx_{list_idx}/shell_{shell_no}.csv')
            except:
                df_shell = pd.DataFrame(columns = df_block_info_shell.columns)
            df_shell = pd.concat([df_shell, df_block_info_shell]).reset_index(drop = True)
            df_shell.to_csv(f'{layer_directory}/tmp_block_info_idx_{list_idx}/shell_{shell_no}.csv', index = False)
        os.remove(f'{layer_directory}/tmp_block_info_idx_{list_idx}/df_index_{d_id}.csv')
        

def model_smoothing(shell_no):
    df_updated_model_shell = pd.DataFrame(data = {'SHELL#': shell_no, 'BLOCK#': all_blocks, out_property_header: 0.})    
    for block_no in all_blocks:
        df_smoothing = pd.read_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv')

        # ONLY ENTER THE SMOOTHING LOOP if there is at least one segment in df_smoothing.
        # otherwise, if there are 0 paths, don't bother smoothing and just use the value of the current block.
        if len(df_smoothing) > 0:
            # make a histogram of all of the sectors in df_smoothing
            unique_path_segment_sectors_hist = np.histogram(df_smoothing['SECTOR'], bins = np.arange(1, mod_input.azimuthal_sectors + 2, 1))[0]
    
            # azimuthal smoothing is an option in mod_input. if it's on, also include azimuthal information in smoothing.
            # if not, just include gaussian weights based on distance to the center of each path segment in df_smoothing.
            if layer_azimuthal_weighting == True:
                # compute the main smoothed perturbation for the center block!
                # this method incorporates a weight based on the azimuthal coverage of the radius
                azimuthal_denominator = 0.
                azimuthal_sector_averages = []
                mean_element_perturbation = 0.

                # loop through each of the six values in unique_path_segment_sectors_hist
                for path_sectors in range(len(unique_path_segment_sectors_hist)):

                    # for each of those values, determine whether or not the value is 0
                    paths_in_sector = unique_path_segment_sectors_hist[path_sectors]

                    # if the value is not 0 (i.e., if there are 1 or more paths in the radius within the current sector) then compute a mean
                    if paths_in_sector > 0:
                        azimuthal_denominator += 1.

                        # slice df_smoothing to have just the path segments in the current path sector
                        df_final_smoothing = df_smoothing.loc[df_smoothing['SECTOR'] == (path_sectors + 1)]
                        pre_az_numerator = 0.
                        pre_az_denominator = 0.

                        for segment_final in range(len(df_final_smoothing)):
                            weights = []
                            segment_final_pert = df_final_smoothing[out_property_header].iloc[segment_final]
                            weights.append(df_final_smoothing['GAUSS_WEIGHT'].iloc[segment_final])
                            if layer_path_length_weighting == True:
                                weights.append(df_final_smoothing['PATH_LENGTH_KM'].iloc[segment_final])
                                
                            for special_weight in layer_special_weights:
                                weight = df_final_smoothing[special_weight].iloc[segment_final]
                                
                                if math.isnan(weight) == False and weight != 0.:
                                    weights.append(weight)

                            if math.isnan(segment_final_pert) == False:
                                pre_az_numerator += (segment_final_pert * np.prod(weight))
                                pre_az_denominator += (np.prod(weights))

                        # compute the final weighted average for the current sector, considering path segmenth length,...
                        # gaussian weight, comprehensive weight, and azimuthal coverage.
                        if math.isnan(pre_az_numerator / pre_az_denominator) == False:
                            azimuthal_sector_averages.append(pre_az_numerator / pre_az_denominator)
                    else:
                        pass

                # now that we have all of the weighted averages for each sector, combine them into one FINAL, TOTAL weighted average for the current block.
                for az_sector_average in azimuthal_sector_averages:
                    mean_element_perturbation += (1. / azimuthal_denominator) * az_sector_average

            # if the azimuthal weighting option is turned off, just compute a weighted average based on the weights indicated in mod_input.
            elif layer_azimuthal_weighting == False:
                smoothing_numerator = 0.
                smoothing_denominator = 0.

                for segment_final in range(len(df_smoothing)):
                    weights = []
                    segment_final_pert = df_smoothing[out_property_header].iloc[segment_final]
                    weights.append(df_smoothing['GAUSS_WEIGHT'].iloc[segment_final])
                    if layer_path_length_weighting == True:
                        weights.append(df_smoothing['PATH_LENGTH_KM'].iloc[segment_final])
                        
                    for special_weight in layer_special_weights:
                        weight = df_smoothing[special_weight].iloc[segment_final]
                        
                        if math.isnan(weight) == False and weight != 0.:
                            weights.append(weight)

                    if math.isnan(segment_final_pert) == False:
                        smoothing_numerator += (segment_final_pert * np.prod(weights))
                        smoothing_denominator += (np.prod(weights))

                if smoothing_denominator != 0.:
                    mean_element_perturbation = smoothing_numerator / smoothing_denominator
        else:
            mean_element_perturbation = float(df_model.loc[(df_model['SHELL#'] == shell_no) & (df_model['BLOCK#'] == block_no)][out_property_header])

        element_idx = df_updated_model_shell.loc[(df_updated_model_shell['SHELL#'] == shell_no) & (df_updated_model_shell['BLOCK#'] == block_no)].index
        df_updated_model_shell.loc[element_idx, out_property_header] = mean_element_perturbation

    df_updated_model_shell.to_csv(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv', index = False)
    
    
def merge_block_files(shell_no):
    df_main_info_file = pd.read_csv(f'{block_centric_path}/shell_{shell_no}.csv')
    df_main_info_file = df_main_info_file.drop(df_main_info_file.index)
    
    for list_index in range(len(paths_lists)):
        if os.path.exists(f'{layer_directory}/tmp_block_info_idx_{list_index}/shell_{shell_no}.csv'):
            df_tmp_info_file = pd.read_csv(f'{layer_directory}/tmp_block_info_idx_{list_index}/shell_{shell_no}.csv')
            df_main_info_file = pd.concat([df_main_info_file, df_tmp_info_file])
            os.remove(f'{layer_directory}/tmp_block_info_idx_{list_index}/shell_{shell_no}.csv')

    df_main_info_file = df_main_info_file.reset_index(drop = True)
    df_main_info_file.to_csv(f'{block_centric_path}/shell_{shell_no}.csv', index = False)


def update_smoothing_block(shell_no):
    df_updated_shell = pd.read_csv(f'{block_centric_path}/shell_{shell_no}.csv')
    df_smoothing_radii = pd.read_csv(f'{layer_directory}/shell_{shell_no}_smoothing_radii.csv')
    for block_no in all_blocks:
        block_center_lat = float(df_blocks.loc[df_blocks['BLOCK#'] == block_no]['CENTER_LAT'])
        block_center_lon = float(df_blocks.loc[df_blocks['BLOCK#'] == block_no]['CENTER_LON'])

        block_smoothing_radius = float(df_smoothing_radii.loc[df_smoothing_radii['BLOCK#'] == block_no]['RADIUS'])
        df_smoothing_block = pd.read_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv')
        neighbors = df_smoothing_block['NEIGHBOR'].unique()
        
        df_smoothing_block = df_smoothing_block.drop(df_smoothing_block.index)
        
        for neighbor in neighbors:
            df_updated_block = df_updated_shell.loc[df_updated_shell['BLOCK#'] == neighbor].rename(columns = {'BLOCK#': 'NEIGHBOR'})
            df_smoothing_block = pd.concat([df_smoothing_block, df_updated_block])
        df_smoothing_block = df_smoothing_block.reset_index(drop = True)
        df_smoothing_block['GAUSS_WEIGHT'] = 0.
            
        for segment_id in range(len(df_smoothing_block)):
            segment_center_lat = df_smoothing_block['CENTER_LAT'].iloc[segment_id]
            segment_center_lon = df_smoothing_block['CENTER_LON'].iloc[segment_id]
            
            dist_to_center = mod_geo.GCP_length(block_center_lat, block_center_lon, segment_center_lat, segment_center_lon)
            gaus_weight = gaussian(block_smoothing_radius, dist_to_center, layer_cutoff_weight)
            df_smoothing_block.loc[segment_id, 'GAUSS_WEIGHT'] = gaus_weight
        
        df_smoothing_block.to_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv', index = False)


def find_smoothing_radii(shell_no):
    df_shell_smoothing_radii = pd.DataFrame(data = {'BLOCK#': list(all_blocks), 'RADIUS': 0.})
    df_shell_elements = pd.read_csv(f'{block_centric_path}/shell_{shell_no}.csv')
    
    for block_no in all_blocks:
        # df_element has the information for the current element of all of the path segments that have gone through that element for all of the update phases
        element_center_lat = float(df_blocks.loc[df_blocks['BLOCK#'] == block_no]['CENTER_LAT'])
        element_center_lon = float(df_blocks.loc[df_blocks['BLOCK#'] == block_no]['CENTER_LON'])
        df_neighbors = pd.read_csv(f'./{mod_input.near_neighbors_directory}/block_{block_no}_neighbors.csv')

        for radius in layer_smoothing_radii:
            # find block centers within the current radius, INCLUDING the current element
            # grab these from the near_neighbors.csv file.
            # df_near_neighbors is a slice of the neighbors.csv file that has only the near neighbors within the current working radius
            df_near_neighbors = df_neighbors.loc[df_neighbors['RADIUS_DEG'] <= radius]
            eligible_blocks = df_near_neighbors['NEIGHBOR']
            
            df_smoothing = pd.DataFrame()
            for neighbor in eligible_blocks:
                df_neighbor = df_shell_elements.loc[df_shell_elements['BLOCK#'] == neighbor].rename(columns = {'BLOCK#': 'NEIGHBOR'})
                df_smoothing = pd.concat([df_smoothing, df_neighbor])
            df_smoothing = df_smoothing.reset_index(drop = True)

            # df_unique_paths is a dataframe that has the total number of segments for each full path within the current radius.
            df_unique_paths = df_smoothing.groupby(['PHASE','PATH_ID']).size().reset_index().rename(columns = {0: 'TOTAL_SEGMENTS'})
            df_unique_paths['SECTOR'] = 0
            
            # test to see if all of those blocks meet the criteria we define
            # first, check if there are enough paths in general in the radius
            total_paths_in_radius = len(df_unique_paths)

            # if there aren't, check the next radius
            if (total_paths_in_radius < layer_total_required_paths) and (radius < layer_max_radius):
                pass

            elif (total_paths_in_radius < layer_total_required_paths) and (radius == layer_max_radius):
                break

            # if there are, then check the azimuthal coverage requirements:
            else:
                for unique_path_idx in range(len(df_unique_paths)):
                    unique_path_phase = df_unique_paths['PHASE'].iloc[unique_path_idx]
                    unique_path_id = df_unique_paths['PATH_ID'].iloc[unique_path_idx]
                    unique_path_azimuth = df_smoothing.loc[(df_smoothing['PHASE'] == unique_path_phase) & (df_smoothing['PATH_ID'] == unique_path_id)]['AZIMUTH'].mean()
                    unique_path_sector = azimuthal_sector(unique_path_azimuth)
                    df_unique_paths.loc[unique_path_idx, 'SECTOR'] = unique_path_sector

                unique_path_sectors_hist = np.histogram(df_unique_paths['SECTOR'], bins = np.arange(1, mod_input.azimuthal_sectors + 2, 1))[0]

                # count how many sectors are covered
                total_sectors_covered = 0
                for count in unique_path_sectors_hist:
                    if count > 0:
                        total_sectors_covered += 1

                # if the total number of sectors covered meets or exceeds the set threshold, OR the radius is at its maximum allowed, then set sufficient coverage to True
                if (total_sectors_covered >= layer_total_required_azimuths) or (radius == layer_max_radius):
                    break

                # otherwise, coverage is still insufficient
                else:
                    pass

        # now that the smoothing radius and the paths/path segments within that smoothing radius for the current block has been determined, find the Gaussian weights for each of those segments
        if len(df_smoothing) > 0:
            # find the Gaussian weight for all segments in the radius to be smoothed (i.e., all segments in df_smoothing)
            df_smoothing['GAUSS_WEIGHT'] = 0.
            for unique_path_segment in range(len(df_smoothing)):

                unique_path_segment_center_lat = df_smoothing['CENTER_LAT'].iloc[unique_path_segment]
                unique_path_segment_center_lon = df_smoothing['CENTER_LON'].iloc[unique_path_segment]

                # find the distance from the center of the current element to the center of the path segment
                dist_to_center = mod_geo.GCP_length(element_center_lat, element_center_lon, unique_path_segment_center_lat, unique_path_segment_center_lon)
                
                # compute the gaussian weight for that path segment and add it to df_smoothing
                gaus_weight = gaussian(radius, dist_to_center, layer_cutoff_weight)
                df_smoothing.loc[unique_path_segment, 'GAUSS_WEIGHT'] = gaus_weight
                
        df_smoothing.to_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv', index = False)

        smoothing_idx = df_shell_smoothing_radii.loc[df_shell_smoothing_radii['BLOCK#'] == block_no].index
        df_shell_smoothing_radii.loc[smoothing_idx, 'RADIUS'] = radius
    df_shell_smoothing_radii.to_csv(f'{layer_directory}/shell_{shell_no}_smoothing_radii.csv', index = False)



########################
# START MODEL UPDATE!! #
########################
start = time.time()

update_path = f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/{job_id}_update'
block_centric_path = f'{update_path}/{mod_input.block_centric_directory}'

try:
    os.mkdir(update_path)
    os.mkdir(block_centric_path)
except:
    pass

layers_complete = 0
# Start looping through each main layer. Each will be iterated through until they reach sufficient variance reduction.
for layer_bottom_shell in mod_input.layer_base_shells:
    layer_start = time.time()
    layer_top_shell = mod_input.layer_top_shells[layers_complete]
    layer_stripping_depth = mod_pandas.get_shell_info(layer_bottom_shell)[3]
    layer_depth_min = mod_pandas.get_shell_info(layer_top_shell)[1]
    freeze_previous_layer = mod_input.freeze_previous_layers[layers_complete]
    phases_to_update = mod_input.update_phases[layers_complete]
    type_of_data_subselection = mod_input.type_of_phase_subselection[layers_complete]
    subselection_of_phases_to_update = mod_input.subselection_of_phase_data_to_use[layers_complete]
    iteration_to_stop_RMS = mod_input.iteration_to_stop_RMS_weighting[layers_complete]
    shells_to_update = list(df_shells.loc[(df_shells['SHELL#'] >= layer_top_shell) & (df_shells['SHELL#'] <= layer_bottom_shell)]['SHELL#'])
    layer_residual_limits = mod_input.residual_limits[layers_complete]
    layer_variance_cutoff_type = mod_input.variance_cutoff_type[layers_complete]
    layer_variance_cutoff = mod_input.variance_cutoff[layers_complete]
    layer_smoothing_radii = mod_input.smoothing_radii[layers_complete]
    layer_max_radius = layer_smoothing_radii[-1]
    layer_total_required_paths = mod_input.total_required_paths[layers_complete]
    layer_total_required_azimuths = mod_input.total_required_azimuths[layers_complete]
    layer_cutoff_weight = mod_input.cutoff_weight[layers_complete]
    layer_azimuthal_weighting = mod_input.azimuthal_weighting[layers_complete]
    layer_path_length_weighting = mod_input.path_length_weighting[layers_complete]
    layer_special_weights = mod_input.special_weights[layers_complete]
    layer_residual_header = mod_input.residual_header[layers_complete]
    layer_dataset_description = mod_input.dataset_description[layers_complete]
    
    layers_complete += 1
    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
        fout.write(f'-----   STARTING UPDATE FOR LAYER {layers_complete} of {len(mod_input.layer_base_shells)}, DEPTH RANGE: {layer_depth_min} km - {layer_stripping_depth} km   -----\n')
    
    # write to model update log
    if os.path.exists(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv'):
    with open(f'./test_log.csv', 'a') as fout:
        fout.write(f'{job_id},')
        fout.write(f'{layers_complete} of {len(mod_input.layer_base_shells)},')
        fout.write(f'{model},')
        fout.write(f'{mod_input.starting_RMS_model_to_use},')
        fout.write(f'{mod_input.reference_model},')
        fout.write(f'{layer_dataset_description},')
        fout.write(f'{layer_residual_header},')
        fout.write(f"[{';'.join(str(i) for i in phases_to_update)}],")
        fout.write(f'{type_of_data_subselection},')
        fout.write(f"[{';'.join(str(i) for i in subselection_of_phases_to_update)}],")
        fout.write(f"[{';'.join(str(i) for i in layer_residual_limits)}],")
        fout.write(f'{layer_depth_min} km - {layer_stripping_depth} km,')
        fout.write(f'{iteration_to_stop_RMS},')
        fout.write(f'{layer_variance_cutoff_type},')
        fout.write(f'{layer_variance_cutoff},')
        fout.write(f'{mod_input.azimuthal_sectors},')
        fout.write(f'{layer_azimuthal_weighting},')
        fout.write(f'{layer_path_length_weighting},')
        fout.write(f"[{';'.join(str(i) for i in layer_special_weights)}],")
        fout.write(f"[{';'.join(str(i) for i in layer_smoothing_radii)}],")
        fout.write(f'{layer_total_required_paths},')
        fout.write(f'{layer_total_required_azimuths},')
        fout.write(f'{layer_cutoff_weight},')

    else:
        with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv', 'w') as fout:
            fout.write('PID,Layer,Input_model,Input_RMS,Reference_model,Dataset_description,Residual_header,Phases_used,Type_of_phase_subselection,Phase_subselection,Residual_limits,Layer_dimensions,Stopped_RMS_weighting_iteration,Variance_cutoff_type,Variance_cutoff,Azimuthal_sectors,Azimuthal_weighting,Path_length_weighting,Special_weights,Smoothing_radii,Total_required_paths,Total_required_sectors,Gaussian_cutoff_weight\n')
        with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv', 'a') as fout:
            fout.write(f'{job_id},')
            fout.write(f'{layers_complete} of {len(mod_input.layer_base_shells)},')
            fout.write(f'{model},')
            fout.write(f'{mod_input.starting_RMS_model_to_use},')
            fout.write(f'{mod_input.reference_model},')
            fout.write(f'{layer_dataset_description},')
            fout.write(f'{layer_residual_header},')
            fout.write(f"[{';'.join(str(i) for i in phases_to_update)}],")
            fout.write(f'{type_of_data_subselection},')
            fout.write(f"[{';'.join(str(i) for i in subselection_of_phases_to_update)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_residual_limits)}],")
            fout.write(f'{layer_depth_min} km - {layer_stripping_depth} km,')
            fout.write(f'{iteration_to_stop_RMS},')
            fout.write(f'{layer_variance_cutoff_type},')
            fout.write(f'{layer_variance_cutoff},')
            fout.write(f'{mod_input.azimuthal_sectors},')
            fout.write(f'{layer_azimuthal_weighting},')
            fout.write(f'{layer_path_length_weighting},')
            fout.write(f"[{';'.join(str(i) for i in layer_special_weights)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_smoothing_radii)}],")
            fout.write(f'{layer_total_required_paths},')
            fout.write(f'{layer_total_required_azimuths},')
            fout.write(f'{layer_cutoff_weight},')

    layer_directory = f'{update_path}/layer_{layers_complete}_shell_{layer_top_shell}_to_{layer_bottom_shell}'
    try:
        os.mkdir(layer_directory)
    except:
        pass

    for shell_no in all_shells:
        with open(f'{block_centric_path}/shell_{shell_no}.csv', 'w') as file:
            file.write(f'BLOCK#,PHASE,PATH_ID,PATH_LENGTH_KM,{out_property_header},AZIMUTH,SECTOR\n')
    
    # loop through all of the paths for all of the phases to be updated in this layer and split them up into smaller lists for multiprocessing.
    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
        fout.write(f'- Partitioning paths into batches for backmapping\n')
    paths_lists = [] # a list of length `cores_requested` with phase-wise sublists to be processed on each CPU core
    paths_list = [] # a single list of ALL paths for ALL phases to be updated.

    for phase_idx in range(len(phases_to_update)):
        phase = phases_to_update[phase_idx]
        phase_subselection = subselection_of_phases_to_update[phase_idx]
        phase_name = phase.split('_')[0]
        all_paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.paths_directory}/{phase_name}_*.csv')
        if type_of_data_subselection == 'proportion':
            subset_of_paths = int(len(all_paths) * phase_subselection)
            paths = random.sample(all_paths, subset_of_paths)
        elif type_of_data_subselection == 'number':
            paths = random.sample(all_paths, phase_subselection)

        # append all paths to the main `paths_list`
        for path in paths:
            paths_list.append(path)
    
    # fill `paths_lists` with individual, semi-evenly divided lists, that each contain further divided lists determined by phase. like: [[[vs_phase_paths (n = 7000)],[vs_2_phase_paths; n = 3000]], [etc]]
    tot_paths = len(paths_list)
    paths_per_core = math.floor(tot_paths / cores_requested)

    # determine the total approximate number of paths that should go in each main list
    start_path_idx = 0
    for core in range(1, cores_requested + 1):
        if core != cores_requested:
            path_list = paths_list[start_path_idx:start_path_idx + paths_per_core]
        elif core == cores_requested:
            path_list = paths_list[start_path_idx:]

        # in the sublist, make further sublists based on phase
        # first, determine which phases are present in this sublist
        path_phases = []
        for path in path_list:
            path_phase = path.split('/')[2]
            
            if path_phase not in path_phases:
                path_phases.append(path_phase)

        # then divide the sublist based on the phases that are present.
        paths_sublist = []
        for phase in path_phases:
            phase_sublist = []
            
            for path in path_list:
                path_phase = path.split('/')[2]
            
                if path_phase == phase:
                    phase_sublist.append(path)
                    
            paths_sublist.append(phase_sublist)
    
        paths_lists.append(paths_sublist)
        start_path_idx += paths_per_core
        
    # make a temporary directory for all of the different jobs that will run for backmapping and path>block centric conversion. these will be merged once the process is complete
    for list_index in range(len(paths_lists)):
        try:
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Making temporary block info files: batch {list_index + 1} of {len(paths_lists)}\n')
            os.mkdir(f'{layer_directory}/tmp_block_info_idx_{list_index}')
        except:
            pass
    
    unexplained_means = []
    unexplained_stds = []
    unexplained_vars = []
    total_paths_itr = []

    ############################################
    # START WHILE LOOP FOR CURRENT MODEL LAYER #
    ############################################
    layer_iteration = 0
    continue_layer_iterating = True
    while continue_layer_iterating == True:
        iteration_start = time.time()
        layer_iteration += 1

        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'\n- ** Begin iteration #{layer_iteration} **\n')

        # make the directory for the current layer and iteration of the model in the main model update directory if it's not already made.
        iteration_directory_path = f'{layer_directory}/iteration_{layer_iteration}'
        try:
            os.mkdir(f'{iteration_directory_path}')
        except:
            pass

        for phase in phases_to_update:
            # make a directory for the forthcoming backmapped paths if it doesn't already exist.
            try:
                os.mkdir(f'./{mod_input.phases_directory}/{phase}/{mod_input.backmapped_paths_directory}_{job_id}')
            except:
                pass

        if iteration_to_stop_RMS != None and iteration_to_stop_RMS <= layer_iteration:
            df_rms[RMS_header] = 1.
            df_rms[f'NORM_{RMS_header}'] = 1.

        ############################################
        # PART 1: BACKMAP RESIDUALS INTO RAYPATHS. #
        ############################################
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start backmapping and converting from path-to-block format\n')
        backmapping_start = time.time()
        
        if __name__ == '__main__':
            process_list = []
            process_idx = 0
            for list_index in range(len(paths_lists)):
                process_idx += 1
                with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                    fout.write(f'  - Starting process {process_idx} of {len(paths_lists)}')
                p = mp.Process(target = backmap_residual, args = (paths_lists[list_index], list_index,))
                p.start()
                process_list.append(p)

            for process in process_list:
                process.join()
        
        backmapping_time = time.time() - backmapping_start

        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished backmapping and converting from path-to-block format; runtime: {backmapping_time / 60} minutes / {(backmapping_time / 60) / 60} hours\n')

        
        ## Merge all files for the different processes into the main block files:
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start merging temporary block files into main block files\n')
        merge_start = time.time()
        
        if __name__ == '__main__':
            merge_process_list = []

            for shell_no in shells_to_update:
                with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                    fout.write(f'  - Starting process for shell {shell_no}; total of {len(shells_to_update)} shells to merge\n')
                p = mp.Process(target = merge_block_files, args = (shell_no,))
                p.start()
                merge_process_list.append(p)
                
            for merge_process in merge_process_list:
                merge_process.join()

        merge_time = time.time() - merge_start
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished merging temporary block files into main block files; runtime: {merge_time} seconds / {merge_time / 60} minutes\n')
        
        
        #### NOW, ALL OF THE BACKMAPPED PATHS ARE SAVED ####
        #### LOOP THROUGH ALL OF THE BACKMAPPED PATHS AND FIND THE PERTURBATION INFORMATION FOR ALL BLOCKS THAT EACH PATH GOES THROUGH ####

        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start computing variance\n')
        residual_start = time.time()
        df_residuals = pd.DataFrame()
        for list_index in range(len(paths_lists)):
            df_tmp_residuals = pd.read_csv(f'{iteration_directory_path}/tmp_residuals_{list_index}.csv')
            df_residuals = pd.concat([df_residuals, df_tmp_residuals])
            os.remove(f'{iteration_directory_path}/tmp_residuals_{list_index}.csv')

        df_residuals = df_residuals.reset_index(drop = True)
        df_residuals.to_csv(f'{iteration_directory_path}/residuals.csv', index = False)
        if layer_iteration == 1 and layers_complete == 1:
            df_residuals.to_csv(f'{update_path}/{model}_starting_residuals_{job_id}.csv', index = False)
        unexplained_diffs = df_residuals['MISFIT'].to_numpy()
        
        unexplained_means.append(np.mean(unexplained_diffs))
        unexplained_stds.append(np.std(unexplained_diffs))
        unexplained_vars.append(np.var(unexplained_diffs))
        total_paths_itr.append(len(df_residuals))
        
        if layer_variance_cutoff_type == 'total iterations' and layer_variance_cutoff == 1: 
            continue_layer_iterating = False
        
        if layer_iteration > 1:
            std1 = unexplained_stds[-1]
            std2 = unexplained_stds[-2]
            std_diff = std2 - std1

            var1 = unexplained_vars[-1]
            var2 = unexplained_vars[-2]
            var_diff = var2 - var1

            if layer_variance_cutoff_type == 'reduction':
                if var_diff < layer_variance_cutoff:
                    continue_layer_iterating = False

            elif layer_variance_cutoff_type == 'total iterations':
                if layer_iteration == layer_variance_cutoff:
                    continue_layer_iterating = False


            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Iteration #{layer_iteration}; variance: {var1}, variance difference: {var_diff}\n')
        residual_time = time.time() - residual_start

        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished computing variance: {residual_time} seconds\n')
        
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- * backmapping COMPLETE for all phases\n')
            fout.write(f'- * beginning smoothing...\n')

        # if this is the first iteration for the current layer, find the appropriate smoothing radius for each block based on the new block information.
        if layer_iteration == 1:
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Start computation of smoothing radii (first iteration only)\n')
            smoothing_radii_start = time.time()
            
            try:
                os.mkdir(f'{layer_directory}/smoothing_block_information')
            except:
                pass
            
            if __name__ == '__main__':
                find_radii_process_list = []
                for shell_no in shells_to_update:
                    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                        fout.write(f'  - Starting process for shell {shell_no}; total: {len(shells_to_update)} shells to update\n')
                    try:
                        os.mkdir(f'{layer_directory}/smoothing_block_information/shell_{shell_no}')
                    except:
                        pass
                    p = mp.Process(target = find_smoothing_radii, args = (shell_no,))
                    p.start()
                    find_radii_process_list.append(p)
              
                for find_radii_process in find_radii_process_list:
                    find_radii_process.join()
            
            smoothing_radii_time = time.time() - smoothing_radii_start

            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Finished finding smoothing radii; runtime: {smoothing_radii_time / 60} minutes / {(smoothing_radii_time / 60) / 60} hours\n')

        else:
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Start updating smoothing blocks\n')
            smoothing_radii_start = time.time()
            
            if __name__ == '__main__':
                smoothing_update_list = []
                for shell_no in shells_to_update:
                    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                        fout.write(f'  - Starting process for shell {shell_no}; total: {len(shells_to_update)} shells to update\n')
                    p = mp.Process(target = update_smoothing_block, args = (shell_no,))
                    p.start()
                    smoothing_update_list.append(p)
                    
                for smoothing_update_process in smoothing_update_list:
                    smoothing_update_process.join()

            smoothing_radii_time = time.time() - smoothing_radii_start
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Finished updating smoothing blocks; runtime: {smoothing_radii_time / 60} minutes / {(smoothing_radii_time / 60) / 60} hours\n')
        
        
        # now, go directly to smoothing. loop through each block in each shell and grab the information in the shell_{shell}/block_{block} file.
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start model smoothing\n')
        smoothing_start = time.time()
        if __name__ == '__main__':
            smoothing_process_list = []
            for shell_no in all_shells:
                if shell_no in shells_to_update:
                    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                        fout.write(f'  - Starting process for shell {shell_no}; total of {len(shells_to_update)} to smooth\n')
                    p = mp.Process(target = model_smoothing, args = (shell_no,))
                    p.start()
                    smoothing_process_list.append(p)
                
                else:
                    df_updated_model_shell = df_model.loc[df_model['SHELL#'] == shell_no].copy().reset_index(drop = True)
                    df_updated_model_shell.to_csv(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv', index = False)

            for smoothing_process in smoothing_process_list:
                smoothing_process.join()

        smoothing_time = time.time() - smoothing_start
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished model smoothing; runtime: {smoothing_time} seconds / {smoothing_time / 60} minutes\n')

        
        df_model = df_model.drop(df_model.index)
        for shell_no in all_shells:
            df_updated_shell = pd.read_csv(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv')
            df_model = pd.concat([df_model, df_updated_shell])
            # os.remove(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv')
        df_model = df_model.reset_index(drop = True)
        df_model.to_csv(f'{iteration_directory_path}/{model}_updated_model_itr_{layer_iteration}.csv', index = False)

        ###############################################################################################
        ### FINAL STEP: RESET THE BLOCK INFO FILES, DF_RMS AND DF_MODEL TO START THE NEXT ITERATION ###
        ###############################################################################################
        # Now, with the updated model (df_model), compute a new RMS profile
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Resetting block info files, starting rms, and starting model for the next iteration\n')
        rms_values = []
        normalized_rms_values = []

        for shell_no in all_shells:
            df_shell_rms = df_model.loc[df_model['SHELL#'] == shell_no]
            rms_value = rms(list(df_shell_rms[out_property_header]))
            df_rms_idx = df_rms.loc[df_rms['SHELL#'] == shell_no].index
            df_rms.loc[df_rms_idx, RMS_header] = rms_value

        max_rms_value = df_rms[RMS_header].max()
        df_rms[f'NORM_{RMS_header}'] = df_rms[RMS_header] / max_rms_value
        df_rms.to_csv(f'{iteration_directory_path}/{model}_updated_rms_itr_{layer_iteration}.csv', index = False)

        # make plottable files for the current iteration of the solution model and the difference between this model and the starting model
        try:
            os.mkdir(f'{iteration_directory_path}/plot_files')
        except:
            pass

        lats = np.arange(mod_input.start_lat, mod_input.final_lat, mod_input.reference_lat)
        lons = np.arange(mod_input.start_lon, mod_input.final_lon, mod_input.reference_lon)
        
        for shell_to_plot in all_shells:
            up_perturbs_file = f'{iteration_directory_path}/plot_files/{model}_shell_{shell_to_plot}_updated_perturbs_plot_ready_{job_id}.csv'
            perturb_diffs_file = f'{iteration_directory_path}/plot_files/{model}_shell_{shell_to_plot}_perturb_diffs_plot_ready_{job_id}.csv'

            df_og_shell_data = df_original_model.loc[df_original_model['SHELL#'] == shell_to_plot]
            df_model_shell_data = df_model.loc[df_model['SHELL#'] == shell_to_plot]
            if layer_iteration == 1:
                radii_file = f'{iteration_directory_path}/plot_files/{model}_shell_{shell_to_plot}_smoothing_radii_plot_ready_{job_id}.csv'
                if shell_to_plot in shells_to_update:
                    df_model_radii = pd.read_csv(f'{layer_directory}/shell_{shell_to_plot}_smoothing_radii.csv')
                else:
                    df_model_radii = pd.DataFrame(columns = ['BLOCK#', 'RADIUS'])
                    df_model_radii['BLOCK#'] = all_blocks
                    df_model_radii['RADIUS'] = 0.
                df_radii = pd.DataFrame(data = {'LON': lons})
            
            df_updated = pd.DataFrame(data = {'LON': lons})
            df_differences = pd.DataFrame(data = {'LON': lons})
            
            # make empty lists that will be filled with numpy arrays for each latitude band, for each model.
            up_dvs = []
            diff_dvs = []
            radii = []

            # loop through all latitudes in the model space
            for lat in lats:
                # make empty lists to fill in the perturbations for the current latitude band
                lat_up_dvs = []
                lat_diff_dvs = []
                lat_radii = []
            
                # loop through all longitudes in the model space
                for lon in lons:
                    # find the block that the current lat/lon pair falls in
                    block = mod_pandas.find_block_id(lat, lon)
                    original_perturb = float(df_og_shell_data.loc[df_og_shell_data['BLOCK#'] == block][out_property_header])
                    updated_perturb = float(df_model_shell_data.loc[df_model_shell_data['BLOCK#'] == block][out_property_header])
                    diff_perturb = updated_perturb - original_perturb
                    if layer_iteration == 1:
                        radius = int(df_model_radii.loc[df_model_radii['BLOCK#'] == block]['RADIUS']) - 0.1
                        lat_radii.append(radius)
                    
                    lat_up_dvs.append(updated_perturb)
                    lat_diff_dvs.append(diff_perturb)
                
                up_dvs.append(np.array(lat_up_dvs))
                diff_dvs.append(np.array(lat_diff_dvs))
                
                df_updated[f'{lat}'] = lat_up_dvs
                df_differences[f'{lat}'] = lat_diff_dvs

                if layer_iteration == 1:
                    radii.append(np.array(lat_radii))
                    df_radii[f'{lat}'] = lat_radii
    
            df_updated.to_csv(up_perturbs_file, index = False)
            df_differences.to_csv(perturb_diffs_file, index = False)

            if layer_iteration == 1:
                df_radii.to_csv(radii_file, index = False)

        iteration_time = time.time() - iteration_start

        # with open(f'{iteration_directory_path}/{model}_update_info.txt', 'w') as fout:
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'total time for iteration: {iteration_time / 60} minutes / {iteration_time / 60 / 60} hours\n')
            fout.write(f'total time for backmapping & converting from path to block format: {backmapping_time / 60} minutes / {(backmapping_time / 60) / 60} hours\n')
            fout.write(f'total time for merging all temporary block information into main block files: {merge_time / 60} minutes\n')
            fout.write(f'total time for compiling all unexplained differences and computing variance: {residual_time / 60} minutes\n')
            fout.write(f'total time for finding smoothing radii, Gaussian weights, and saving smoothing files OR updating smoothing files: {smoothing_radii_time / 60} minutes / {(smoothing_radii_time / 60) / 60} hours\n')
            fout.write(f'total time for smoothing: {smoothing_time / 60} minutes\n')
            fout.write(f'mean time for backmapping per path (total paths: {tot_paths}): {backmapping_time / tot_paths} seconds\n')
            fout.write(f'mean time for smoothing per shell: {smoothing_time / len(shells_to_update)} seconds\n')
            fout.write(f'mean time for smoothing per block: {(smoothing_time / len(shells_to_update)) / len(all_blocks)} seconds\n')
            fout.write(f'mean misfit (dataset residual - predicted residual): {unexplained_means[-1]} seconds \n')
            fout.write(f'misfit standard deviation: {unexplained_stds[-1]}\n')
            fout.write(f'misfit variance: {unexplained_vars[-1]}\n\n\n\n')

        
        for phase in phases_to_update:
            shutil.rmtree(f'./{mod_input.phases_directory}/{phase}/{mod_input.backmapped_paths_directory}_{job_id}')

    # this level is where a model will just have completed a full set of iterations and crossed the variance threshold. save here (df_model and df_rms) for individual layers.
    for list_index in range(len(paths_lists)):
        shutil.rmtree(f'{layer_directory}/tmp_block_info_idx_{list_index}') ## LAYER STRIPPING DEBUG
    
    df_layer_variance = pd.DataFrame(data = {'layer': layers_complete, 'iteration': list(range(1, len(unexplained_means) + 1, 1)), 'total_paths': total_paths_itr, 'misfit_mean': unexplained_means, 'misfit_std': unexplained_stds, 'misfit_var': unexplained_vars})
    df_variance = pd.concat([df_variance, df_layer_variance])
    df_variance = df_variance.reset_index(drop = True)
    
    df_model.to_csv(f'{layer_directory}/{model}_updated_model.csv', index = False)
    df_rms.to_csv(f'{layer_directory}/{model}_updated_rms.csv', index = False)

    layer_end = time.time()
    layer_time = layer_end - layer_start

df_model.to_csv(f'{update_path}/{model}_final_updated_model_{job_id}.csv', index = False)
df_rms.to_csv(f'{update_path}/{model}_final_updated_rms_{job_id}.csv', index = False)
df_variance.to_csv(f'{update_path}/{model}_variance_reduction_{job_id}.csv', index = False)
df_residuals.to_csv(f'{update_path}/{model}_final_residuals_{job_id}.csv', index = False)
shutil.rmtree(f'{block_centric_path}')
    
## this level is where a model will be fully complete. save it here. if the last layer in the layer stripping variable is set to the whole mantle,
#### this saved file will be the same as the last one from the previous indent.
end = time.time()
update_time = end - start
