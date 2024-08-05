import multiprocessing as mp
import mod_refmodels
import pandas as pd
import numpy as np
import mod_database
import mod_input
import mod_track
import mod_geo
import random
import shutil
import math
import glob
import time
import sys
import csv
import os

job_id = os.getpid()
start_time = time.time()
start_subj = f'Process began (PID: {job_id}); part4.model_update.py'
start_text = f'Process {job_id} began\nModel: {mod_input.input_model}\n'
try:
    mod_track.SendMsg(start_subj, start_text)
except:
    pass

cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places
sectors = mod_input.azimuthal_sectors
sector_extent = 180. / mod_input.azimuthal_sectors
sector_ids = list(range(1, sectors + 1))
model = mod_input.input_model
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
df_rms = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.starting_RMS_model_to_use}_update/{mod_input.starting_RMS_model_to_use}_RMS.csv')

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
    return int(df_az_slice['ID'].iloc[0])

# Find the standard deviation of the Gaussian function that will be used for weighting, and then define the Gaussian as a function of radius:
def gaussian(total_R, x, gaussian_cutoff_weight):
    gaussian_numerator = ((total_R) ** 2.)
    gaussian_denominator = np.log(gaussian_cutoff_weight)
    gaussian_standard_deviation = np.sqrt((gaussian_numerator / gaussian_denominator) * -0.5)
    exp_n = ((x) ** 2.)
    exp_d = (2 * gaussian_standard_deviation ** 2.)
    g = np.exp(-(exp_n / exp_d))
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

# initialize an empty dataframe for variance reduction
df_variance = pd.DataFrame(columns = ['layer', 'iteration', 'total_paths', 'misfit_mean', 'misfit_std', 'misfit_var'])


# define functions for multiprocessing/parallelization
def remove_mean(lists_of_paths, list_idx):
    df_main_data = pd.DataFrame(columns = ['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'PATH_ID', 'PHASE', layer_residual_header, 'DATASET', 'DT_PRED'])
    for list_of_paths in lists_of_paths:
        phase = str(list_of_paths[0].split('/')[2])
        phase_name = str(phase.split('_')[0])
        dataset = '_'.join(phase.split('_')[1:])
        df_data = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_master_data.csv')
        df_data = df_data[['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'PATH_ID', 'PHASE', layer_residual_header]].copy()
        df_data['PHASE'] = phase_name
        df_data['DATASET'] = dataset
        df_data['DT_PRED'] = np.NaN
        df_data = np.array(df_data).T
        
        for path in list_of_paths:
            remove_mean = True
            tot_time = 0.
            df_p = pd.read_csv(path)
            path_id = int(path.split('_')[-3])
            dataset_idx = np.where(df_data[3] == path_id)[0][0]
            # check if there are any segments in df_p that have a length of zero and drop these from the dataframe
            df_p = df_p.loc[df_p['SEG_DIST_KM'] != 0.].reset_index(drop = True)
            df_p = np.array(df_p)
            path_max_depth = df_p.T[1].max()

            # check if the current path is above the layer stripping depth. continue only if it is.
            if path_max_depth <= layer_stripping_depth and path_max_depth > layer_depth_min:

                # find the reference model (e.g., prem) total travel time for the path. this will just be the total time in df_p.
                ref_model_time = df_p.T[25, -1]

                # grab the dataset residual for the current path
                dataset_residual = df_data[5, dataset_idx]
                
                if layer_residual_limits[0] == True:
                    if layer_residual_limits[1] <= dataset_residual <= layer_residual_limits[2]:
                        remove_mean = True
                    else:
                        remove_mean = False
                
                elif layer_residual_limits[0] == False:
                    remove_mean = True
                    
                if remove_mean == True:
                    ## LOOP 1: compute the travel time of the current path through the starting model (df_model).
                    for segment in range(len(df_p)):
                        # for each path segment in the current working path, find its attributes
                        dist_km = df_p[segment, 17]
                        shell_no = df_p[segment, 19]
                        block_no = df_p[segment, 20]
                        mid_depth = df_p[segment, 21]
                        prem_val = df_p[segment, 22]

                        # df_model_slice is the single line in the model dataframe for the element that the current path segment is in.
                        # from df_model_slice, find the velocity perturbation of that element (model_element_perturb)
                        model_element_perturb = float(df_model.loc[(df_model['SHELL#'] == shell_no) & (df_model['BLOCK#'] == block_no)][out_property_header].iloc[0])

                        # find the velocity value for the path segment through the starting model.
                        seg_v = prem_val + ((model_element_perturb / 100.) * prem_val)
                        df_p[segment, 22] = seg_v

                        # use that velocity value and the length of the path (in km) to find the travel time for just this path segment
                        seg_time = dist_km / seg_v
                        df_p[segment, 18] = seg_time

                        # add that path segment travel time to the cumulative `tot_time` variable
                        tot_time += seg_time
                        df_p[segment, 25]

                    # update df_p, which is now the path through the model
                    df_p[:, 24] = (df_p[:, 25] / df_p[-1, 25]) * 100

                    # total starting model time - total reference model time (residual between starting model and reference model)
                    predicted_residual = tot_time - ref_model_time
                    df_data[-1, dataset_idx] = predicted_residual
        df_data = df_data.T
        df_data = pd.DataFrame(df_data, columns = ['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'PATH_ID', 'PHASE', layer_residual_header, 'DATASET', 'DT_PRED'])
        df_main_data = pd.concat([df_main_data, df_data]).reset_index(drop = True)
    df_main_data['EQ_DATE'] = df_main_data['EQ_DATE'].astype(int)
    df_main_data['PATH_ID'] = df_main_data['PATH_ID'].astype(int)
    df_main_data = df_main_data.dropna(subset = ['DT_PRED']).reset_index(drop = True)
    df_main_data.to_csv(f'{layer_directory}/tmp_pred_residuals_{list_idx}.csv', index = False)


def backmap_residual(lists_of_paths, list_idx):
    df_idx = 0
    block_info_cols = ['SHELL#', 'BLOCK#', 'PHASE', 'PATH_ID', 'PATH_LENGTH_KM', 'CENTER_LAT', 'CENTER_LON', out_property_header, 'AZIMUTH', 'SECTOR']
    data_cols = ['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'PATH_ID', layer_residual_header]
    special_weight_ids = []
    special_weight_idx = 5
    for special_weight in layer_special_weights:
        block_info_cols.append(special_weight)
        data_cols.append(special_weight)
        special_weight_ids.append(special_weight_idx)
        special_weight_idx += 1
    df_block_info = pd.DataFrame(columns = block_info_cols)
    df_residuals = pd.DataFrame(columns = ['PHASE', 'PATH_ID', 'REF_pred_time', 'ORIG_DT', 'STARTING_pred_time', 'STARTING_DT', 'MISFIT'])
    
    for list_of_paths in lists_of_paths:
        phase = str(list_of_paths[0].split('/')[2])
        phase_name = str(phase.split('_')[0])
        dataset = '_'.join(phase.split('_')[1:])
        df_data = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_master_data.csv').copy()
        df_data = np.array(df_data[data_cols])

        # apply the mean removal
        if layer_remove_mean == True:
            for line in range(len(df_data)):
                eq_date = df_data[line, 0]
                eq_lat = df_data[line, 1]
                eq_lon = df_data[line, 2]
                # find the corresponding mean correction in df_events
                try:
                    df_events_idx = np.where((df_events[0] == eq_date) & (df_events[1] == eq_lat) & (df_events[2] == eq_lon) & (df_events[3] == dataset))[0][0]
                    mean_corr = df_events[-1, df_events_idx]
                    df_data[line, 4] = df_data[line, 4] - mean_corr
                except:
                    pass
                
        df_data = df_data.T       
        for path in list_of_paths:
            backmap = True
            tot_time = 0.
            df_p = pd.read_csv(path)
            path_cols = list(df_p.columns)
            path_id = int(path.split('_')[-3])
            dataset_idx = np.where(df_data[3] == path_id)[0][0]
            # check if there are any segments in df_p that have a length of zero and drop these from the dataframe
            df_p = df_p.loc[df_p['SEG_DIST_KM'] != 0.].reset_index(drop = True)
            df_p = np.array(df_p)
            path_max_depth = df_p.T[1].max()

            # check if the current path is above the layer stripping depth. continue only if it is.
            if path_max_depth <= layer_stripping_depth and path_max_depth > layer_depth_min:

                # find the reference model (e.g., prem) total travel time for the path. this will just be the total time in df_p.
                ref_model_time = df_p.T[25, -1]

                # grab the dataset residual for the current path
                dataset_residual = df_data[4, dataset_idx]
                
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
                        dist_km = df_p[segment, 17]
                        shell_no = df_p[segment, 19]
                        block_no = df_p[segment, 20]
                        mid_depth = df_p[segment, 21]
                        prem_val = df_p[segment, 22]

                        # df_model_slice is the single line in the model dataframe for the element that the current path segment is in.
                        # from df_model_slice, find the velocity perturbation of that element (model_element_perturb)
                        model_element_perturb = float(df_model.loc[(df_model['SHELL#'] == shell_no) & (df_model['BLOCK#'] == block_no)][out_property_header].iloc[0])

                        # find the velocity value for the path segment through the starting model.
                        seg_v = prem_val + ((model_element_perturb / 100.) * prem_val)
                        df_p[segment, 22] = seg_v

                        # use that velocity value and the length of the path (in km) to find the travel time for just this path segment
                        seg_time = dist_km / seg_v
                        df_p[segment, 18] = seg_time

                        # add that path segment travel time to the cumulative `tot_time` variable
                        tot_time += seg_time
                        df_p[segment, 25]

                    # update df_p, which is now the path through the model
                    df_p[:, 24] = (df_p[:, 25] / df_p[-1, 25]) * 100

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
                    df_p = pd.DataFrame(df_p, columns = path_cols)
                    df_p['START_SHELL#'] = df_p['START_SHELL#'].astype(int)
                    df_p['START_BLOCK#'] = df_p['START_BLOCK#'].astype(int)
                    df_p['END_SHELL#'] = df_p['END_SHELL#'].astype(int)
                    df_p['END_BLOCK#'] = df_p['END_BLOCK#'].astype(int)
                    df_p['SEG_SHELL#'] = df_p['SEG_SHELL#'].astype(int)
                    df_p['SEG_BLOCK#'] = df_p['SEG_BLOCK#'].astype(int)

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
                        norm_rms_update = float(df_layer_rms.loc[df_layer_rms['SHELL#'] == segment_shell_no][f'NORM_{RMS_header}'].iloc[0])
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
                    df_p_update_np = np.array(df_p_update)
                    
                    # Isolate the unique block and shell IDS in each path and condense to representative block values
                    df_path_elements = np.array(df_p_update.groupby(['SEG_SHELL#','SEG_BLOCK#']).size().reset_index())
                    # loop through each unique combination
                    for pair in range(len(df_path_elements)):
                        shell_no = df_path_elements[pair, 0]
                        block_no = df_path_elements[pair, 1]

                        # find all of the segments in the current block
                        df_update_indices = np.where((df_p_update_np.T[19] == shell_no) & (df_p_update_np.T[20] == block_no))[0]
                        df_slice = df_p_update_np[df_update_indices]

                        updated_velocity = df_slice[0, 22]
                        mid_depth = df_slice[0, 21]
                        prem_velocity = mod_refmodels.prem_vel(mod_input.data_wave_type, mid_depth)
                        percent_perturbation = (100. * (updated_velocity / prem_velocity)) - 100.

                        # if we're in the first iteration, then add all new data. if not, just update the velocity perturbation column.
                        total_length_km = np.sum(df_slice.T[17])
                        total_length_deg = np.sum(df_slice.T[16])
                        start_lat = df_slice[0, 3]
                        start_lon = df_slice[0, 4]
                        end_lat = df_slice[-1, 10]
                        end_lon = df_slice[-1, 11]
                        azimuth = np.mean(df_slice.T[15])
                        az_sector = azimuthal_sector(azimuth)
                        segment_center = total_length_deg / 2.
                        segment_center_coords = mod_geo.GCP_point(start_lat, start_lon, end_lat, end_lon, total_length_deg, segment_center)

                        # add all of that information to the last row of the block file
                        new_block_info = [shell_no, block_no, phase, path_id, total_length_km, segment_center_coords[0], segment_center_coords[1], percent_perturbation, azimuth, az_sector]

                        for special_weight_id in special_weight_ids:
                            new_block_info.append(df_data[special_weight_id, dataset_idx])
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
    df_model_arr = np.array(df_model)
    df_updated_model_shell = np.array(pd.DataFrame(data = {'SHELL#': shell_no, 'BLOCK#': all_blocks, out_property_header: 0.}))
    cols = ['SECTOR', out_property_header, 'GAUSS_WEIGHT']
    special_weight_id = 3
    special_weight_ids = []
    for special_weight in layer_special_weights:
        cols.append(special_weight)
        special_weight_ids.append(special_weight_id)
        special_weight_id += 1
        
    for block_no in all_blocks:
        df_smoothing = pd.read_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv')
        df_smoothing = df_smoothing[cols]
        df_smoothing = np.array(df_smoothing)

        # ONLY ENTER THE SMOOTHING LOOP if there is at least one segment in df_smoothing.
        # otherwise, if there are 0 paths, don't bother smoothing and just use the value of the current block.
        if len(df_smoothing) > 0:
            # make a histogram of all of the sectors in df_smoothing
            unique_path_segment_sectors_hist = np.histogram(df_smoothing.T[0], bins = np.arange(1, mod_input.azimuthal_sectors + 2, 1))[0]
        
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
                        df_final_smoothing = df_smoothing[np.isin(df_smoothing.T[0], (path_sectors + 1))]
                        pre_az_numerator = 0.
                        pre_az_denominator = 0.

                        for segment_final in range(len(df_final_smoothing)):
                            weights = []
                            segment_final_pert = df_final_smoothing[segment_final, 1]
                            weights.append(df_final_smoothing[segment_final, 2])

                            for special_weight_idx in special_weight_ids:
                                weight = df_final_smoothing[segment_final, special_weight_idx]
                                
                                if math.isnan(weight) == False and weight != 0.:
                                    weights.append(weight)

                            if math.isnan(segment_final_pert) == False:
                                pre_az_numerator += (segment_final_pert * np.prod(weights))
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
                    segment_final_pert = df_smoothing[segment_final, 1]
                    weights.append(df_smoothing[segment_final, 2])
                    
                    for special_weight_idx in special_weight_ids:
                        weight = df_smoothing[segment_final, special_weight_idx]
                        
                        if math.isnan(weight) == False and weight != 0.:
                            weights.append(weight)

                    if math.isnan(segment_final_pert) == False:
                        smoothing_numerator += (segment_final_pert * np.prod(weights))
                        smoothing_denominator += (np.prod(weights))
        
                if smoothing_denominator != 0.:
                    mean_element_perturbation = smoothing_numerator / smoothing_denominator
        else:
            mean_element_perturbation_idx = np.where((df_model_arr.T[0] == shell_no) & (df_model_arr.T[1] == block_no))[0][0]
            mean_element_perturbation = df_model_arr[mean_element_perturbation_idx, -1]
    
        element_idx = np.where((df_updated_model_shell.T[0] == shell_no) & (df_updated_model_shell.T[1] == block_no))[0][0]
        df_updated_model_shell[element_idx, -1] = mean_element_perturbation

    df_updated_model_shell = pd.DataFrame(df_updated_model_shell, columns = ['SHELL#', 'BLOCK#', out_property_header])
    df_updated_model_shell['SHELL#'] = df_updated_model_shell['SHELL#'].astype(int)
    df_updated_model_shell['BLOCK#'] = df_updated_model_shell['BLOCK#'].astype(int)
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
    df_smoothing_radii = np.array(pd.read_csv(f'{layer_directory}/shell_{shell_no}_smoothing_radii.csv'))
    cols = list(df_updated_shell.columns)
    cols.append('GAUSS_WEIGHT')
    
    for block_no in all_blocks:
        # find the center lat/lon of the current block
        block_center_lat = mod_database.get_block_info(block_no)[6]
        block_center_lon = mod_database.get_block_info(block_no)[7]
    
        # find the block's smoothing radius
        block_smoothing_radius_idx = np.where(df_smoothing_radii.T[0] == block_no)[0]
        block_smoothing_radius = df_smoothing_radii[block_smoothing_radius_idx, 1][0]
    
        # open the file with all of the raypath segments to be smoothed for this block
        df_smoothing_block = pd.read_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv')
        neighbors = df_smoothing_block['NEIGHBOR'].unique()
        df_smoothing_block = df_smoothing_block.drop(df_smoothing_block.index)
    
        for neighbor in neighbors:
            df_updated_block = df_updated_shell.loc[df_updated_shell['BLOCK#'] == neighbor].rename(columns = {'BLOCK#': 'NEIGHBOR'})
            df_smoothing_block = pd.concat([df_smoothing_block, df_updated_block])
        df_smoothing_block = df_smoothing_block.reset_index(drop = True)
        df_smoothing_block['GAUSS_WEIGHT'] = 0.
        df_smoothing_block = np.array(df_smoothing_block)
            
        for segment_id in range(len(df_smoothing_block)):
            segment_center_lat = df_smoothing_block[segment_id, 4]
            segment_center_lon = df_smoothing_block[segment_id, 5]
            
            dist_to_center = mod_geo.GCP_length(block_center_lat, block_center_lon, segment_center_lat, segment_center_lon)
            gaus_weight = gaussian(block_smoothing_radius, dist_to_center, layer_gaussian_cutoff_weight)
            if gaus_weight == 0.:
                print(gaus_weight)
            df_smoothing_block[segment_id, -1] = gaus_weight
    
        df_smoothing_block = pd.DataFrame(df_smoothing_block, columns = cols).rename(columns = {'BLOCK#': 'NEIGHBOR'})
        df_smoothing_block['NEIGHBOR'] = df_smoothing_block['NEIGHBOR'].astype(int)
        df_smoothing_block['PATH_ID'] = df_smoothing_block['PATH_ID'].astype(int)
        df_smoothing_block['SECTOR'] = df_smoothing_block['SECTOR'].astype(int)
        
        df_smoothing_block.to_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv', index = False)


def find_smoothing_radii(shell_no):
    df_shell_elements = pd.read_csv(f'{block_centric_path}/shell_{shell_no}.csv')
    df_shell_smoothing_radii = np.array(pd.DataFrame(data = {'BLOCK#': list(all_blocks), 'RADIUS': 0.}))
    df_shell_elements['GAUSS_WEIGHT'] = 0.
    cols = list(df_shell_elements.columns)
    
    for block_no in all_blocks:
        # access the main neighbors file for this block.
        element_center_lat = mod_database.get_block_info(block_no)[6]
        element_center_lon = mod_database.get_block_info(block_no)[7]
        df_neighbors = np.array(pd.read_csv(f'./{mod_input.near_neighbors_directory}/block_{block_no}_neighbors.csv'))

        radius_found = False
        for radius in layer_smoothing_radii:
            # find near neighbors that fall within the current working radius and make a dataframe with all of the paths that sample those neighbors (df_smoothing)
            near_neighbor_indices = np.where(df_neighbors.T[1] <= radius)[0]
            eligible_blocks = df_neighbors[near_neighbor_indices, 0]
            df_smoothing = df_shell_elements.loc[df_shell_elements['BLOCK#'].isin(eligible_blocks)].reset_index(drop = True)

            # if the current radius is equal to the maximum radius, just define df_smoothing based on that radius and don't faff with the other criteria
            if radius == layer_max_radius:
                radius_found = True

            # if not, check coverage criteria
            else:
                # df_unique_paths is a dataframe that has the total number of segments for each full path within the current radius.
                df_unique_paths = df_smoothing[['PHASE', 'PATH_ID', 'AZIMUTH']].groupby(['PHASE','PATH_ID']).mean()['AZIMUTH']
                
                # first, check if there are enough total paths in the radius:
                # if there aren't, move on to the next radius
                if len(df_unique_paths) < layer_total_required_paths:
                    pass
                
                else:
                    sectors_covered = []

                    # find the azimuthal sector of each unique path azimuth
                    for unique_az in df_unique_paths:
                        unique_path_sector = azimuthal_sector(unique_az)

                        # if that sector isn't yet represented, append it to sectors_covered
                        if unique_path_sector not in sectors_covered:
                            sectors_covered.append(unique_path_sector)

                            # if the number of sectors covered meets the criteria, break
                            if len(sectors_covered) >= layer_total_required_azimuths:
                                radius_found = True
                                break
                            
            if radius_found == True:
                break

        # now that the smoothing radius and the paths/path segments within that smoothing radius for the current block has been determined, find the Gaussian weights for each of those segments
        if len(df_smoothing) > 0:
            # find the Gaussian weight for all segments in the radius to be smoothed (i.e., all segments in df_smoothing)
            df_smoothing = np.array(df_smoothing)
            
            for unique_path_segment in range(len(df_smoothing)):
                unique_path_segment_center_lat = df_smoothing[unique_path_segment, 4]
                unique_path_segment_center_lon = df_smoothing[unique_path_segment, 5]
    
                # find the distance from the center of the current element to the center of the path segment
                dist_to_center = mod_geo.GCP_length(element_center_lat, element_center_lon, unique_path_segment_center_lat, unique_path_segment_center_lon)
                
                # compute the gaussian weight for that path segment and add it to df_smoothing
                gaus_weight = gaussian(radius, dist_to_center, layer_gaussian_cutoff_weight)
                df_smoothing[unique_path_segment, -1] = gaus_weight
    
        df_smoothing = pd.DataFrame(df_smoothing, columns = cols).rename(columns = {'BLOCK#': 'NEIGHBOR'})
        df_smoothing['NEIGHBOR'] = df_smoothing['NEIGHBOR'].astype(int)
        df_smoothing['PATH_ID'] = df_smoothing['PATH_ID'].astype(int)
        df_smoothing['SECTOR'] = df_smoothing['SECTOR'].astype(int)
        df_smoothing.to_csv(f'{layer_directory}/smoothing_block_information/shell_{shell_no}/block_{block_no}.csv', index = False)
        
        smoothing_idx = np.where(df_shell_smoothing_radii.T[0] == block_no)
        df_shell_smoothing_radii[smoothing_idx, 1] = radius
        
    df_shell_smoothing_radii = pd.DataFrame(df_shell_smoothing_radii, columns = ['BLOCK#', 'RADIUS'])
    df_shell_smoothing_radii['BLOCK#'] = df_shell_smoothing_radii['BLOCK#'].astype(int)
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
    layer_stripping_depth = mod_database.get_shell_info(layer_bottom_shell)[3]
    layer_depth_min = mod_database.get_shell_info(layer_top_shell)[1]
    layer_remove_mean = mod_input.remove_residual_mean[layers_complete]
    freeze_previous_layer = mod_input.freeze_previous_layers[layers_complete]
    phases_to_update = mod_input.update_phases[layers_complete]
    type_of_data_subselection = mod_input.type_of_phase_subselection[layers_complete]
    subselection_of_phases_to_update = mod_input.subselection_of_phase_data_to_use[layers_complete]
    iteration_to_stop_RMS = mod_input.iteration_to_stop_RMS_weighting[layers_complete]
    shells_to_update = list(df_shells.loc[(df_shells['SHELL#'] >= layer_top_shell) & (df_shells['SHELL#'] <= layer_bottom_shell)]['SHELL#'])
    layer_residual_limits = mod_input.residual_limits[layers_complete]
    layer_cutoff_type = mod_input.cutoff_type[layers_complete]
    layer_cutoff = mod_input.cutoff[layers_complete]
    layer_smoothing_radii = mod_input.smoothing_radii[layers_complete]
    layer_max_radius = layer_smoothing_radii[-1]
    layer_total_required_paths = mod_input.total_required_paths[layers_complete]
    layer_total_required_azimuths = mod_input.total_required_azimuths[layers_complete]
    layer_gaussian_cutoff_weight = mod_input.gaussian_cutoff_weight[layers_complete]
    layer_azimuthal_weighting = mod_input.azimuthal_weighting[layers_complete]
    layer_special_weights = mod_input.special_weights[layers_complete]
    layer_residual_header = mod_input.residual_header[layers_complete]
    layer_dataset_description = mod_input.dataset_description[layers_complete]
    
    layers_complete += 1
    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
        fout.write(f'\n-----\n')
        fout.write(f'-----   STARTING UPDATE FOR LAYER {layers_complete} of {len(mod_input.layer_base_shells)}, DEPTH RANGE: {layer_depth_min} km - {layer_stripping_depth} km   -----\n')
        fout.write(f'-----\n')
    
    # write to model update log
    if os.path.exists(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv'):
        with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv', 'a') as fout:
            fout.write(f'{job_id},')
            fout.write(f'{layers_complete} of {len(mod_input.layer_base_shells)},')
            fout.write(f'{model},')
            fout.write(f'{mod_input.starting_RMS_model_to_use},')
            fout.write(f'{mod_input.reference_model},')
            fout.write(f'{layer_dataset_description},')
            fout.write(f'{layer_residual_header},')
            fout.write(f'{layer_remove_mean},')
            fout.write(f"[{';'.join(str(i) for i in phases_to_update)}],")
            fout.write(f'{type_of_data_subselection},')
            fout.write(f"[{';'.join(str(i) for i in subselection_of_phases_to_update)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_residual_limits)}],")
            fout.write(f'{layer_depth_min} km - {layer_stripping_depth} km,')
            fout.write(f'{iteration_to_stop_RMS},')
            fout.write(f'{layer_cutoff_type},')
            fout.write(f'{layer_cutoff},')
            fout.write(f'{mod_input.azimuthal_sectors},')
            fout.write(f'{layer_azimuthal_weighting},')
            fout.write(f"[{';'.join(str(i) for i in layer_special_weights)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_smoothing_radii)}],")
            fout.write(f'{layer_total_required_paths},')
            fout.write(f'{layer_total_required_azimuths},')
            fout.write(f'{layer_gaussian_cutoff_weight}\n')

    else:
        with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv', 'w') as fout:
            fout.write('PID,Layer,Input_model,Input_RMS,Reference_model,Dataset_description,Residual_header,Remove_mean,Phases_used,Type_of_phase_subselection,Phase_subselection,Residual_limits,Layer_dimensions,Stopped_RMS_weighting_iteration,Cutoff_type,Cutoff,Azimuthal_sectors,Azimuthal_weighting,Special_weights,Smoothing_radii,Total_required_paths,Total_required_sectors,Gaussian_cutoff_weight\n')
        with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model}_update/model_update_log_{model}.csv', 'a') as fout:
            fout.write(f'{job_id},')
            fout.write(f'{layers_complete} of {len(mod_input.layer_base_shells)},')
            fout.write(f'{model},')
            fout.write(f'{mod_input.starting_RMS_model_to_use},')
            fout.write(f'{mod_input.reference_model},')
            fout.write(f'{layer_dataset_description},')
            fout.write(f'{layer_residual_header},')
            fout.write(f'{layer_remove_mean},')
            fout.write(f"[{';'.join(str(i) for i in phases_to_update)}],")
            fout.write(f'{type_of_data_subselection},')
            fout.write(f"[{';'.join(str(i) for i in subselection_of_phases_to_update)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_residual_limits)}],")
            fout.write(f'{layer_depth_min} km - {layer_stripping_depth} km,')
            fout.write(f'{iteration_to_stop_RMS},')
            fout.write(f'{layer_cutoff_type},')
            fout.write(f'{layer_cutoff},')
            fout.write(f'{mod_input.azimuthal_sectors},')
            fout.write(f'{layer_azimuthal_weighting},')
            fout.write(f"[{';'.join(str(i) for i in layer_special_weights)}],")
            fout.write(f"[{';'.join(str(i) for i in layer_smoothing_radii)}],")
            fout.write(f'{layer_total_required_paths},')
            fout.write(f'{layer_total_required_azimuths},')
            fout.write(f'{layer_gaussian_cutoff_weight}\n')

    layer_directory = f'{update_path}/layer_{layers_complete}_shell_{layer_top_shell}_to_{layer_bottom_shell}'
    try:
        os.mkdir(layer_directory)
    except:
        pass

    for shell_no in all_shells:
        with open(f'{block_centric_path}/shell_{shell_no}.csv', 'w') as file:
            file.write(f'BLOCK#,PHASE,PATH_ID,PATH_LENGTH_KM,CENTER_LAT,CENTER_LON,{out_property_header},AZIMUTH,SECTOR\n')
    
    # loop through all of the paths for all of the phases to be updated in this layer and split them up into smaller lists for multiprocessing.
    with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
        fout.write(f'- Partitioning paths into batches for backmapping\n')
    paths_lists = [] # a list of length `cores_requested` with phase-wise sublists to be processed on each CPU core
    paths_list = [] # a single list of ALL paths for ALL phases to be updated.

    for phase_idx in range(len(phases_to_update)):
        phase = phases_to_update[phase_idx]
        phase_subselection = subselection_of_phases_to_update[phase_idx]
        phase_name = phase.split('_')[0]
        all_paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.resampled_directory}/{phase_name}_*.csv')
        if type_of_data_subselection == 'proportion':
            subset_of_paths = int(len(all_paths) * phase_subselection)
            paths = random.sample(all_paths, subset_of_paths)
        elif type_of_data_subselection == 'number':
            if len(all_paths) < phase_subselection:
                paths = all_paths
            else:
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
            os.mkdir(f'{layer_directory}/tmp_block_info_idx_{list_index}')
        except:
            pass
    
    unexplained_means = []
    unexplained_stds = []
    unexplained_vars = []
    total_paths_itr = []

    ## PART 0: remove the mean from the residuals of each event based on the starting tomography model
    if layer_remove_mean == True:
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start removing mean from residuals for a proxy event relocation\n')
        rm_mean_start = time.time()
        
        if __name__ == '__main__':
            process_list = []
            process_idx = 0
            for list_index in range(len(paths_lists)):
                process_idx += 1
                p = mp.Process(target = remove_mean, args = (paths_lists[list_index], list_index,))
                p.start()
                process_list.append(p)
    
            for process in process_list:
                process.join()
        
        df_pred = pd.DataFrame(columns = ['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'PATH_ID', 'PHASE', layer_residual_header, 'DATASET', 'DT_PRED'])
        for list_index in range(len(paths_lists)):
            df_pred_idx = pd.read_csv(f'{layer_directory}/tmp_pred_residuals_{list_index}.csv')
            df_pred = pd.concat([df_pred, df_pred_idx])
            os.remove(f'{layer_directory}/tmp_pred_residuals_{list_index}.csv')
    
        df_events = df_pred.groupby(['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'DATASET']).size().reset_index().drop(columns = [0])
        df_events['MEAN_DT'] = 0.
        df_events = np.array(df_events)
        for line in range(len(df_events)):
            eq_date = df_events[line, 0]
            eq_lat = df_events[line, 1]
            eq_lon = df_events[line, 2]
            dataset = df_events[line, 3]
            df_pred_slice = df_pred.loc[(df_pred['EQ_DATE'] == eq_date) & (df_pred['EQ_LAT'] == eq_lat) & (df_pred['EQ_LON'] == eq_lon) & (df_pred['DATASET'] == dataset)]
            df_events[line, -1] = df_pred_slice['DT_PRED'].mean()
    
        df_events = pd.DataFrame(df_events, columns = ['EQ_DATE', 'EQ_LAT', 'EQ_LON', 'DATASET', 'MEAN_DT'])
        df_events.to_csv(f'{layer_directory}/mean_event_residuals.csv', index = False)
        df_events = np.array(df_events).T
        
        rm_mean_time = time.time() - rm_mean_start
    
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished removing mean; runtime: {mod_track.runtime(rm_mean_time)}\n')

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
        previous_iteration_directory_path = f'{layer_directory}/iteration_{layer_iteration - 1}'
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
                p = mp.Process(target = backmap_residual, args = (paths_lists[list_index], list_index,))
                p.start()
                process_list.append(p)

            for process in process_list:
                process.join()
        
        backmapping_time = time.time() - backmapping_start

        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished backmapping and converting from path-to-block format; runtime: {mod_track.runtime(backmapping_time)}\n')
        
        ## Merge all files for the different processes into the main block files:
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Start merging temporary block files into main block files\n')
        merge_start = time.time()
        
        if __name__ == '__main__':
            merge_process_list = []

            for shell_no in shells_to_update:
                p = mp.Process(target = merge_block_files, args = (shell_no,))
                p.start()
                merge_process_list.append(p)
                
            for merge_process in merge_process_list:
                merge_process.join()

        merge_time = time.time() - merge_start
        with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
            fout.write(f'- Finished merging temporary block files into main block files; runtime: {mod_track.runtime(merge_time)}\n')
        
        
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

        if layer_iteration == 1:
            df_residuals.to_csv(f'{layer_directory}/{model}_starting_residuals_{job_id}.csv', index = False)
        else:
            df_residuals.to_csv(f'{previous_iteration_directory_path}/residuals.csv', index = False)
        
        unexplained_diffs = df_residuals['MISFIT'].to_numpy()
        unexplained_means.append(np.mean(unexplained_diffs))
        unexplained_stds.append(np.std(unexplained_diffs))
        unexplained_vars.append(np.var(unexplained_diffs))
        total_paths_itr.append(len(df_residuals))
        
        if layer_iteration > 1:
            std1 = unexplained_stds[-1]
            std2 = unexplained_stds[-2]
            std_diff = std2 - std1

            var1 = unexplained_vars[-1]
            var2 = unexplained_vars[-2]
            var_diff = var2 - var1

            if layer_cutoff_type == 'reduction':
                if var_diff < layer_cutoff:
                    continue_layer_iterating = False

            elif layer_cutoff_type == 'total iterations':
                if layer_iteration == layer_cutoff + 1:
                    continue_layer_iterating = False

        if continue_layer_iterating == False:
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Finished computing variance; threshold reached and layer update ended\n')
            shutil.rmtree(iteration_directory_path)
            
        elif continue_layer_iterating == True:
            residual_time = time.time() - residual_start
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Finished computing variance; threshold not yet reached\n')
            
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
                    fout.write(f'- Finished finding smoothing radii; runtime: {mod_track.runtime(smoothing_radii_time)}\n')
    
            else:
                with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                    fout.write(f'- Start updating smoothing blocks\n')
                smoothing_radii_start = time.time()
                
                if __name__ == '__main__':
                    smoothing_update_list = []
                    for shell_no in shells_to_update:
                        p = mp.Process(target = update_smoothing_block, args = (shell_no,))
                        p.start()
                        smoothing_update_list.append(p)
                        
                    for smoothing_update_process in smoothing_update_list:
                        smoothing_update_process.join()
    
                smoothing_radii_time = time.time() - smoothing_radii_start
                with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                    fout.write(f'- Finished updating smoothing blocks; runtime: {mod_track.runtime(smoothing_radii_time)}\n')            
        
            # now, go directly to smoothing. loop through each block in each shell and grab the information in the shell_{shell}/block_{block} file.
            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'- Start model smoothing\n')
            smoothing_start = time.time()
            if __name__ == '__main__':
                smoothing_process_list = []
                for shell_no in all_shells:
                    if shell_no in shells_to_update:
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
                fout.write(f'- Finished model smoothing; runtime: {mod_track.runtime(smoothing_time)}\n')
    
            
            df_model = df_model.drop(df_model.index)
            for shell_no in all_shells:
                df_updated_shell = pd.read_csv(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv')
                df_model = pd.concat([df_model, df_updated_shell])
                os.remove(f'{iteration_directory_path}/{model}_updated_shell_{shell_no}.csv')
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
                        block = mod_database.find_block_id(lat, lon)
                        original_perturb = float(df_og_shell_data.loc[df_og_shell_data['BLOCK#'] == block][out_property_header].iloc[0])
                        updated_perturb = float(df_model_shell_data.loc[df_model_shell_data['BLOCK#'] == block][out_property_header].iloc[0])
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

            with open(f'{update_path}/pt7_{job_id}_update_log.txt', 'a') as fout:
                fout.write(f'total time for iteration: {mod_track.runtime(iteration_time)}\n')
                fout.write(f'total time for backmapping & converting from path to block format: {mod_track.runtime(backmapping_time)}\n')
                fout.write(f'total time for merging all temporary block information into main block files: {mod_track.runtime(merge_time)}\n')
                fout.write(f'total time for compiling all unexplained differences and computing variance: {mod_track.runtime(residual_time)}\n')
                fout.write(f'total time for finding smoothing radii, Gaussian weights, and saving smoothing files OR updating smoothing files: {mod_track.runtime(smoothing_radii_time)}\n')
                fout.write(f'total time for smoothing: {mod_track.runtime(smoothing_time)}\n')
                fout.write(f'mean time for backmapping per path (total paths: {tot_paths}): {mod_track.runtime(backmapping_time / tot_paths)}\n')
                fout.write(f'mean time for smoothing per shell: {mod_track.runtime(smoothing_time / len(shells_to_update))}\n')
                fout.write(f'mean time for smoothing per block: {mod_track.runtime((smoothing_time / len(shells_to_update)) / len(all_blocks))}\n')
            
        for phase in phases_to_update:
            shutil.rmtree(f'./{mod_input.phases_directory}/{phase}/{mod_input.backmapped_paths_directory}_{job_id}')

    # this level is where a model will just have completed a full set of iterations and crossed the variance threshold. save here (df_model and df_rms) for individual layers.
    for list_index in range(len(paths_lists)):
        shutil.rmtree(f'{layer_directory}/tmp_block_info_idx_{list_index}') ## LAYER STRIPPING DEBUG

    # merge all smoothing block files:
    df_smoothing_radii_master = pd.DataFrame(columns = ['SHELL#', 'BLOCK#', 'RADIUS'])
    for shell in shells_to_update:
        shell_list = []
        df_shell_smoothing = pd.read_csv(f'{layer_directory}/shell_{shell}_smoothing_radii.csv')
        for i in range(len(df_shell_smoothing['BLOCK#'])):
            shell_list.append(shell)
        df_shell_smoothing.insert(0, 'SHELL#', shell_list)
        df_smoothing_radii_master = pd.concat([df_smoothing_radii_master, df_shell_smoothing])
        os.remove(f'{layer_directory}/shell_{shell}_smoothing_radii.csv')
    df_smoothing_radii_master.to_csv(f'{layer_directory}/smoothing_radii.csv', index = False)

    # df_layer_variance = pd.DataFrame(data = {'layer': layers_complete, 'iteration': list(range(1, len(unexplained_means) + 1, 1)), 'total_paths': total_paths_itr, 'misfit_mean': unexplained_means, 'misfit_std': unexplained_stds, 'misfit_var': unexplained_vars})
    df_layer_variance = pd.DataFrame(data = {'layer': layers_complete, 'iteration': list(range(0, len(unexplained_means), 1)), 'total_paths': total_paths_itr, 'misfit_mean': unexplained_means, 'misfit_std': unexplained_stds, 'misfit_var': unexplained_vars})
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

end_subj = f'Process ended (PID: {job_id}); part4.model_update.py'
end_text = f'Process {job_id} complete;\nRuntime: {mod_track.runtime(time.time() - start_time)}'
try:
    mod_track.SendMsg(end_subj, end_text)
except:
    pass
