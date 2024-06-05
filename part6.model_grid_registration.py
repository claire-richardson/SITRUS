import numpy as np
import pandas as pd
import mod_input
import mod_geo
import mod_database
import mod_refmodels
import netCDF4
import glob
import time
import os
import multiprocessing as mp


# define relevant input variables from mod_input, including file names and paths
cdp = mod_input.computed_decimal_places
df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
all_shells = df_shells['SHELL#']
all_blocks = df_blocks['BLOCK#']
model_name = mod_input.input_model
lat_header = mod_input.lat_header
lon_header = mod_input.lon_header
depth_header = mod_input.depth_header
prop_header = mod_input.property_header
convert = mod_input.convert_vel_to_perturb


if mod_input.data_wave_type == 'S':
    out_property_header = 'dVs_%'
    RMS_header = 'RMS_dVs'
elif mod_input.data_wave_type == 'P':
    out_property_header = 'dVp_%'
    RMS_header = 'RMS_dVp'

def rms(property_list):
    n = len(property_list)
    x = 0
    for p in property_list:
        x += (p**2)
    return np.sqrt(x / n)

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'w') as fout:
    fout.write(f'Starting interpolation for model {model_name}\n')

# define the input model to interpolate
df_input_model = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/{model_name}.csv')
# compute the voigt average, if necessary
if mod_input.voigt[0] == True:
    df_input_model[prop_header] = np.sqrt(((2. * (df_input_model[mod_input.voigt[1]] ** 2)) + (df_input_model[mod_input.voigt[2]] ** 2)) / 3.)
df_model = df_input_model[[lat_header, lon_header, depth_header, prop_header]].copy()

if df_model[lon_header].max() >= 180.:
    for line in range(len(df_model)):
        model_long = df_model[lon_header].iloc[line]
        if model_long >= 180.:
            df_model.loc[line, lon_header] = model_long - 360.

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- finding all shell and block IDs for model nodes\n')

df_model['SHELL#'] = 0.
df_model['BLOCK#'] = 0.
df_model = np.array(df_model)
for line in range(len(df_model)):
    model_lat = df_model[line, 0]
    model_lon = df_model[line, 1]
    model_depth = df_model[line, 2]
    if model_depth < 0:
        pass
    else:
        model_shell = mod_database.find_shell_id(model_depth)
        try:
            model_block = mod_database.find_block_id(model_lat, model_lon)
        except:
            print(model_lat, model_lon, model_depth)
        df_model[line, 4] = model_shell
        df_model[line, 5] = model_block
    
with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- defining the new dataframe\n')

for shell in all_shells:
    model_length = np.zeros(len(all_blocks))
    idx = 0
    shell_info = mod_database.get_shell_info(shell)
    shell_mid = shell_info[2]
    df_registered_model = pd.DataFrame(data = {'SHELL#': shell, 'SHELL_MID': shell_mid, 'BLOCK#': model_length, 'BLOCK_LAT_MID': model_length, 'BLOCK_LON_MID': model_length, out_property_header: 0.})
    df_registered_model = np.array(df_registered_model)
    for block in all_blocks:
        block_info = mod_database.get_block_info(block)
        block_lat_mid = block_info[6]
        block_lon_mid = block_info[7]
        df_registered_model[idx, 2] = block
        df_registered_model[idx, 3] = block_lat_mid
        df_registered_model[idx, 4] = block_lon_mid
        idx += 1
    df_registered_model = pd.DataFrame(df_registered_model, columns = ['SHELL#', 'SHELL_MID', 'BLOCK#', 'BLOCK_LAT_MID', 'BLOCK_LON_MID', out_property_header])
    df_registered_model.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell}.csv', index = False)

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- interpolating the model\n')

def interpolate_model(shell_no):
    # reregister the model to our grid system by either finding the value at the center of our blocks or making an average of all model values that fall within each block
    df_registered_shell = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell_no}.csv')
    shell_mid = float(df_registered_shell['SHELL_MID'].iloc[0])
    for line in range(len(df_registered_shell)):
        block = df_registered_shell[line, 2]
        block_lat_mid = df_registered_shell[line, 3]
        block_lon_mid = df_registered_shell[line, 4]
        # first, check to see if there's a model value at this exact point.
        model_indices = np.where((df_model.T[0] == block_lat_mid) & (df_model.T[1] == block_lon_mid) & (df_model.T[2] == shell_mid))[0]
        if len(model_indices) == 1:
            df_registered_shell[line, 5] = df_model[model_indices][0][3]
        # if there isn't a model value at this exact point, calculate an average of all model values that fall into our block
        else:
            model_indices = np.where((df_model.T[4] == shell_no) & (df_model.T[5] == block))[0]
            if len(model_indices) == 0:
                df_registered_shell[line, 5] = 0.
            else:
                df_registered_shell[line, 5] = np.mean(df_model[model_indices, 3])

    if convert == True:
        update_ref_val = mod_refmodels.prem_vel(mod_input.data_wave_type, shell_mid)
        for idx in range(len(df_registered_shell)):
            block = df_registered_shell[idx, 2]
            vel = df_registered_shell[idx, 5]
            if vel == 0.:
                pass
            else:
                updated_perturbation = ((vel / update_ref_val) * 100.) - 100.
            df_registered_shell[idx, 5] = updated_perturbation

    df_registered_shell = pd.DataFrame(df_registered_shell, columns = ['SHELL#', 'SHELL_MID', 'BLOCK#', 'BLOCK_LAT_MID', 'BLOCK_LON_MID', out_property_header])
    df_registered_shell.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell}.csv', index = False)


model_registration_start = time.time()
if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for shell in all_shells:
        process_idx += 1

        p = mp.Process(target = interpolate_model, args = (shell,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

model_registration_time = time.time() - model_registration_start
with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- interpolation complete; runtime: {model_registration_time / 60} minutes / {(model_registration_time / 60) / 60} hours\n')

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- preparing final model\n')
# finally, concatenate all shells into the final registered model and compute RMS
df_registered_model = pd.DataFrame(columns = ['SHELL#', 'SHELL_MID', 'BLOCK#', 'BLOCK_LAT_MID', 'BLOCK_LON_MID', out_property_header])
for shell in all_shells:
    df_registered_shell = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell}.csv')
    df_registered_model = pd.concat([df_registered_model, df_registered_shell])
    os.remove(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell}.csv')

df_registered_model = df_registered_model[['SHELL#', 'BLOCK#', out_property_header]].copy()
df_registered_model['SHELL#'] = df_registered_model['SHELL#'].astype(int)
df_registered_model['BLOCK#'] = df_registered_model['BLOCK#'].astype(int)

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- computing the RMS profile\n')
## RMS PERTURBATION COMPUTATION ##
s = []
r = []
n = []

for shell in all_shells:
    df_slice = df_registered_model.loc[df_registered_model['SHELL#'] == shell]
    v = rms(list(df_slice[out_property_header]))
    s.append(shell)
    r.append(v)

max_value = max(r)

for val in r:
    w = val / max_value
    n.append(w)

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt6_log.txt', 'a') as fout:
    fout.write(f'- saving the registered model and the RMS output files\n')
# save the rms file
df_rms = pd.DataFrame(data = {'SHELL#': s, RMS_header: r, f'NORM_{RMS_header}': n})
df_rms.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/{model_name}_RMS.csv', index = False)

# save the interpolated file
df_registered_model.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/{model_name}_grid_registered.csv', index = False)
    














