import os
import sys
import glob
import time
import shutil
import getopt
import mod_geo
import netCDF4
import mod_track
import mod_input
import numpy as np
import pandas as pd
import mod_database
import mod_refmodels
from netCDF4 import Dataset
import multiprocessing as mp
from datetime import datetime, timezone

pid = os.getpid()
start_time = time.time()
start_subj = f'Process began (PID: {pid}); part3.model_processing.py'
start_text = f'Process {pid} began\nModel: {mod_input.input_model}\nConvert from NC to CSV: {mod_input.convert_nc_to_csv}'
try:
    mod_track.SendMsg(start_subj, start_text)
except:
    pass
    
def cd(path):
    os.chdir(os.path.expanduser(path))

if mod_input.convert_nc_to_csv == True:
    
    models_name = mod_input.tomography_model_directory
    input_csv = f'{mod_input.input_model}.csv'
    input_nc = f'{mod_input.input_model}.nc'
    
    if os.path.exists(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}') == True:
        pass
    else:
        # make the .csv
        os.mkdir(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/')
        cd(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}')
        
        SCRIPT = os.path.basename(sys.argv[0])
        VERSION = 'V.2020.273'
    
        DEBUG = False
    
        GEOCSV_VERSION = 'GeoCSV2.0'
    
        # root _directory
        ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    
        # optional directory to store netCDF and GeoCSV files
        DATA_DIR = os.path.join(ROOT_DIR, models_name)
    
        
        LON_VARIABLE = 'longitude'  # must match the netCDF file's longitude variable
        LAT_VARIABLE = 'latitude'  # must match the netCDF file's longitude variable
        DEPTH_VARIABLE = 'depth'  # must match the netCDF file's depth variable
        VALID_MODES = {'depth': 'km', 'single': ''}
        DELIMITER = '|'
        NETCDF_FILE_NAME = input_nc #None
        BASE_NAME = None
        VIEW_HEADER = False
    
        # output mode (depth | lat | lon) as individual files based on depth, lat, lon
        # or (single) as a single file
        OUTPUT_MODE = 'single'
    
        def get_variable_attributes(model_data, header, variable, variable_name):
            """add  variable attributes to the header
    
            Keyword arguments:
            model_data: Dataset instance of the model_file
            header: GeoCSV header variable for the model
            variable : variable to add
            variable_name: variable name used to represent this variable
    
            Return values:
            the geoCSV header for the model
            """
            header.append('# {}_column: {}\n'.format(variable, variable_name))
            for attr, value in vars(model_data.variables[variable]).items():
                if '_range' in attr:
                    header.append('# {}_{}: {},{}\n'.format(variable, attr, value[0], value[1]))
                else:
                    header.append('# {}_{}: {}\n'.format(variable, attr, value))
            return header
    
    
        def get_model_header(model_file, model_data):
            """create GeoCSV header for the model
    
            Keyword arguments:
            model_file: the netCDF model file name
            model_data: Dataset instance of the model_file
    
            Return values:
            the geoCSV header for the model
            """
            header = list()
            # GeoCSV header
            header.append('# dataset: {}\n'.format(GEOCSV_VERSION))
            header.append('# created: {} UTC ({})\n'.format(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), SCRIPT))
            header.append('# netCDF_file: {}\n'.format(os.path.basename(model_file)))
            header.append('# delimiter: {}\n'.format(DELIMITER))
    
            # global attributes
            history_done = False
            history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")} Converted to GeoCSV by {SCRIPT} ,' \
                f'{VERSION} ' \
                f'from {NETCDF_FILE_NAME}'
            for attr, value in vars(model_data).items():
                if isinstance(value, str):
                    value = value.replace('\n', '; ')
                if attr.lower() == 'history':
                    value = f'{history}; {value}'
                    history_done = True
                header.append(f'# global_{attr}: {value}\n')
    
            if not history_done:
                header.append(f'# global_history: {history}\n')
    
            # variable s
            header = get_variable_attributes(model_data, header, LAT_VARIABLE, 'latitude')#, LAT_VARIABLE)
            header = get_variable_attributes(model_data, header, LON_VARIABLE, 'longitude')#, LON_VARIABLE)
            header = get_variable_attributes(model_data, header, DEPTH_VARIABLE, 'depth')#, DEPTH_VARIABLE)
    
            return ''.join(header)
    
    
        def get_var_header(model_data, var):
            """create GeoCSV header for a variable
    
            Keyword arguments:
            model_data: Dataset instance of the model_file
            var: the netCDF model file variable
    
            Return values:
            the geoCSV header for the variable
            """
            header = list()
            header = get_variable_attributes(model_data, header, var, var)
            return ''.join(header)
    
        def check_netcdf_file():
            """check the input netCDF model file and make sure it exist and extract info
    
            Return values:
            this_file:  file name with no extension
            """
            # check the model file and extract necessary information
            # must be in the argument list
            if NETCDF_FILE_NAME is None:
                print('[ERR] the netCDF model file name is required', flush=True)
                # usage_csv()
                sys.exit(1)
    
            # user may provide full path
            elif os.path.isfile(NETCDF_FILE_NAME):
                model_file_name = NETCDF_FILE_NAME
                base_file_name, ext = os.path.splitext(model_file_name)
    
            # user may place it under the data directory
            elif os.path.isfile(os.path.join(DATA_DIR, NETCDF_FILE_NAME)):
                model_file_name = os.path.join(DATA_DIR, NETCDF_FILE_NAME)
                base_file_name, ext = os.path.splitext(model_file_name)
    
            # could not find the file
            else:
                print('[ERR] could not find the netCDF model file {}'.format(NETCDF_FILE_NAME), flush=True)
                # usage_csv()
                sys.exit(1)
    
            return model_file_name, base_file_name
    
        def make_model_geocsv():
            """create GeoCSV file from a netCDF model file
    
            Keyword arguments:
            model_file: the netCDF model file name
            """
            model_file, base_file_name = check_netcdf_file()
    
            data_header = list()
            model_data = Dataset(model_file)
    
            try:
                # conversion to string is done to preserve precision
                lat = list()
                lon = list()
                depth = list()
                for this_value in model_data.variables[LAT_VARIABLE][:]:
                    lat.append("{}".format(str(this_value)))
                for this_value in model_data.variables[LON_VARIABLE][:]:
                    lon.append("{}".format(str(this_value)))
                for this_value in model_data.variables[DEPTH_VARIABLE][:]:
                    depth.append("{}".format(str(this_value)))
            except Exception as ex:
                print('\n[ERR] the expected variables ({}, {}, {}) not in the variable list: {}\n'.format(
                    LAT_VARIABLE, LON_VARIABLE, DEPTH_VARIABLE, str(list(model_data.variables.keys()))))
                sys.exit(1)
    
            emcin = {}
    
            # make sure this is a 3D netCDF file
            var_3d = list()
            for var in model_data.variables.keys():
                if len(model_data.variables[var].shape) == 3:
                    var_3d.append(var)
            if len(var_3d) <= 0:
                print('\n[ERR] not a 3D netCDF file\n\n', flush=True)
                sys.exit(1)
    
            output_data = list()
            for k, this_depth in enumerate(depth):
                if OUTPUT_MODE == 'single' and k == 0:
                    data_header = list()
                    output_file = '{}.csv'.format(base_file_name)
                    fp = open(output_file, 'w')
                    # print('[INFO] Output file: {}'.format(output_file), flush=True)
                    fp.write(get_model_header(model_file, model_data))
                    data_header.append('{}{}{}{}{}'.format(LAT_VARIABLE, DELIMITER, LON_VARIABLE,
                                                           DELIMITER, DEPTH_VARIABLE))
    
                if DEBUG:
                    print('[INFO] Processing depth: {}'.format(this_depth), flush=True)
    
                index = [-1, -1, -1]
                for i, this_lat in enumerate(lat):
                    for j, this_lon in enumerate(lon):
                        if OUTPUT_MODE == 'single':
                            output_data.append('{}{}{}{}{}'.format(str(this_lat), DELIMITER, str(this_lon), DELIMITER,
                                                                   str(this_depth)))
                        for var in model_data.variables.keys():
                            depth_index = None
                            lat_index = None
                            lon_index = None
                            if var.encode('ascii', 'ignore').decode("utf-8") not in [LAT_VARIABLE, LON_VARIABLE,
                                                                                     DEPTH_VARIABLE]:
                                if ((OUTPUT_MODE == 'single' and (not i and not j and not k)) or
                                        (OUTPUT_MODE == 'depth' and (not i and not j))):
                                    fp.write(get_var_header(model_data, var))
                                    data_header.append('{}{}'.format(DELIMITER, var))
                                # find the variable ordering
                                if lat_index is None:
                                    for l in range(len(model_data.variables[var].dimensions)):
                                        if model_data.variables[var].dimensions[l].encode('ascii', 'ignore').decode(
                                               "utf-8") == DEPTH_VARIABLE:
                                            depth_index = l
                                        elif model_data.variables[var].dimensions[l].encode('ascii', 'ignore').decode(
                                                "utf-8") == LON_VARIABLE:
                                            lon_index = l
                                        else:
                                            lat_index = l
    
                                    if var not in emcin.keys():
                                        try:
                                            emcin[var] = model_data.variables[var][:]
                                        except Exception as err:
                                            print('\n[ERR] problem reading variable "{}"'.format(var))
                                            print('{0}\n'.format(err))
                                            sys.exit(2)
    
                                index[depth_index] = k
                                index[lat_index] = i
                                index[lon_index] = j
                                # nan values, otherwise we write string to preserve the precision
                                if str(emcin[var][index[0]][index[1]][index[2]]) == '--':
                                    output_data.append('{}{}'.format(DELIMITER,
                                                                     float(emcin[var][index[0]][index[1]][index[2]])))
                                else:
                                    # conversion to string is done to preserve precision
                                    output_data.append('{}{}'.format(DELIMITER,
                                                                     str(emcin[var][index[0]][index[1]][index[2]])))
                        output_data.append('\n')
    
            if OUTPUT_MODE == 'single':
                fp.write('{}\n'.format(''.join(data_header)))
                fp.write(''.join(output_data))
                fp.close()
        make_model_geocsv()
    
        cd('../../')
        shutil.move(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{input_csv}', f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}')
    
        # now that the model is converted, remove all of the comments at the beginning of the file
        df_model = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}', sep = '|', comment = '#')
        df_model.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}', index = False)


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
shells_to_register = all_shells.loc[~all_shells.isin(mod_input.zero_shells)]


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

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'w') as fout:
    fout.write(f'Starting grid registration for model {model_name}\n')

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

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
    fout.write(f'- finding all shell and block IDs for model nodes\n')

df_model['SHELL#'] = 0
df_model['BLOCK#'] = 0
df_model = np.array(df_model)
for line in range(len(df_model)):
    model_lat = df_model[line, 0]
    model_lon = df_model[line, 1]
    model_depth = df_model[line, 2]
    if model_depth < 0:
        pass
    else:
        model_shell = mod_database.find_shell_id(model_depth)
        model_block = mod_database.find_block_id(model_lat, model_lon)
        df_model[line, 4] = model_shell
        df_model[line, 5] = model_block

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
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

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
    fout.write(f'- interpolating the model\n')

def interpolate_model(shell_no):
    # reregister the model to our grid system by either finding the value at the center of our blocks or making an average of all model values that fall within each block
    df_registered_shell = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell_no}.csv')
    shell_mid = float(df_registered_shell['SHELL_MID'].iloc[0])
    df_registered_shell = np.array(df_registered_shell)
    for line in range(len(df_registered_shell)):
        block = df_registered_shell[line, 2]
        block_lat_mid = df_registered_shell[line, 3]
        block_lon_mid = df_registered_shell[line, 4]

        # first, check to see if there's a model value at this exact point.
        # model_indices = np.where((df_model.T[0] == block_lat_mid) & (df_model.T[1] == block_lon_mid) & (df_model.T[2] == shell_mid))[0]
        # if len(model_indices) == 1:
        #     df_registered_shell[line, 5] = df_model[model_indices][0][3]
        # # if there isn't a model value at this exact point, calculate an average of all model values that fall into our block
        # else:
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
                updated_perturbation = 0.
            else:
                updated_perturbation = ((vel / update_ref_val) * 100.) - 100.
            df_registered_shell[idx, 5] = updated_perturbation

    df_registered_shell = pd.DataFrame(df_registered_shell, columns = ['SHELL#', 'SHELL_MID', 'BLOCK#', 'BLOCK_LAT_MID', 'BLOCK_LON_MID', out_property_header])
    df_registered_shell.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/tmp_shell_{shell}.csv', index = False)


model_registration_start = time.time()
if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for shell in shells_to_register:
        process_idx += 1

        p = mp.Process(target = interpolate_model, args = (shell,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

model_registration_time = time.time() - model_registration_start


with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
    fout.write(f'- interpolation complete; runtime: {mod_track.runtime(model_registration_time)}\n')

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
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

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
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

with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
    fout.write(f'- saving the registered model and the RMS output files\n')

# save the rms file
df_rms = pd.DataFrame(data = {'SHELL#': s, RMS_header: r, f'NORM_{RMS_header}': n})
df_rms.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/{model_name}_RMS.csv', index = False)

# save the interpolated file
df_registered_model.to_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/{model_name}_grid_registered.csv', index = False)





with open(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/pt3_log_model_processing.txt', 'a') as fout:
    fout.write(f'- making the original model plot files\n')
    
try:
    os.mkdir(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/original_model_plot_files')
except:
    pass

lats = np.arange(mod_input.start_lat, mod_input.final_lat, mod_input.reference_lat)
lons = np.arange(mod_input.start_lon, mod_input.final_lon, mod_input.reference_lon)
plot_file_path = f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{model_name}_update/original_model_plot_files'

for shell_to_plot in all_shells:
    og_perturbs_file = f'{plot_file_path}/{model_name}_shell_{shell_to_plot}_original_perturbs_plot_ready_{mod_input.reference_lat}deg_lat_by_{mod_input.reference_lon}deg_lon.csv'
    df_og_shell_data = df_registered_model.loc[df_registered_model['SHELL#'] == shell_to_plot]
    df_og = pd.DataFrame(data = {'LON': lons})

    # make empty lists that will be filled with numpy arrays for each latitude band, for each model.
    og_dvs = []
    
    # loop through all latitudes in the model space
    for lat in lats:
        # make empty lists to fill in the perturbations for the current latitude band
        lat_og_dvs = []
    
        # loop through all longitudes in the model space
        for lon in lons:
            # find the block that the current lat/lon pair falls in
            block = mod_database.find_block_id(lat, lon)
            original_perturb = float(df_og_shell_data.loc[df_og_shell_data['BLOCK#'] == block][out_property_header].iloc[0])
            lat_og_dvs.append(original_perturb)
        
        og_dvs.append(np.array(lat_og_dvs))
        df_og[f'{lat}'] = lat_og_dvs
    df_og.to_csv(og_perturbs_file, index = False)

end_subj = f'Process ended (PID: {pid}); part3.model_processing.py'
end_text = f'Process {pid} complete;\nRuntime: {mod_track.runtime(time.time() - start_time)}'
try:
    mod_track.SendMsg(end_subj, end_text)
except:
    pass
