import numpy as np
import pandas as pd
import mod_input
import mod_geo
import mod_database
import mod_refmodels
import netCDF4
import glob
import os


# define relevant input variables from mod_input, including file names and paths
cdp = mod_input.computed_decimal_places
input_csv = f'{mod_input.input_model}.csv'
output_csv = f'{mod_input.input_model}_grid_registered.csv'
df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
lat_header = mod_input.lat_header
lon_header = mod_input.lon_header
depth_header = mod_input.depth_header
prop_header = mod_input.property_header
if mod_input.data_wave_type == 'S':
    out_property_header = 'dVs_%'
    RMS_header = 'RMS_dVs'
elif mod_input.data_wave_type == 'P':
    out_property_header = 'dVp_%'
    RMS_header = 'RMS_dVp'


# define the input model to interpolate
df_input_model = pd.read_csv(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}', sep = mod_input.delimiter, comment = '#')

# compute the voigt average, if necessary
if mod_input.voigt[0] == True:
    df_input_model[prop_header] = np.sqrt(((2. * (df_input_model[mod_input.voigt[1]] ** 2)) + (df_input_model[mod_input.voigt[2]] ** 2)) / 3.)

df_model = df_input_model[[lat_header, lon_header, depth_header, prop_header]]

# condition the input model so there are not lat/lon/depth gaps.
# longitude
max_lon = max(df_model[lon_header])
if max_lon > 180.:
    df_model[lon_header] = df_model[lon_header] - 180.
max_lon = max(df_model[lon_header])
min_lon = min(df_model[lon_header])

max_lat = max(df_model[lat_header])
min_lat = min(df_model[lat_header])

min_depth = min(df_model[depth_header])
if min_depth > min(df_shells['DEPTH_MID']):
    df_depth_add = df_model.loc[df_model[depth_header] == min_depth].copy()
    df_depth_add.loc[:, depth_header] = 0.
    df_model = pd.concat([df_depth_add, df_model])

max_depth = max(df_model[depth_header])
if max_depth < max(df_shells['DEPTH_MID']):
    df_depth_add = df_model.loc[df_model[depth_header] == max_depth].copy()
    df_depth_add.loc[:, depth_header] = max(df_shells['DEPTH_MAX'])
    df_model = pd.concat([df_model, df_depth_add])
df_model.reset_index()

# make lists of block and shell ids
shell_ids = list(df_shells['SHELL#'])
block_ids = list(df_blocks['BLOCK#'])
lat_bands = list(df_blocks['LAT_BAND_ID'].unique())


# make lists of unique values of the three dimensions of the input model
model_depths = list(df_model[depth_header].unique())
model_lats = list(df_model[lat_header].unique())
model_lons = list(df_model[lon_header].unique())


# make dataframe for interpolated values
new_shells = []
new_blocks = []
new_props = []

for shell in shell_ids:
    for block in block_ids:
        new_shells.append(shell)
        new_blocks.append(block)
        new_props.append(0.)
df_interp = pd.DataFrame(data = {'SHELL#': new_shells, 'BLOCK#': new_blocks, out_property_header: new_props})


## MAKE DATAFRAMES FOR MODEL BLOCKS ##
# depth
model_depth_mins = []
model_depth_maxs = []
for depth in range(len(model_depths) - 1):
    depth_2 = depth + 1
    model_depth_mins.append(model_depths[depth])
    model_depth_maxs.append(model_depths[depth_2])
df_model_depths = pd.DataFrame(data = {'depth_min': model_depth_mins, 'depth_max': model_depth_maxs})

# latitude
model_lat_mins = []
model_lat_maxs = []
for lat in range(len(model_lats) - 1):
    lat_2 = lat + 1
    model_lat_mins.append(model_lats[lat])
    model_lat_maxs.append(model_lats[lat_2])
df_model_lats = pd.DataFrame(data = {'lat_min': model_lat_mins, 'lat_max': model_lat_maxs})

# longitude
model_lon_mins = []
model_lon_maxs = []
for lon in range(len(model_lons) - 1):
    lon_2 = lon + 1
    model_lon_mins.append(model_lons[lon])
    model_lon_maxs.append(model_lons[lon_2])
df_model_lons = pd.DataFrame(data = {'lon_min': model_lon_mins, 'lon_max': model_lon_maxs})

# make dataframe for coordinate system overlap
coord_block_ids = []
center_lats = []
center_lons = []
tomog_lat_mins = []
tomog_lat_maxs = []
tomog_lon_mins = []
tomog_lon_maxs = []

for band in lat_bands:
    blocks = df_blocks.loc[df_blocks['LAT_BAND_ID'] == band]['BLOCK#']
    for block in blocks:
        block_info = mod_database.get_block_info(block)
        center_lat = block_info[6]
        center_lon = block_info[7]
        coord_block_ids.append(block_info[0])
        center_lats.append(center_lat)
        center_lons.append(center_lon)
        
        # if the center lon is less than the minimum longitude in the model or greater than the maximum lon in the model:
        if (center_lon < min_lon) or (center_lon >= max_lon):
            lon_min = max_lon
            lon_max = min_lon
        else:
            df_tomog_lons = df_model_lons.loc[(df_model_lons['lon_min'] <= center_lon) & (df_model_lons['lon_max'] > center_lon)]
            lon_min = float(df_tomog_lons['lon_min'])
            lon_max = float(df_tomog_lons['lon_max'])
        tomog_lon_mins.append(lon_min)
        tomog_lon_maxs.append(lon_max)
        
        if center_lat < min_lat:
            lat_min = min_lat #-90.
            lat_max = min_lat
        elif center_lat >= max_lat:
            lat_min = max_lat
            lat_max = max_lat #90.
        else:
            df_tomog_lats = df_model_lats.loc[(df_model_lats['lat_min'] <= center_lat) & (df_model_lats['lat_max'] > center_lat)]
            lat_min = float(df_tomog_lats['lat_min'])
            lat_max = float(df_tomog_lats['lat_max'])
        tomog_lat_mins.append(lat_min)
        tomog_lat_maxs.append(lat_max)

df_lateral = pd.DataFrame(data = {'BLOCK#': coord_block_ids, 'CENTER_LAT': center_lats, 'CENTER_LON': center_lons, 'MODEL_LAT_MIN': tomog_lat_mins, 'MODEL_LAT_MAX': tomog_lat_maxs, 'MODEL_LON_MIN': tomog_lon_mins, 'MODEL_LON_MAX': tomog_lon_maxs})


# shell, center_depth, model_depth_min, model_depth_max
center_depths = []
top_depths = []
bottom_depths = []
for shell in shell_ids:
    shell_info = mod_database.get_shell_info(shell)
    center_depth = shell_info[2]
    top_depth = float(df_model_depths.loc[(df_model_depths['depth_min'] <= center_depth) & (df_model_depths['depth_max'] > center_depth)]['depth_min'])
    bottom_depth = float(df_model_depths.loc[(df_model_depths['depth_min'] <= center_depth) & (df_model_depths['depth_max'] > center_depth)]['depth_max'])

    center_depths.append(center_depth)
    top_depths.append(top_depth)
    bottom_depths.append(bottom_depth)
    
df_radial = pd.DataFrame(data = {'SHELL#': shell_ids, 'CENTER_DEPTH': center_depths, 'MODEL_DEPTH_MIN': top_depths, 'MODEL_DEPTH_MAX': bottom_depths})


## MAIN REREGISTRATION CODE ##
# for shell in shell_ids:
for s in range(len(df_radial)):
    shell = df_radial['SHELL#'].iloc[s]
    if shell == shell_ids[0]:
        for block in block_ids:
            element_idx = df_interp.loc[(df_interp['BLOCK#'] == block) & (df_interp['SHELL#'] == shell)].index
            df_interp.loc[element_idx, out_property_header] = 0.
    else:
        # grab the center depth from our depth shells. this is the target depth
        # find the values for the min and max depths from the model
        depth_mid = df_radial['CENTER_DEPTH'].iloc[s]
        model_depth_min = float(df_radial['MODEL_DEPTH_MIN'].iloc[s])
        model_depth_max = float(df_radial['MODEL_DEPTH_MAX'].iloc[s])
        df_model_ztop = df_model.loc[df_model[depth_header] == model_depth_min].copy()
        df_model_zbottom = df_model.loc[df_model[depth_header] == model_depth_max].copy()

        for block in range(len(df_lateral)):
            block_id = df_lateral['BLOCK#'].iloc[block]
            center_lat = float(df_lateral['CENTER_LAT'].iloc[block])
            center_lon = float(df_lateral['CENTER_LON'].iloc[block])
            lat_min = float(df_lateral['MODEL_LAT_MIN'].iloc[block])
            lat_max = float(df_lateral['MODEL_LAT_MAX'].iloc[block])
            lon_min = float(df_lateral['MODEL_LON_MIN'].iloc[block])
            lon_max = float(df_lateral['MODEL_LON_MAX'].iloc[block])

            # find the velocity values at all eight vertices of the model block and interpolate
            # top surface/depth min:
            v1 = float(df_model_ztop.loc[(df_model_ztop[lat_header] == lat_min) & (df_model_ztop[lon_header] == lon_min)][prop_header])
            v2 = float(df_model_ztop.loc[(df_model_ztop[lat_header] == lat_min) & (df_model_ztop[lon_header] == lon_max)][prop_header])
            v3 = float(df_model_ztop.loc[(df_model_ztop[lat_header] == lat_max) & (df_model_ztop[lon_header] == lon_min)][prop_header])
            v4 = float(df_model_ztop.loc[(df_model_ztop[lat_header] == lat_max) & (df_model_ztop[lon_header] == lon_max)][prop_header])

            # interpolate to target coordinates
            i1 = np.interp(center_lon, [lon_min, lon_max], [v1, v2])
            i2 = np.interp(center_lon, [lon_min, lon_max], [v3, v4])
            i3 = np.interp(center_lat, [lat_min, lat_max], [i1, i2])

            # bottom surface/depth_max
            v5 = float(df_model_zbottom.loc[(df_model_zbottom[lat_header] == lat_min) & (df_model_zbottom[lon_header] == lon_min)][prop_header])
            v6 = float(df_model_zbottom.loc[(df_model_zbottom[lat_header] == lat_min) & (df_model_zbottom[lon_header] == lon_max)][prop_header])
            v7 = float(df_model_zbottom.loc[(df_model_zbottom[lat_header] == lat_max) & (df_model_zbottom[lon_header] == lon_min)][prop_header])
            v8 = float(df_model_zbottom.loc[(df_model_zbottom[lat_header] == lat_max) & (df_model_zbottom[lon_header] == lon_max)][prop_header])

            # interpolate to target coordinates
            i4 = np.interp(center_lon, [lon_min, lon_max], [v5, v6])
            i5 = np.interp(center_lon, [lon_min, lon_max], [v7, v8])
            i6 = np.interp(center_lat, [lat_min, lat_max], [i4, i5])

            # interpolate to get the velocity at the target depth
            i7 = np.interp(depth_mid, [model_depth_min, model_depth_max], [i3, i6])
                
            # append values to new lists
            element_idx = df_interp.loc[(df_interp['BLOCK#'] == block_id) & (df_interp['SHELL#'] == shell)].index
            df_interp.loc[element_idx, out_property_header] = i7

shells = df_interp['SHELL#'].unique()
blocks = df_interp['BLOCK#'].unique()
mid_depths = list(df_shells['DEPTH_MID'])

# convert values to % perturb. relative to the reference model
for shell in shells[1:]:
    df_update_slice = df_interp.loc[df_interp['SHELL#'] == shell].copy()
    update_mid_depth = float(df_shells.loc[df_shells['SHELL#'] == shell]['DEPTH_MID'])
    update_ref_val = mod_refmodels.prem_vel(mod_input.data_wave_type, update_mid_depth)

    for b in range(len(df_update_slice)):
        block = df_update_slice['BLOCK#'].iloc[b]
        vel = df_update_slice[out_property_header].iloc[b]
        updated_perturbation = ((vel / update_ref_val) * 100.) - 100.
        
        element_idx = df_interp.loc[(df_interp['BLOCK#'] == block) & (df_interp['SHELL#'] == shell)].index
        df_interp.loc[element_idx, out_property_header] = updated_perturbation


## RMS PERTURBATION COMPUTATION ##
def rms(property_list):
    n = len(property_list)
    x = 0
    for p in property_list:
        x += (p**2)
    return np.sqrt(x / n)

s = []
r = []
n = []

for shell in shells:
    df_slice = df_interp.loc[df_interp['SHELL#'] == shell]
    v = rms(list(df_slice[out_property_header]))
    s.append(shell)
    r.append(v)

max_value = max(r)

for val in r:
    w = val / max_value
    n.append(w)

# save the rms file
df_rms = pd.DataFrame(data = {'SHELL#': s, RMS_header: r, f'NORM_{RMS_header}': n})
np.savetxt(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{mod_input.input_model}_RMS.csv', df_rms, fmt = f'%i,%1.{cdp}f,%1.{cdp}f', delimiter = ',', header = f'SHELL#,{RMS_header},NORM_{RMS_header}', comments = '')


# save the interpolated file
col_names = f'SHELL#,BLOCK#,{out_property_header}'
np.savetxt(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{output_csv}', df_interp, fmt = f'%i,%i,%1.{cdp}f', delimiter = ',', header = col_names, comments = '')


