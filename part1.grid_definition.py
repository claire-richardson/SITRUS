# import libraries and packages:
import numpy as np
import pandas as pd
import mod_input
import mod_geo
import os

# define starting block parameters and output files:
depth_bounds = mod_input.shell_bounds
depth_mins = depth_bounds[:-1]
depth_maxs = depth_bounds[1:]
total_radius = mod_input.total_radius
shell_numbers = list(range(1, len(depth_bounds)))

ref_lat = mod_input.reference_lat
ref_lon = mod_input.reference_lon
start_lat = mod_input.start_lat
final_lat = mod_input.final_lat
start_lon = mod_input.start_lon
final_lon = mod_input.final_lon
total_lat = final_lat - start_lat
total_lon = final_lon - start_lon
lat_steps = int(total_lat / ref_lat)
lats = list(range(start_lat, final_lat + ref_lat, ref_lat))

shell_file = mod_input.shell_file
block_file = mod_input.block_file
cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places

#################
# DEPTH SHELLS: #
#################
rad_data = {'SHELL#': shell_numbers, 'DEPTH_MIN': depth_mins, 'DEPTH_MAX': depth_maxs}
df_rad = pd.DataFrame(data = rad_data)
df_rad['DEPTH_DIFF'] = df_rad['DEPTH_MAX'] - df_rad['DEPTH_MIN']
df_rad['DEPTH_DIFF_MID'] = df_rad['DEPTH_DIFF'] / 2
df_rad['DEPTH_CENTER'] = df_rad['DEPTH_MIN'] + df_rad['DEPTH_DIFF_MID']
df_rad['TOP_RADIUS'] = total_radius - df_rad['DEPTH_MIN']
df_rad['MID_RADIUS'] = total_radius - df_rad['DEPTH_CENTER']
df_rad['BOTTOM_RADIUS'] = total_radius - df_rad['DEPTH_MAX']
df_rad['DEPTH_MID'] = df_rad['DEPTH_MIN'] + df_rad['DEPTH_DIFF_MID']


# reconfigure data and save as csv.
df_shell = pd.DataFrame()
df_shell['SHELL#'] = df_rad['SHELL#']
df_shell['DEPTH_MIN'] = df_rad['DEPTH_MIN']
df_shell['DEPTH_MID'] = df_rad['DEPTH_MID']
df_shell['DEPTH_MAX'] = df_rad['DEPTH_MAX']
df_shell['DEPTH_DIFF'] = df_rad['DEPTH_DIFF']
df_shell['TOP_RADIUS'] = df_rad['TOP_RADIUS']
df_shell['MID_RADIUS'] = df_rad['MID_RADIUS']
df_shell['BOTTOM_RADIUS'] = df_rad['BOTTOM_RADIUS']

np.savetxt(shell_file, df_shell, fmt = f'%i,%1.{rdp}f,%1.{rdp}f,%1.{rdp}f,%1.{rdp}f,%1.{rdp}f,%1.{rdp}f,%1.{rdp}f', delimiter = ',', header = 'SHELL#,DEPTH_MIN,DEPTH_MID,DEPTH_MAX,DEPTH_DIFF,TOP_RADIUS,MID_RADIUS,BOTTOM_RADIUS', comments = '')


##############
# 2D BLOCKS: #
##############
# define data lists
radius = []
equatorial_circumference = []
degree_in_km = []
block_area = []
circumference_per_lat = []
latitude = []

# make calculations for approximate equal area blocks:
for depth in depth_mins:
    # radius:
    rad = total_radius - depth
    radius.append(rad)

    # gridblock area:
    area = (np.pi/180) * rad**2 * (np.sin(np.deg2rad(ref_lat)) - np.sin(np.deg2rad(0))) * (ref_lon - 0)
    block_area.append(area)

    # circumference at equator:
    eq_circ = 2 * np.pi * rad
    equatorial_circumference.append(eq_circ)

    # degrees of lon in km at the equator:
    deg_in_km = eq_circ / total_lon
    degree_in_km.append(deg_in_km)

    # circumference at other latitudes:
    for lat in range(start_lat, final_lat, ref_lat):
        latitude.append(lat)
        absolute_lat = np.absolute(lat)
        if lat != 0:
            scale_lat = 90 - absolute_lat
        if lat == 0:
            scale_lat = 90
        scale_factor = scale_lat / 90
        circ = scale_factor * eq_circ
        circumference_per_lat.append(circ)

    # block info:
    block_data = {'LATITUDE': latitude, 'CIRC': circumference_per_lat}
    df_block = pd.DataFrame(data = block_data)
    df_block['LATITUDE'] = df_block['LATITUDE'] + 360
    df_block['BAND_AREA'] = (2 * np.pi * rad**2) * (np.sin(np.deg2rad(df_block['LATITUDE'] + ref_lat)) - np.sin(np.deg2rad(df_block['LATITUDE'])))
    df_block['LATITUDE'] = df_block['LATITUDE'] - 360
    df_block['#_OF_BLOCKS'] = df_block['BAND_AREA'] / area
    df_block['#_OF_ROUND_BLOCKS'] = df_block['#_OF_BLOCKS'].round(decimals = 0)
    df_block['#_OF_ROUND_BLOCKS'] = df_block['#_OF_ROUND_BLOCKS'].astype(int)
    df_block['BLOCK_AREA'] = df_block['BAND_AREA'] / df_block['#_OF_ROUND_BLOCKS']
    
# depth info:
depth_data = {'RADIUS': radius, 'EQ_CIRC': equatorial_circumference, 'DEG_IN_KM': degree_in_km, 'BLOCK_AREA': block_area}
df_depth = pd.DataFrame(data=depth_data)

# block file output:
block_number = []
lat_band_id = []
lat_min = []
lat_max = []
lon_min = []
lon_max = []
center_lat = []
center_lon = []

df_tmp = pd.DataFrame()
df_tmp['#_OF_ROUND_BLOCKS'] = df_block['#_OF_ROUND_BLOCKS'].head(n = lat_steps)
df_tmp['LON_DIVISIONS'] = total_lon / df_tmp['#_OF_ROUND_BLOCKS']
lon_divisions = df_tmp['LON_DIVISIONS'].to_list()
blocks = df_tmp['#_OF_ROUND_BLOCKS'].to_list()

# initialize lat band id
lat_band = 1
# lat min, max, center: 
for value in blocks:
    for num in range(value):
        lat_min.append(start_lat)
        lat_max.append(start_lat + ref_lat)
        center_lat.append(start_lat + (ref_lat / 2))
        lat_band_id.append(lat_band)
    lat_band += 1
    start_lat += ref_lat

start_lat = mod_input.start_lat

# lon min, max, center:
for value in blocks:
    division = total_lon / value
    lon = start_lon
    for num in range(value):
        lon_min.append(lon)
        lon_max.append(lon + division)
        center_lon.append(lon + (division / 2))
        lon += division
        
# block number:
for value in range(len(lat_min)):
    block_no = value + 1
    block_number.append(block_no)

block_dim_data = {'BLOCK#': block_number, 'LAT_BAND_ID': lat_band_id, 'LAT_MIN': lat_min, 'LAT_MAX': lat_max, 'LON_MIN': lon_min, 'LON_MAX': lon_max, 'CENTER_LAT': center_lat, 'CENTER_LON': center_lon}
df_block_dim = pd.DataFrame(data = block_dim_data)
np.savetxt(block_file, df_block_dim, fmt = f'%i, %i, %1.{cdp}f, %1.{cdp}f, %1.{cdp}f, %1.{cdp}f, %1.{cdp}f, %1.{cdp}f', delimiter = ',', header = 'BLOCK#,LAT_BAND_ID,LAT_MIN,LAT_MAX,LON_MIN,LON_MAX,CENTER_LAT,CENTER_LON', comments = '')

## NEAR NEIGHBORS
if mod_input.make_near_neighbors_files == True:
    os.mkdir(f'./{mod_input.near_neighbors_directory}')

    for b in range(len(df_block_dim)):
        b_no = b + 1
        print(f'working on block number {b_no} of {len(df_block_dim)}')
        b_center_lat = df_block_dim.loc[df_block_dim['BLOCK#'] == b_no]['CENTER_LAT']
        b_center_lon = df_block_dim.loc[df_block_dim['BLOCK#'] == b_no]['CENTER_LON']
        eligible_blocks = []
        radius_degrees = []

        for block in range(len(df_block_dim)):
            block_no = block + 1
            block_center_lat = df_block_dim.loc[df_block_dim['BLOCK#'] == block_no]['CENTER_LAT']
            block_center_lon = df_block_dim.loc[df_block_dim['BLOCK#'] == block_no]['CENTER_LON']
            distance_to_center = mod_geo.GCP_length(b_center_lat, b_center_lon, block_center_lat, block_center_lon)

            if distance_to_center <= mod_input.max_near_neighbors_radius:
                eligible_blocks.append(block_no)
                radius_degrees.append(distance_to_center)

        df_near_neighbors = pd.DataFrame(data = {'NEIGHBOR': eligible_blocks, 'RADIUS_DEG': radius_degrees})
        np.savetxt(f'{mod_input.near_neighbors_directory}/block_{b_no}_neighbors.csv', df_near_neighbors, fmt = f'%i,%1.{cdp}f', delimiter = ',', header = 'NEIGHBOR,RADIUS_DEG', comments = '')
