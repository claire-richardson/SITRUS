##########
#
# Pandas functions to index the 3D mesh in SITRUS updates
#
# READ ME:
#
# This module contains various functions to index the 3D mesh defined by the first script in
# the SITRUS package. These functions depend on the spatial parameters defined in the 
# `shell_dimensions.csv` and `block_dimensions.csv` files, which are made by the first SITRUS
# script. This module calls the depth of the core-mantle boundary from `mod_input.py, which
# should be equal to the value at the base of the lowermost depth shell in
# `shell_dimensions.csv`. The user should make sure that that value is correct/consistent.
#
# the functions in this module include:
# 1.  get_shell_info: get ID and spatial dimensions of a given shell ID
# 2.  find shell_id: find the shell ID of a given depth
# 3.  slice_shell_top: find the shell ID of a depth shell with a given top depth bound
# 4.  slice_shell_bottom: find the shell ID of depth shell with a given bottom depth bound
# 5.  get_block_info: get IDs and spatial dimensions of a given block ID
# 6.  find_block_id: find the block ID of a given lat/lon pair
# 7.  find_boundary_type: find the type of boundary in the 3D mesh that a given point on a raypath falls on
# 8.  slice_block_mins: find the block ID of a block from its minimum latitude and longitude bounds
# 9.  slice_block_min_max: find the block ID of a block from its minimum latitude bound and maximum longitude bound
# 10. slice_block_min_mid: find the block ID of a block from its minimum latitude bound and a longitude value between its minimum and maximum longitude bounds
# 11. slice_block_max_mid: find the block ID of a block from its maximum latitude bound and a longitude value between its minimum and maximum longitude bounds
# 12. slice_block_band_lon: find the block ID of a block from its latitude band ID and a longitude value between its minimum and maximum longitude bounds
# 
##########

import mod_input
import pandas as pd

CMB_depth = mod_input.core_mantle_boundary
shells = mod_input.shell_file
blocks = mod_input.block_file
df_shells = pd.read_csv(shells)
df_blocks = pd.read_csv(blocks)

#### SHELL INFORMATION ####
def get_shell_info(shell_id):
    '''
    Function to find the dimensional information of the given depth shell. Calls information from `shells` (default: `shell_dimensions.csv`)
    ======
    Inputs:
    - `shell_id`: Shell ID number of the shell of interest ('SHELL#' in `shell_dimensions.csv')
    ======
    Outputs:
    - list of 8 values describing the shell dimensions and attributes, including:
        - shell_id: Shell ID number [format: int; unit: unitless]
        - depth_min: Minimum depth of the shell, i.e., the depth value at the top of the shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - depth_mid: Middle depth of the shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - depth_max: Maximum depth of the shell, i.e., the depth value at the base of the shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - depth_diff: Thickness of the shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - top_radius: Radius of the top of the depth shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - mid_radius: Radius of the middle of the depth shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
        - bottom_radius: Radius of the base of the depth shell [format: float; unit: km or same as convention in `shell_dimensions.csv`]
    ======
    '''
    df_slice = df_shells.loc[df_shells['SHELL#'] == shell_id]
    shell_id = int(df_slice['SHELL#'])
    depth_min = float(df_slice['DEPTH_MIN'])
    depth_mid = float(df_slice['DEPTH_MID'])
    depth_max = float(df_slice['DEPTH_MAX'])
    depth_diff = float(df_slice['DEPTH_DIFF'])
    top_radius = float(df_slice['TOP_RADIUS'])
    mid_radius = float(df_slice['MID_RADIUS'])
    bottom_radius = float(df_slice['BOTTOM_RADIUS'])
    return [shell_id, depth_min, depth_mid, depth_max, depth_diff, top_radius, mid_radius, bottom_radius]

def find_shell_id(z):
    '''
    Function to find the shell ID that a given depth belongs to
    ======
    Inputs:
    - `z`: depth [format: float or int; unit: km or same as convention in `shell_dimensions.csv`]
    Outputs:
    - shell ID of depth `z` [format: int; unit: unitless]
    ======
    '''
    if z == CMB_depth:
        df_slice = df_shells.loc[df_shells['DEPTH_MAX'] == z]
    else:
        df_slice = df_shells.loc[(df_shells['DEPTH_MIN'] <= z) & (df_shells['DEPTH_MAX'] > z)]
    shell_id = int(df_slice['SHELL#'])
    return shell_id

def slice_shell_top(top_depth):
    '''
    Function to find a shell ID by its top depth
    ======
    Inputs:
    - `top_depth`: the top depth of the depth shell you wish to find [format: float or int; unit: km or same as convention in `shell_dimensions.csv`]
    ======
    Outputs:
    - shell ID with top depth ('DEPTH_MIN') equal to `top_depth` [format: int; unit: unitless]
    '''
    df_slice = df_shells.loc[df_shells['DEPTH_MIN'] == top_depth]
    shell_id = int(df_slice['SHELL#'])
    return shell_id

def slice_shell_bottom(bottom_depth):
    '''
    Function to find a shell ID by its bottom depth
    ======
    Inputs:
    - `bottom_depth`: the bottom depth of the depth shell you wish to find [format: float or int; unit: km or same as convention in `shell_dimensions.csv`]
    ======
    Outputs:
    - shell ID with bottom depth ('DEPTH_MAX') equal to `bottom_depth` [format: int; unit: unitless]
    '''
    df_slice = df_shells.loc[df_shells['DEPTH_MAX'] == bottom_depth]
    shell_id = int(df_slice['SHELL#'])
    return shell_id

#### BLOCK INFORMATION ####
def get_block_info(block_id):
    '''
    Function to find the dimensional information of the given block. Calls information from `blocks` (default: `block_dimensions.csv`)
    ======
    Inputs:
    - `block_id`: Block ID number of the block of interest ('BLOCK#' in `block_dimensions.csv')
    ======
    Outputs:
    - list of 8 values describing the block dimensions and attributes, including:
        - block_id: Block ID number [format: int; unit: unitless]
        - lat_id: ID number of the latitude band of the block [format: int; unit: unitless]
        - lat_min: Minimum latitude of the block, i.e., the latitude value at the southern bound of the block [format: float; unit: degrees]
        - lat_max: Maximum latitude of the block, i.e., the latitude value at the northern bound of the block [format: float; unit: degrees]
        - lon_min: Minimum longitude of the block, i.e., the longitude value at the western bound of the block [format: float; unit: degrees]
        - lon_max: Maximum longitude of the block, i.e., the longitude value at the eastern bound of the block [format: float; unit: degrees]
        - center_lat: Latitude at the center of the block [format: float; unit: degrees]
        - center_lon: Longitude at the center of the block [format: float; unit: degrees]
    ======
    '''
    df_slice = df_blocks.loc[df_blocks['BLOCK#'] == block_id]
    block_id = int(df_slice['BLOCK#'])
    lat_id = int(df_slice['LAT_BAND_ID'])
    lat_min = float(df_slice['LAT_MIN'])
    lat_max = float(df_slice['LAT_MAX'])
    lon_min = float(df_slice['LON_MIN'])
    lon_max = float(df_slice['LON_MAX'])
    center_lat = float(df_slice['CENTER_LAT'])
    center_lon = float(df_slice['CENTER_LON'])
    return [block_id, lat_id, lat_min, lat_max, lon_min, lon_max, center_lat, center_lon]

def find_block_id(lat, lon):
    '''
    Function to find the block ID that a given latitude and longitude pair belong to
    ======
    Inputs:
    - `lat`: latitude [format: float or int; unit: degrees]
    - `lon`: longitude [format: float or int; unit: degrees]
    Outputs:
    - shell ID of depth `z` [format: int; unit: unitless]
    ======
    '''
    # if lon >= 180.:
    #     diff = lon - 180.
    #     lon = 180. - diff
    # elif lon < -180.:
    #     diff = -180. - lon
    #     lon = 180. - diff
    if lat != 90.:
        df_slice = df_blocks.loc[(df_blocks['LAT_MIN'] <= lat) & (df_blocks['LAT_MAX'] > lat) & (df_blocks['LON_MIN'] <= lon) & (df_blocks['LON_MAX'] > lon)]
    elif lat == 90.:
        df_slice = df_blocks.loc[(df_blocks['LAT_MAX'] == 90.) & (df_blocks['LON_MIN'] <= lon) & (df_blocks['LON_MAX'] > lon)]
    block_id = int(df_slice['BLOCK#'])
    return block_id

def find_boundary_type(lat, lon, z):
    '''
    Function to find the type of boundary that a point along a raypath falls on in the 3D block mesh defined by `shells` and `blocks`.
    Can be one of three types, or a combination of any of these, or None:
    - `T`: indicates the point falls on a latitude boundary
    - `N`: indicates the point falls on a longitude boundary
    - `Z`: indicates the point falls on a depth boundary
    - `O`: indicates the point does not fall on a boundary
    ======
    Inputs:
    - `lat`: latitude of the point of interest [format: float or int; unit: degrees]
    - `lon`: longitude of the point of interest [format: float or int; unit: degrees]
    - `z`: depth of the point of interest [format: float or int; unit: km or same as convention in `shell_dimensions.csv`]
    ======
    Outputs:
    - boundary type [format: string; unit: unitless]
    '''
    # if lon >= 180.:
    #     diff = lon - 180.
    #     lon = 180. - diff
    # elif lon < -180.:
    #     diff = -180. - lon
    #     lon = 180. - diff
    if z == CMB_depth:
        df_shell_slice = df_shells.loc[df_shells['DEPTH_MAX'] == z]
    else:
        df_shell_slice = df_shells.loc[(df_shells['DEPTH_MIN'] <= z) & (df_shells['DEPTH_MAX'] > z)]
    shell_id = int(df_shell_slice['SHELL#'])
    depth_min = float(df_shell_slice['DEPTH_MIN'])
    depth_max = float(df_shell_slice['DEPTH_MAX'])
    if lat == 90.:
        df_block_slice = df_blocks.loc[(df_blocks['LAT_MAX'] == lat) & (df_blocks['LON_MIN'] <= lon) & (df_blocks['LON_MAX'] > lon)]
    else:
        df_block_slice = df_blocks.loc[(df_blocks['LAT_MIN'] <= lat) & (df_blocks['LAT_MAX'] > lat) & (df_blocks['LON_MIN'] <= lon) & (df_blocks['LON_MAX'] > lon)] 
    block_id = int(df_block_slice['BLOCK#'])
    lat_min = float(df_block_slice['LAT_MIN'])
    lon_min = float(df_block_slice['LON_MIN'])
    if z == depth_min or z == depth_max:
        depth_boundary = 'Z'
    else:
        depth_boundary = 'O'
    if lat == lat_min and lon != lon_min:
        block_boundary = 'T'
    elif lon == lon_min and lat != lat_min:
        block_boundary = 'N'
    elif lat == lat_min and lon == lon_min:
        block_boundary = 'TN'
    else:
        block_boundary = 'O'
    if depth_boundary != 'O' and block_boundary == 'O':
        boundary = depth_boundary
    elif depth_boundary == 'O' and block_boundary != 'O':
        boundary = block_boundary
    elif depth_boundary == 'O' and block_boundary == 'O':
        boundary = 'O'
    else:
        boundary = block_boundary + depth_boundary
    return boundary

def slice_block_mins(lat_min, lon_min):
    '''
    Function to find the block ID from its minimum bounds
    ======
    Inputs:
    - `lat_min`: minimum latitude/southern bound of the block of interest [format: float or int; unit: degrees]
    - `lon_min`: minimum longitude/western bound of the block of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - block ID of block defined by `lat_min` and `lon_min` [format: int; unit: unitless]
    '''
    # if lon_min == 180.:
    #     lon_min = -180.
    df_slice = df_blocks.loc[(df_blocks['LAT_MIN'] == lat_min) & (df_blocks['LON_MIN'] == lon_min)]
    block_id = int(df_slice['BLOCK#'])
    return block_id

def slice_block_min_max(lat_min, lon_max):
    '''
    Function to find the block ID from its minimum latitude bound and maximum longitude bound
    ======
    Inputs:
    - `lat_min`: minimum latitude/southern bound of the block of interest [format: float or int; unit: degrees]
    - `lon_max`: maximum longitude/eastern bound of the block of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - block ID of block defined by `lat_min` and `lon_max` [format: int; unit: unitless]
    '''
    # if lon_max == -180.:
    #     lon_max = 180.
    df_slice = df_blocks.loc[(df_blocks['LAT_MIN'] == lat_min) & (df_blocks['LON_MAX'] == lon_max)]
    block_id = int(df_slice['BLOCK#'])
    return block_id

def slice_block_min_mid(lat_min, lon_mid):
    '''
    Function to find the block ID from its minimum latitude bound and a middle longitude
    ======
    Inputs:
    - `lat_min`: minimum latitude/southern bound of the block of interest [format: float or int; unit: degrees]
    - `lon_mid`: a longitude between the western and eastern bounds of the block of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - block ID of block defined by `lat_min` and `lon_mid` [format: int; unit: unitless]
    '''
    # if lon_mid < -180.:
    #     lon_mid == lon_mid + 360.
    df_slice = df_blocks.loc[(df_blocks['LAT_MIN'] == lat_min) & (df_blocks['LON_MIN'] <= lon_mid) & (df_blocks['LON_MAX'] > lon_mid)]
    block_id = int(df_slice['BLOCK#'])
    return block_id

def slice_block_max_mid(lat_max, lon_mid):
    '''
    Function to find the block ID from its maximum latitude bound and a middle longitude
    ======
    Inputs:
    - `lat_max`: maximum latitude/northern bound of the block of interest [format: float or int; unit: degrees]
    - `lon_mid`: a longitude between the western and eastern bounds of the block of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - block ID of block defined by `lat_max` and `lon_mid` [format: int; unit: unitless]
    '''
    # if lon_mid < -180.:
    #     lon_mid == lon_mid + 360.
    df_slice = df_block.loc[(df_blocks['LAT_MAX'] == lat_max) & (df_blocks['LON_MIN'] <= lon_mid) & (df_blocks['LON_MAX'] > lon_mid)]
    block_id = int(df_slice['BLOCK#'])
    return block_id

def slice_block_band_lon(band, lon_mid):
    '''
    Function to find the block ID from its latitude band and a middle longitude
    ======
    Inputs:
    - `band`: ID of the latitude band of the block of interest [format: int; unit: unitless]
    - `lon_mid`: a longitude between the western and eastern bounds of the block of interest [format: float or int; unit: degrees]
    ======
    Outputs:
    - block ID of block defined by `band` and `lon_mid` [format: int; unit: unitless]
    '''
    # if lon_mid < -180.:
    #     lon_mid == lon_mid + 360.
    df_slice = df_blocks.loc[(df_blocks['LAT_BAND_ID'] == band) & (df_blocks['LON_MIN'] <= lon_mid) & (df_blocks['LON_MAX'] > lon_mid)]
    block_id = int(df_slice['BLOCK#'])
    return block_id
