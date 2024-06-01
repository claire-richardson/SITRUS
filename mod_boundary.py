##########
#
# Functions to calculate points along raypaths at mesh boundaries
#
# READ ME:
#
# !!!!! DO NOT MODIFY THESE FUNCTIONS !!!!!
#
# This module contains four functions written to be imported and used by the SITRUS software
# package. The three main functions compute a new point to be added to the raypath, wherever
# the raypath crosses a block boundary in the mesh defined by the first part of the SITRUS
# software. The fourth function resamples the path to specifications defined in
# `mod_input.py`. The functions in this module also depend on `mod_database.py` and
# `mod_geo.py`.
#
# the functions in this module include: 
# 1.  different_shell_same_block: Function to find the boundary between two points along a raypath that have the same Block ID but different Shell IDs
# 2.  same_shell_different_block: Function to find the boundary between two points along a raypath that have different Block IDs but the same Shell ID
# 3.  different_shell_different_block: Function to find the boundary between two points along a raypath that have different Block IDs and different Shell IDs
# 4.  resample: Function to resample a raypath to have points distributed at a specific target length
#
##########

import math
import mod_input
import mod_database
import mod_geo

depth_bounds = mod_input.shell_bounds
cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places

###############################################
#### PART OF THE CODE FOR BOUNDARY FINDING ####
###############################################

def different_shell_same_block(az, dist1, depth1, lat1, lon1, shell1, block1, type1, dist2, depth2, lat2, lon2, shell2, block2, type2):
    bd_depths = []
    bd_dists = []
    bd_lats = []
    bd_lons = []
    bd_shells = []
    bd_blocks = []
    bd_types = []
    delta_dist = dist2 - dist1
    # if the path is downgoing
    if round(depth1, rdp) < round(depth2, rdp):
        shell_boundaries = list(range(shell1 + 1, shell2 + 1))
        upgoing = False
        downgoing = True
    # if the path is upgoing
    elif round(depth1, rdp) > round(depth2, rdp):
        shell_boundaries = list(range(shell1, shell2, -1))
        upgoing = True
        downgoing = False
    # add information for each boundary
    for bound in shell_boundaries:
        if downgoing == True and 'Z' in type2 and bound == shell_boundaries[-1]:
            add_boundary_info = False
        elif upgoing == True and 'Z' in type1 and bound == shell_boundaries[0]:
            add_boundary_info = False
        else:
            add_boundary_info = True
            ref_shell_info = mod_database.get_shell_info(bound)
            bound_depth = ref_shell_info[1]
            if downgoing == True:
                bound_ratio = (bound_depth - depth1) / (depth2 - depth1)
            elif upgoing == True:
                bound_ratio = (depth1 - bound_depth) / (depth1 - depth2)
            bound_dist = dist1 + (delta_dist * bound_ratio)
            bound_pt = mod_geo.pt_from_dist(az, lat1, lon1, (delta_dist * bound_ratio))
            bound_lat = bound_pt[0]
            bound_lon = bound_pt[1]
            bound_type = mod_database.find_boundary_type(bound_lat, bound_lon, bound_depth)
        if add_boundary_info == True:
            bd_depths.append(round(bound_depth, cdp))
            bd_dists.append(round(bound_dist, cdp))
            bd_lats.append(round(bound_lat, cdp))
            bd_lons.append(round(bound_lon, cdp))
            bd_shells.append(bound)
            bd_blocks.append(block1)
            bd_types.append(bound_type)

    for n in range(len(bd_lons)):
        if bd_lons[n] >= 180.:
            bd_lons[n] -= 360.
        if bd_lons[n] < -180:
            bd_lons[n] += 360
    
    bd_info = [bd_depths, bd_dists, bd_lats, bd_lons, bd_shells, bd_blocks, bd_types]
    return bd_info




def same_shell_different_block(az, dist1, depth1, lat1, lon1, shell1, block1, type1, dist2, depth2, lat2, lon2, shell2, block2, type2):
    bd_depths = []
    bd_dists = []
    bd_lats = []
    bd_lons = []
    bd_shells = []
    bd_blocks = []
    bd_types = []

    block1_info = mod_database.get_block_info(block1)
    block2_info = mod_database.get_block_info(block2)
    band1 = block1_info[1]
    band2 = block2_info[1]

    # if the path stays in the same latitude band, i.e., the path crosses only a longitude boundary or starts on a pole:
    if band1 == band2:
        # initialize while loop
        ref_lat = lat1
        ref_lon = lon1
        ref_block = block1
        while True:
            # if the path crosses a pole, no need to add a new point.
            if (round(az, rdp) == 180. and round(lat1, rdp) == 90.) or (round(az, rdp) == 0. and round(lat1, rdp) == -90.):
                break
            # if it doesn't cross a pole, find the longitude boundary
            else:
                ref_block_info = mod_database.get_block_info(ref_block)
                # if the path is heading east
                if 0. < round(az, rdp) < 180.:
                    bound_lon = ref_block_info[5]
                    if bound_lon == 180.:
                        bound_lon = -180.    
                    bound_pt = mod_geo.known_lon(az, lat1, lon1, bound_lon)
                    # if lon2 is on the found longitude boundary, break (special case).
                    if round(lon2, rdp) == round(bound_pt[1], rdp):
                        break
                    bound_block = mod_database.slice_block_mins(ref_block_info[2], ref_block_info[5])
                    ref_block = bound_block
                    bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                    bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                    bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                    bd_depths.append(round(bound_depth, cdp))
                    bd_dists.append(round(bound_dist, cdp))
                    bd_lats.append(round(bound_pt[0], cdp))
                    bd_lons.append(round(bound_pt[1], cdp))
                    bd_shells.append(shell1)
                    bd_blocks.append(bound_block)
                    bd_types.append(bound_type)
                    # if the new reference block is the same as block2, then break.
                    if ref_block == block2:
                        break
                # if the path is heading west
                elif 180. < round(az, rdp) < 360.:
                    if round(ref_lon, rdp) == round(block2_info[5], rdp) or (round(ref_lon, rdp) == -180. and round(block2_info[5], rdp) == 180.):
                        break
                    bound_lon = ref_block_info[4]
                    bound_pt = mod_geo.known_lon(az, lat1, lon1, bound_lon)
                    bound_block = ref_block
                    bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                    bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                    bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                    bd_depths.append(round(bound_depth, cdp))
                    bd_dists.append(round(bound_dist, cdp))
                    bd_lats.append(round(bound_pt[0], cdp))
                    bd_lons.append(round(bound_pt[1], cdp))
                    bd_shells.append(shell1)
                    bd_blocks.append(bound_block)
                    bd_types.append(bound_type)
                    if ref_block_info[4] == -180.:
                        ref_block = mod_database.slice_block_min_max(ref_block_info[2], 180.)
                    else:
                        ref_block = mod_database.slice_block_min_max(ref_block_info[2], ref_block_info[4])
                # next iteration for while loop
                ref_lat = bound_pt[0]
                ref_lon = bound_pt[1]
    # if the path crosses at least one latitude boundary
    else:
        # if the path is going north and crosses at least one lat boundary:
        if band1 < band2:
            lat_bands = list(range(band1, band2 + 1))
            lat_boundaries = lat_bands[1:]
        # if the path is going south and crosses at least one lat boundary:
        elif band1 > band2:
            lat_bands = list(range(band1, band2 - 1, -1))
            lat_boundaries = lat_bands[:-1]
        # deal with longitude changes for points on the poles
        if (round(az, rdp) == 180. and round(lat1, rdp) == 90.) or (round(az, rdp) == 0. and round(lat1, rdp) == -90.):
            if lon1 < 0.:
                lon1 += 180.                
            elif lon1 >= 0.:
                lon1 -= 180.
        # initialize loops
        lat_ref_lat = lat1
        lat_ref_lon = lon1
        # begin loops
        for band in lat_bands:
            lat_ref_block = mod_database.slice_block_band_lon(band, lat_ref_lon)
            lat_ref_block_info = mod_database.get_block_info(lat_ref_block)
            lat_ref_lon_min = lat_ref_block_info[4]
            lat_ref_lon_max = lat_ref_block_info[5]
            # if the path is going directly south
            if round(az, rdp) == 180.:
                lat_bound_pt = [lat_ref_block_info[2], lat_ref_lon]
                lon_check = False
                if lat_ref_block == block2 or (round(lat_bound_pt[0], rdp) == round(lat1, rdp) and round(lat_bound_pt[1], rdp) == round(lon1, rdp)):
                    add_lat_boundary = False
                else:
                    lat_bound_block = mod_database.find_block_id(lat_bound_pt[0], lat_bound_pt[1])
                    lat_bound_dist = dist1 + (lat1 - lat_bound_pt[0])                
                    lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                    lat_bound_shell = shell1
                    add_lat_boundary = True
            # if the path is going directly north
            elif round(az, rdp) == 0.:
                lat_bound_pt = [lat_ref_block_info[3], lat_ref_lon]
                lon_check = False
                if round(lat_bound_pt[0], rdp) >= round(lat2, rdp):
                    add_lat_boundary = False
                else:
                    lat_bound_block = mod_database.find_block_id(lat_bound_pt[0], lat_bound_pt[1])
                    lat_bound_dist = dist1 + (abs(lat1 - lat_bound_pt[0]))
                    lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                    lat_bound_shell = shell1
                    add_lat_boundary = True
            # if the path is going south and the code is in the first band and the first point is on the bottom of the lat band    
            elif band1 > band2 and band == lat_bands[0] and round(lat1, rdp) == round(lat_ref_block_info[2], rdp):
                lon_check = False
                add_lat_boundary = False
                lat_bound_pt = [lat1, lon1]
            # if the latitude reference point is also on a longitude line (accounts for corners)
            elif round(lat_ref_lon, rdp) == round(lat_ref_lon_min, rdp) and 180. < round(az, rdp) < 360.:
                if lat_ref_lon_min == -180.:
                    test_block = mod_database.slice_block_min_max(lat_ref_block_info[2], 180.)
                else:
                    test_block = mod_database.slice_block_min_max(lat_ref_block_info[2], lat_ref_lon_min)
                test_block_info = mod_database.get_block_info(test_block)
                if block2 == test_block:
                    lon_check = False
                    add_lat_boundary = False
                    lat_bound_pt = [lat2, lon2]
                elif block2 != test_block and band == band2:
                    lon_check = True
                    add_lat_boundary = False
                    lat_bound_pt = [lat2, lon2]
                elif band != band2:
                    # if the path is going to the northwest
                    if round(az, rdp) > 270.:
                        lat_bound_lat = lat_ref_block_info[3]
                        if round(lat_bound_lat, rdp) == round(lat2, rdp):
                            add_lat_boundary = False
                            lat_bound_pt = [lat2, lon2]
                        else:
                            add_lat_boundary = True
                            lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                    # if the path is going to the southwest
                    elif round(az, rdp) < 270.:
                        lat_bound_lat = lat_ref_block_info[2]
                        if round(lat_bound_lat, rdp) <= round(lat2, rdp):
                            add_lat_boundary = False
                            lat_bound_pt = [lat2, lon2]
                        else:
                            add_lat_boundary = True
                            lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                    lat_bound_block = mod_database.slice_block_min_mid(lat_bound_pt[0], lat_bound_pt[1])
                    lat_bound_shell = shell1
                    lat_bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, lat_bound_pt[0], lat_bound_pt[1])
                    lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                    if round(test_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) < round(test_block_info[5], rdp):
                        lon_check = False
                    else:
                        lon_check = True
            else:
                if band != band2:
                    if band < band2:
                        lat_bound_lat = lat_ref_block_info[3]
                        if lat_bound_lat == lat2:
                            add_lat_boundary = False
                            lat_bound_pt = [lat2, lon2]
                        else:
                            add_lat_boundary = True
                            lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                    elif band > band2:
                        lat_bound_lat = lat_ref_block_info[2]
                        if round(lat_bound_lat, rdp) <= round(lat2, rdp):
                            add_lat_boundary = False
                            lat_bound_pt = [lat2, lon2]
                        else:
                            add_lat_boundary = True
                            lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                    if add_lat_boundary == True:
                        lat_bound_block = mod_database.slice_block_min_mid(lat_bound_pt[0], lat_bound_pt[1])
                        lat_bound_shell = shell1
                        lat_bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, lat_bound_pt[0], lat_bound_pt[1])
                        lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                elif band == band2:
                    lat_bound_pt = [lat2, lon2]
                    add_lat_boundary = False
                if round(lat_ref_lon_min, rdp) <= round(lat_bound_pt[1], rdp) < round(lat_ref_lon_max, rdp):
                    lon_check = False
                else:
                    lon_check = True
            if lon_check == True:
                # if the path crosses a longitude boundary before reaching the top of the latitude band
                # initialize the while loop
                ref_block = lat_ref_block
                ref_lat = lat_ref_lat
                ref_lon = lat_ref_lon
                # longitude boundary to the east (path is eastbound)
                if 0. < round(az, rdp) < 180.:
                    while True:
                        ref_block_info = mod_database.get_block_info(ref_block)
                        if (round(ref_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) <= round(ref_block_info[5], rdp)) or (ref_block_info[5] == 180. and round(lat_bound_pt[1], rdp) == -180.):
                            break
                        else:
                            bound_pt = mod_geo.known_lon(az, lat1, lon1, ref_block_info[5])
                            bound_block = mod_database.find_block_id(bound_pt[0], bound_pt[1])
                            ref_block = bound_block
                            bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                            bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                            bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                            bd_depths.append(round(bound_depth, cdp))
                            bd_dists.append(round(bound_dist, cdp))
                            bd_lats.append(round(bound_pt[0], cdp))
                            bd_lons.append(round(bound_pt[1], cdp))
                            bd_shells.append(shell1)
                            bd_blocks.append(bound_block)
                            bd_types.append(bound_type)
                            # next iteration for while loop
                            ref_lat = bound_pt[0]
                            ref_lon = bound_pt[1]
                # longitude boundary to the west (path is westbound)
                elif 180. < round(az, rdp) < 360.:
                    while True:
                        ref_block_info = mod_database.get_block_info(ref_block)
                        if round(ref_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) < round(ref_block_info[5], rdp):
                            break
                        else:
                            bound_pt = mod_geo.known_lon(az, lat1, lon1, ref_block_info[4])
                            bound_block = ref_block
                            ref_block = mod_database.slice_block_min_max(ref_block_info[2], ref_block_info[4])
                            bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                            bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                            bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                            bd_depths.append(round(bound_depth, cdp))
                            bd_dists.append(round(bound_dist, cdp))
                            bd_lats.append(round(bound_pt[0], cdp))
                            bd_lons.append(round(bound_pt[1], cdp))
                            bd_shells.append(shell1)
                            bd_blocks.append(bound_block)
                            bd_types.append(bound_type)
                            ref_lat = bound_pt[0]
                            ref_lon = bound_pt[1]
            # add the latitude boundary
            if add_lat_boundary == True:
                lat_bound_type = mod_database.find_boundary_type(lat_bound_pt[0], lat_bound_pt[1], lat_bound_depth)
                bd_depths.append(round(lat_bound_depth, cdp))
                bd_dists.append(round(lat_bound_dist, cdp))
                bd_lats.append(round(lat_bound_pt[0], cdp))
                bd_lons.append(round(lat_bound_pt[1], cdp))
                bd_shells.append(lat_bound_shell)
                bd_blocks.append(lat_bound_block)
                bd_types.append(lat_bound_type)
            lat_ref_lat = lat_bound_pt[0]
            lat_ref_lon = lat_bound_pt[1]
    
    for n in range(len(bd_lons)):
        if bd_lons[n] >= 180.:
            bd_lons[n] -= 360.
        if bd_lons[n] < -180:
            bd_lons[n] += 360

    bd_info = [bd_depths, bd_dists, bd_lats, bd_lons, bd_shells, bd_blocks, bd_types]
    return bd_info


def different_shell_different_block(az, dist1, depth1, lat1, lon1, shell1, block1, type1, dist2, depth2, lat2, lon2, shell2, block2, type2):
    bd_depths = []
    bd_dists = []
    bd_lats = []
    bd_lons = []
    bd_shells = []
    bd_blocks = []
    bd_types = []
    delta_dist = dist2 - dist1
    # if the path is downgoing
    if round(depth1, rdp) < round(depth2, rdp):
        shell_boundaries = list(range(shell1 + 1, shell2 + 1))
        depth_bands = list(range(shell1, shell2 + 1))
        upgoing = False
        downgoing = True
    # if the path is upgoing
    elif round(depth1, rdp) > round(depth2, rdp):
        shell_boundaries = list(range(shell1, shell2, -1))
        depth_bands = list(range(shell1, shell2 - 1, -1))
        upgoing = True
        downgoing = False

    # initialize the shell finding for-loop
    init_block = block1
    init_lat = lat1
    init_lon = lon1
    for depth_band in depth_bands:
        init_shell_info = mod_database.get_shell_info(depth_band)
        init_block_info = mod_database.get_block_info(init_block)
        # if the path is downgoing and the last point is in depth bounds, AND this is the last loop
        # check if there is a block boundary between
        if downgoing == True and 'Z' in type2 and depth_band == depth_bands[-1]:
            add_shell_boundary = False
            block_check = False
        elif downgoing == True and 'Z' in type2 and depth_band == depth_bands[-2]:
            shell_bound_info = [depth2, dist2, lat2, lon2, shell2, block2]
            add_shell_boundary = False
            if init_block != block2:
                block_check = True
                depth_boundary_lat = lat2
                depth_boundary_lon = lon2
                depth_boundary_block = block2
            else:
                block_check = False
        # if the path is upgoing and the first point is in depth bounds and this is the first loop,
        elif upgoing == True and 'Z' in type1 and depth_band == depth_bands[0]:
            shell_bound_info = [depth1, dist1, lat1, lon1, shell1, block1]
            block_check = False
            add_shell_boundary = False
        elif depth2 not in depth_bounds and depth_band == depth_bands[-1]:
            shell_bound_info = [depth2, dist2, lat2, lon2, shell2, block2]
            depth_boundary_depth = depth2
            depth_boundary_dist = dist2
            depth_boundary_block = block2
            depth_boundary_lat = lat2
            depth_boundary_lon = lon2
            add_shell_boundary = False
            if init_block != block2:
                block_check = True
            else:
                block_check = False
        else:
            if downgoing == True:                
                depth_boundary_depth = init_shell_info[3]
                depth_boundary_ratio = (depth_boundary_depth - depth1) / (depth2 - depth1)
                depth_boundary_shell = mod_database.slice_shell_top(depth_boundary_depth)
            elif upgoing == True:
                depth_boundary_depth = init_shell_info[1]
                depth_boundary_ratio = (depth1 - depth_boundary_depth) / (depth1 - depth2)
                depth_boundary_shell = depth_band
            depth_boundary_dist = dist1 + (delta_dist * depth_boundary_ratio)
            depth_boundary_pt = mod_geo.pt_from_dist(az, lat1, lon1, (delta_dist * depth_boundary_ratio))
            depth_boundary_lat = depth_boundary_pt[0]
            depth_boundary_lon = depth_boundary_pt[1]

            if (round(init_block_info[2], rdp) <= round(depth_boundary_lat, rdp) < round(init_block_info[3], rdp)) and (round(init_block_info[4], rdp) <= round(depth_boundary_lon, rdp) < round(init_block_info[5], rdp) or (round(depth_boundary_lon, rdp) == -180. and init_block_info[5] == 180.)):
                depth_boundary_block = init_block
                block_check = False
            else:
                depth_boundary_block = mod_database.find_block_id(depth_boundary_lat, depth_boundary_lon)
                block_check = True
            shell_bound_info = [depth_boundary_depth, depth_boundary_dist, depth_boundary_lat, depth_boundary_lon, depth_boundary_shell, depth_boundary_block]

            if round(depth_boundary_depth, rdp) == round(depth2, rdp):
                add_shell_boundary = False
            else:
                add_shell_boundary = True

        if block_check == True:
            block1_info = mod_database.get_block_info(init_block)
            block2_info = mod_database.get_block_info(depth_boundary_block)
            band1 = block1_info[1]
            band2 = block2_info[1]

            # if the path stays in the same latitude band
            if band1 == band2:
                # initialize while loop
                ref_lat = init_lat
                ref_lon = init_lon
                ref_block = init_block
                while True:
                    # if the path crosses a pole, no need to add a new point
                    if (round(az, rdp) == 180. and round(lat1, rdp) == 90.) or (round(az, rdp) == 0. and round(lat1, rdp) == -90.):
                        break
                    # if it doesn't cross a pole, find the longitude boundary
                    else:
                        ref_block_info = mod_database.get_block_info(ref_block)
                        # if the path is heading east
                        if 0. < round(az, rdp) < 180.:
                            bound_lon = ref_block_info[5]
                            if bound_lon == 180.:
                                bound_lon = -180.    
                            bound_pt = mod_geo.known_lon(az, lat1, lon1, bound_lon)
                            # if the shell boundary lon is on the found longitude boundary, break (special case).
                            if round(depth_boundary_lon, rdp) == round(bound_pt[1], rdp):
                                break   
                            bound_block = mod_database.slice_block_mins(ref_block_info[2], ref_block_info[5])
                            ref_block = bound_block
                            bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                            bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                            bound_shell = mod_database.find_shell_id(bound_depth)
                            bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                            bd_depths.append(round(bound_depth, cdp))
                            bd_dists.append(round(bound_dist, cdp))
                            bd_lats.append(round(bound_pt[0], cdp))
                            bd_lons.append(round(bound_pt[1], cdp))
                            bd_shells.append(bound_shell)
                            bd_blocks.append(bound_block)
                            bd_types.append(bound_type)
                            # if the new reference block is the same as block2, then break.
                            if ref_block == depth_boundary_block:
                                break
                        # if the path is heading west
                        elif 180. < round(az, rdp) < 360.:
                            if round(ref_lon, rdp) == round(block2_info[5], rdp) or (round(ref_lon, rdp) == -180. and block2_info[5] == 180.):
                                break
                            bound_lon = ref_block_info[4]
                            bound_pt = mod_geo.known_lon(az, lat1, lon1, bound_lon)
                            bound_block = ref_block
                            bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                            bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                            bound_shell = mod_database.find_shell_id(bound_depth)
                            bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                            bd_depths.append(round(bound_depth, cdp))
                            bd_dists.append(round(bound_dist, cdp))
                            bd_lats.append(round(bound_pt[0], cdp))
                            bd_lons.append(round(bound_pt[1], cdp))
                            bd_shells.append(bound_shell)
                            bd_blocks.append(bound_block)
                            bd_types.append(bound_type)
                            if ref_block_info[4] == -180.:
                                ref_block = mod_database.slice_block_min_max(ref_block_info[2], 180.)
                            else:
                                ref_block = mod_database.slice_block_min_max(ref_block_info[2], ref_block_info[4])
                        # next iteration for while loop
                        ref_lat = bound_pt[0]
                        ref_lon = bound_pt[1]
            # if the path crosses at least one latitude boundary
            else:
                # if the path is going north and crosses at least one lat boundary:
                if band1 < band2:
                    lat_bands = list(range(band1, band2 + 1))
                    lat_boundaries = lat_bands[1:]
                # if the path is going south and crosses at least one lat boundary:
                elif band1 > band2:
                    lat_bands = list(range(band1, band2 - 1, -1))
                    lat_boundaries = lat_bands[:-1]
                # deal with longitude changes for points on the poles
                if (round(az, rdp) == 180. and round(init_lat, rdp) == 90.) or (round(az, rdp) == 0. and round(init_lat, rdp) == -90.):
                    if init_lon < 0.:
                        init_lon += 180.                
                    elif init_lon >= 0.:
                        init_lon -= 180.
                # initialize loops
                lat_ref_lat = init_lat
                lat_ref_lon = init_lon
                # begin loops
                for band in lat_bands:
                    lat_ref_block = mod_database.slice_block_band_lon(band, lat_ref_lon)
                    lat_ref_block_info = mod_database.get_block_info(lat_ref_block)
                    lat_ref_lon_min = lat_ref_block_info[4]
                    lat_ref_lon_max = lat_ref_block_info[5]
                    # if the path is going directly south
                    if round(az, rdp) == 180.:
                        lat_bound_pt = [lat_ref_block_info[2], lat_ref_lon]
                        lon_check = False
                        if lat_ref_block == block2 or (round(lat_bound_pt[0], rdp) == round(lat1, rdp) and round(lat_bound_pt[1], rdp) == round(lon1, rdp)):
                            add_lat_boundary = False
                        else:
                            lat_bound_block = mod_database.find_block_id(lat_bound_pt[0], lat_bound_pt[1])
                            lat_bound_dist = dist1 + (lat1 - lat_bound_pt[0])                
                            lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                            lat_bound_shell = mod_database.find_shell_id(lat_bound_depth)
                            add_lat_boundary = True
                    # if the path is going directly north
                    elif round(az, rdp) == 0.:
                        lat_bound_pt = [lat_ref_block_info[3], lat_ref_lon]
                        lon_check = False
                        if round(lat_bound_pt[0], rdp) >= round(depth_boundary_lat, rdp):
                            add_lat_boundary = False
                        else:
                            lat_bound_block = mod_database.find_block_id(lat_bound_pt[0], lat_bound_pt[1])
                            lat_bound_dist = dist1 + (abs(lat1 - lat_bound_pt[0]))
                            lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                            lat_bound_shell = mod_database.find_shell_id(lat_bound_depth)
                            add_lat_boundary = True
                    # if the path is going south and the code is in the first band and the first point is on the bottom of the lat band    
                    elif band1 > band2 and band == lat_bands[0] and round(lat_ref_lat, rdp) == round(lat_ref_block_info[2], rdp):    
                        lon_check = False
                        add_lat_boundary = False
                        lat_bound_pt = [init_lat, init_lon]
                    # if the latitude reference point is also on a longitude line (accounts for corners)
                    elif round(lat_ref_lon, rdp) == round(lat_ref_lon_min, rdp) and 180. < round(az, rdp) < 360.:
                        if lat_ref_lon_min == -180.:
                            test_block = mod_database.slice_block_min_max(lat_ref_block_info[2], 180.)
                        else:
                            test_block = mod_database.slice_block_min_max(lat_ref_block_info[2], lat_ref_lon_min)
                        test_block_info = mod_database.get_block_info(test_block)
                        if depth_boundary_block == test_block:
                            lon_check = False
                            add_lat_boundary = False
                            lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                        elif depth_boundary_block != test_block and band == band2:
                            lon_check = True
                            add_lat_boundary = False
                            lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                        elif band != band2:
                            # if the path is going to the northwest
                            if round(az, rdp) > 270.:
                                lat_bound_lat = lat_ref_block_info[3]
                                if round(lat_bound_lat, rdp) == round(depth_boundary_lat, rdp):
                                    add_lat_boundary = False
                                    lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                                else:
                                    add_lat_boundary = True
                                    lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                            # if the path is going to the southwest
                            elif round(az, rdp) < 270.:
                                lat_bound_lat = lat_ref_block_info[2]
                                if round(lat_bound_lat, rdp) <= round(depth_boundary_lat, rdp):
                                    add_lat_boundary = False
                                    lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                                else:
                                    add_lat_boundary = True
                                    lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                            lat_bound_block = mod_database.slice_block_min_mid(lat_bound_pt[0], lat_bound_pt[1])
                            lat_bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, lat_bound_pt[0], lat_bound_pt[1])
                            lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                            lat_bound_shell = mod_database.find_shell_id(lat_bound_depth)
                            if round(test_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) < round(test_block_info[5], rdp):
                                lon_check = False
                            else:
                                lon_check = True
                    else:
                        if band != band2:
                            # if the path is going north
                            if band < band2:
                                lat_bound_lat = lat_ref_block_info[3]
                                if round(lat_bound_lat, rdp) == round(depth_boundary_lat, rdp):
                                    add_lat_boundary = False
                                    lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                                else:
                                    add_lat_boundary = True
                                    lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)
                            # if the path is giong south
                            elif band > band2:
                                lat_bound_lat = lat_ref_block_info[2]
                                if round(lat_bound_lat, rdp) <= round(depth_boundary_lat, rdp):
                                    add_lat_boundary = False
                                    lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                                else:
                                    add_lat_boundary = True
                                    lat_bound_pt = mod_geo.known_lat(az, lat1, lon1, lat_bound_lat)

                            if add_lat_boundary == True:
                                lat_bound_block = mod_database.slice_block_min_mid(lat_bound_pt[0], lat_bound_pt[1])
                                lat_bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, lat_bound_pt[0], lat_bound_pt[1])
                                lat_bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, lat_bound_dist)
                                lat_bound_shell = mod_database.find_shell_id(lat_bound_depth)
                        elif band == band2:
                            lat_bound_pt = [depth_boundary_lat, depth_boundary_lon]
                            add_lat_boundary = False
                        if round(lat_ref_lon_min, rdp) <= round(lat_bound_pt[1], rdp) < round(lat_ref_lon_max, rdp):
                            lon_check = False
                        else:
                            lon_check = True
                    # if the path crosses a longitude boundary before reaching the top of the latitude band
                    # longitude boundary to the east (path is eastbound)
                    if lon_check == True:
                        # initialize the while loop
                        ref_block = lat_ref_block
                        ref_lat = lat_ref_lat
                        ref_lon = lat_ref_lon
                        if 0. < round(az, rdp) < 180.:
                            while True:
                                ref_block_info = mod_database.get_block_info(ref_block)
                                if (round(ref_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) <= round(ref_block_info[5], rdp)) or (ref_block_info[5] == 180. and round(lat_bound_pt[1], rdp) == -180.):
                                    break
                                else:
                                    bound_pt = mod_geo.known_lon(az, lat1, lon1, ref_block_info[5])
                                    bound_block = mod_database.find_block_id(bound_pt[0], bound_pt[1])
                                    ref_block = bound_block
                                    bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                                    bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                                    bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                                    bound_shell = mod_database.find_shell_id(bound_depth)
                                    bd_depths.append(round(bound_depth, cdp))
                                    bd_dists.append(round(bound_dist, cdp))
                                    bd_lats.append(round(bound_pt[0], cdp))
                                    bd_lons.append(round(bound_pt[1], cdp))
                                    bd_shells.append(bound_shell)
                                    bd_blocks.append(bound_block)
                                    bd_types.append(bound_type)
                                    # next iteration for while loop
                                    ref_lat = bound_pt[0]
                                    ref_lon = bound_pt[1]
                        # longitude boundary to the west (path is westbound)
                        elif 180. < round(az, rdp) < 360.:
                            while True:
                                ref_block_info = mod_database.get_block_info(ref_block)
                                if round(ref_block_info[4], rdp) <= round(lat_bound_pt[1], rdp) < round(ref_block_info[5], rdp):
                                    break
                                else:
                                    bound_pt = mod_geo.known_lon(az, lat1, lon1, ref_block_info[4])
                                    bound_block = ref_block
                                    ref_block = mod_database.slice_block_min_max(ref_block_info[2], ref_block_info[4])
                                    bound_dist = dist1 + mod_geo.GCP_length(lat1, lon1, bound_pt[0], bound_pt[1])
                                    bound_depth = mod_geo.new_depth(depth1, depth2, dist1, dist2, bound_dist)
                                    bound_type = mod_database.find_boundary_type(bound_pt[0], bound_pt[1], bound_depth)
                                    bound_shell = mod_database.find_shell_id(bound_depth)
                                    bd_depths.append(round(bound_depth, cdp))
                                    bd_dists.append(round(bound_dist, cdp))
                                    bd_lats.append(round(bound_pt[0], cdp))
                                    bd_lons.append(round(bound_pt[1], cdp))
                                    bd_shells.append(bound_shell)
                                    bd_blocks.append(bound_block)
                                    bd_types.append(bound_type)
                                    ref_lat = bound_pt[0]
                                    ref_lon = bound_pt[1]    
                    # add the latitude boundary
                    if add_lat_boundary == True:
                        lat_bound_type = mod_database.find_boundary_type(lat_bound_pt[0], lat_bound_pt[1], lat_bound_depth)
                        bd_depths.append(round(lat_bound_depth, cdp))
                        bd_dists.append(round(lat_bound_dist, cdp))
                        bd_lats.append(round(lat_bound_pt[0], cdp))
                        bd_lons.append(round(lat_bound_pt[1], cdp))
                        bd_shells.append(lat_bound_shell)
                        bd_blocks.append(lat_bound_block)
                        bd_types.append(lat_bound_type)
                    lat_ref_lat = lat_bound_pt[0]
                    lat_ref_lon = lat_bound_pt[1]
        if add_shell_boundary == True:
            bound_type = mod_database.find_boundary_type(shell_bound_info[2], shell_bound_info[3], shell_bound_info[0])
            bd_depths.append(round(shell_bound_info[0], cdp))
            bd_dists.append(round(shell_bound_info[1], cdp))
            bd_lats.append(round(shell_bound_info[2], cdp))
            bd_lons.append(round(shell_bound_info[3], cdp))
            bd_shells.append(shell_bound_info[4])
            bd_blocks.append(shell_bound_info[5])
            bd_types.append(bound_type)
        init_block = shell_bound_info[5]
        init_lat = shell_bound_info[2]
        init_lon = shell_bound_info[3]

    for n in range(len(bd_lons)):
        if bd_lons[n] >= 180.:
            bd_lons[n] -= 360.
        if bd_lons[n] < -180:
            bd_lons[n] += 360
    
    bd_info = [bd_depths, bd_dists, bd_lats, bd_lons, bd_shells, bd_blocks, bd_types]
    return bd_info

    
def resample(target_len, path_len, az, delta_dist, p1_depth, p1_dist, p1_lat, p1_lon, p2_depth):    
    ratio = target_len / path_len
    delta_depth = abs(p2_depth - p1_depth)
    seg_depth = delta_depth * ratio
    seg_dist = delta_dist * ratio
    seg_extent = path_len * ratio
    new_coords = mod_geo.pt_from_dist(az, p1_lat, p1_lon, seg_dist)
    new_lat = new_coords[0]
    new_lon = new_coords[1]                
    if p1_depth < p2_depth:
        new_depth = p1_depth + seg_depth
    elif p1_depth > p2_depth:
        new_depth = p1_depth - seg_depth
    elif p1_depth == p2_depth:
        new_depth = p1_depth
    new_shell = mod_database.find_shell_id(new_depth)
    new_block = mod_database.find_block_id(new_lat, new_lon)
    midpoint_depth = mod_geo.midpt_depth(p1_depth, new_depth)
    new_dist = p1_dist + seg_dist
    return [seg_extent, midpoint_depth, new_depth, new_dist, seg_dist, new_lat, new_lon, new_shell, new_block]