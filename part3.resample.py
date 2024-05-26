## Import libraries, modules ##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mod_refmodels
import mod_input
import mod_pandas
import mod_geo
import mod_boundary
import glob
import psycopg2

## Define functions ##
def add_pt_info(path_len, mid_depth, p_azimuth, p_depth, p_dist, p_seg_dist, p_lat, p_lon, p_shell, p_block):
    km_dist = path_len + cumulative_lengths[-1]
    if mod_input.physical_property == 'vs':
        physical_property = mod_refmodels.prem_vs(mid_depth)
    if mod_input.physical_property == 'vp':
        physical_property = mod_refmodels.prem_vp(mid_depth)
    seg_time = path_len / physical_property
    time = seg_time + cumulative_times[-1]
    p_bound_type = mod_pandas.find_boundary_type(p_lat, p_lon, p_depth)
    mid_pt = mod_geo.pt_from_dist(p_azimuth, pt_lats[-1], pt_lons[-1], (p_seg_dist / 2))
    seg_shell = mod_pandas.find_shell_id(mid_depth)
    seg_block = mod_pandas.find_block_id(mid_pt[0], mid_pt[1])
    pt_depths.append(p_depth)
    pt_epidists.append(p_dist)
    pt_lats.append(p_lat)
    pt_lons.append(p_lon)
    pt_shell_nos.append(p_shell)
    pt_block_nos.append(p_block)
    pt_boundary_types.append(p_bound_type)
    azimuths.append(p_azimuth)
    segment_epidists.append(p_seg_dist)
    segment_lengths.append(path_len)
    segment_times.append(seg_time)
    segment_shell_nos.append(seg_shell)
    segment_block_nos.append(seg_block)
    segment_mid_depths.append(mid_depth)
    segment_properties.append(physical_property)
    cumulative_times.append(time)
    cumulative_epidists.append(p_dist)
    cumulative_lengths.append(km_dist)
    
    
def add_bound_info(depth, dist, lat, lon, shell, block, bound_type):
    bound_depths.append(depth)
    bound_dists.append(dist)
    bound_lats.append(lat)
    bound_lons.append(lon)
    bound_shells.append(shell)
    bound_blocks.append(block)
    bound_types.append(bound_type)
    
## Start code ##
cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places
phase = mod_input.phase
file_mod = mod_input.dataset_name_mod
paths = glob.glob(f'./{mod_input.phases_directory}/{phase}{file_mod}/{mod_input.orig_paths_directory}/{phase}_*.csv')
tot_paths = len(paths)
discontinuities = mod_input.discontinuities
depth_bounds = mod_input.shell_bounds
target_length = mod_input.target_path_length
target_length_tol = mod_input.target_path_length_tolerance
lower_lim = target_length - target_length_tol
upper_lim = target_length + target_length_tol
itr = 0

# begin master loop to loop through all raypaths
for path in paths:
    itr += 1
    try:
        df_p = pd.read_csv(path)
        path_id = str(path).split(sep = '_')[-5]
        ray_param = df_p['p'][1]
        df_p = df_p.drop(columns = ['p', 'TIME'])
        max_depth = df_p['DEPTH'].max()

        # first, find any boundaries.
        # initialize lists:
        bound_dists = []
        bound_depths = []
        bound_lats = []
        bound_lons = []
        bound_shells = []
        bound_blocks = []
        bound_types = []

        # loop through each pair of points in the initial path file to find which ones have boundaries between them
        for point in range(len(df_p) - 1):
            point_2 = point + 1

            # attributes of a point on the raypath (source side):
            ps_depth = df_p['DEPTH'].iloc[point]
            ps_epidist = df_p['DIST'].iloc[point]
            ps_lat = df_p['LAT'].iloc[point]
            ps_lon = df_p['LON'].iloc[point]
            ps_shell = df_p['SHELL#'].iloc[point]
            ps_block = df_p['BLOCK#'].iloc[point]
            ps_bound_type = df_p['BOUND'].iloc[point]

            # attributes of the following point on the raypath (receiver side):
            pr_depth = df_p['DEPTH'].iloc[point_2]
            pr_epidist = df_p['DIST'].iloc[point_2]
            pr_lat = df_p['LAT'].iloc[point_2]
            pr_lon = df_p['LON'].iloc[point_2]
            pr_shell = df_p['SHELL#'].iloc[point_2]
            pr_block = df_p['BLOCK#'].iloc[point_2]
            pr_bound_type = df_p['BOUND'].iloc[point_2]

            # check for an inflection point, define a new intermediate point if necessary
            azimuth = mod_geo.azimuth(ps_lat, ps_lon, pr_lat, pr_lon)
            bazimuth = mod_geo.azimuth(pr_lat, pr_lon, ps_lat, ps_lon)
            inflection_info = mod_geo.inflection_finder(azimuth, bazimuth, ps_lat, ps_lon, pr_lat, pr_lon, ps_epidist, pr_epidist)
            
            # if there is an inflection point, retrieve its attributes and include it in the boundary finding process
            if inflection_info[0] == True:
                pi_epidist = inflection_info[1]
                pi_lat = inflection_info[2]
                pi_lon = inflection_info[3]
                pi_depth = mod_geo.new_depth(ps_depth, pr_depth, ps_epidist, pr_epidist, pi_epidist)
                pi_shell = mod_pandas.find_shell_id(pi_depth)
                pi_block = mod_pandas.find_block_id(pi_lat, pi_lon)
                pi_bound_type = mod_pandas.find_boundary_type(pi_lat, pi_lon, pi_depth)
                # add the new attributes to the original taup_path pair (in between them in the lists) to be included in boundary finding
                point_depths = [ps_depth, pi_depth, pr_depth]
                point_epidists = [ps_epidist, pi_epidist, pr_epidist]
                point_lats = [ps_lat, pi_lat, pr_lat]
                point_lons = [ps_lon, pi_lon, pr_lon]
                point_shells = [ps_shell, pi_shell, pr_shell]
                point_blocks = [ps_block, pi_block, pr_block]
                point_bound_types = [ps_bound_type, pi_bound_type, pr_bound_type]
                # add the new point to the boundary dataframe (even if it isn't a boundary, it should be retained)
                add_bound_info(pi_depth, pi_epidist, pi_lat, pi_lon, pi_shell, pi_block, pi_bound_type)

            # if there is no inflection point, move forward as usual with the original taup_path point pair
            elif inflection_info[0] == False:
                point_depths = [ps_depth, pr_depth]
                point_epidists = [ps_epidist, pr_epidist]
                point_lats = [ps_lat, pr_lat]
                point_lons = [ps_lon, pr_lon]
                point_shells = [ps_shell, pr_shell]
                point_blocks = [ps_block, pr_block]
                point_bound_types = [ps_bound_type, pr_bound_type]

            # start the main boundary finding loop
            for pair in range(len(point_depths) - 1):
                pair_2 = pair + 1
                # define the attributes of the start point
                p1_depth = point_depths[pair]
                p1_epidist = point_epidists[pair]
                p1_lat = point_lats[pair]
                p1_lon = point_lons[pair]
                p1_shell = point_shells[pair]
                p1_block = point_blocks[pair]
                p1_bound_type = point_bound_types[pair]

                # define the attributes of the end point
                p2_depth = point_depths[pair_2]
                p2_epidist = point_epidists[pair_2]
                p2_lat = point_lats[pair_2]
                p2_lon = point_lons[pair_2]
                p2_shell = point_shells[pair_2]
                p2_block = point_blocks[pair_2]
                p2_bound_type = point_bound_types[pair_2]

                # calculate the new azimuth and back azimuth
                az = mod_geo.azimuth(p1_lat, p1_lon, p2_lat, p2_lon)

                # if the boundaries are in the same block but cross at least one depth boundary:
                if p1_shell != p2_shell and p1_block == p2_block:
                    bound_info = mod_boundary.different_shell_same_block(az, p1_epidist, p1_depth, p1_lat, p1_lon, p1_shell, p1_block, p1_bound_type, p2_epidist, p2_depth, p2_lat, p2_lon, p2_shell, p2_block, p2_bound_type)
                    for bound in range(len(bound_info[0])):
                        bound_depth = bound_info[0][bound]
                        bound_dist = bound_info[1][bound]
                        bound_lat = bound_info[2][bound]
                        bound_lon = bound_info[3][bound]
                        bound_shell = bound_info[4][bound]
                        bound_block = bound_info[5][bound]
                        bound_type = bound_info[6][bound]
                        add_bound_info(bound_depth, bound_dist, bound_lat, bound_lon, bound_shell, bound_block, bound_type)

                # if the boundaries are in the same shell but cross at least one block boundary
                elif p1_shell == p2_shell and p1_block != p2_block:
                    bound_info = mod_boundary.same_shell_different_block(az, p1_epidist, p1_depth, p1_lat, p1_lon, p1_shell, p1_block, p1_bound_type, p2_epidist, p2_depth, p2_lat, p2_lon, p2_shell, p2_block, p2_bound_type)
                    for bound in range(len(bound_info[0])):
                        bound_depth = bound_info[0][bound]
                        bound_dist = bound_info[1][bound]
                        bound_lat = bound_info[2][bound]
                        bound_lon = bound_info[3][bound]
                        bound_shell = bound_info[4][bound]
                        bound_block = bound_info[5][bound]
                        bound_type = bound_info[6][bound]
                        add_bound_info(bound_depth, bound_dist, bound_lat, bound_lon, bound_shell, bound_block, bound_type)

                # if the boundaries cross at least one shell boundary and at least one block boundary
                elif p1_shell != p2_shell and p1_block != p2_block:
                    bound_info = mod_boundary.different_shell_different_block(az, p1_epidist, p1_depth, p1_lat, p1_lon, p1_shell, p1_block, p1_bound_type, p2_epidist, p2_depth, p2_lat, p2_lon, p2_shell, p2_block, p2_bound_type)
                    for bound in range(len(bound_info[0])):
                        bound_depth = bound_info[0][bound]
                        bound_dist = bound_info[1][bound]
                        bound_lat = bound_info[2][bound]
                        bound_lon = bound_info[3][bound]
                        bound_shell = bound_info[4][bound]
                        bound_block = bound_info[5][bound]
                        bound_type = bound_info[6][bound]
                        add_bound_info(bound_depth, bound_dist, bound_lat, bound_lon, bound_shell, bound_block, bound_type)
        # add everything to the same dataframe:
        df_bounds = pd.DataFrame(data = {'DIST': bound_dists, 'DEPTH': bound_depths, 'LAT': bound_lats, 'LON': bound_lons, 'SHELL#': bound_shells, 'BLOCK#': bound_blocks, 'BOUND': bound_types})
        df_all = pd.concat([df_bounds, df_p], sort = False)
        df_all.sort_values(inplace = True, by = ['DIST'])
        df_all.reset_index(inplace = True, drop = True)

        ###################################
        # END OF BOUNDARY FINDING PORTION #
        ###################################
        #   START OF RESAMPLING PORTION   #
        ###################################
        
        # resampling for path segment format outcome:
        # initialize lists to be filled with new point data
        pt_depths = []
        pt_epidists = []
        pt_lats = []
        pt_lons = []
        pt_shell_nos = []
        pt_block_nos = []
        pt_boundary_types =[]

        azimuths = []
        segment_epidists = []
        segment_lengths = []
        segment_times = []
        segment_shell_nos = []
        segment_block_nos = []
        segment_mid_depths = []
        segment_properties = [] 

        cumulative_times = []
        cumulative_epidists = []
        cumulative_lengths = []
        
        # Put initial values in point and cumulative lists
        pt_depths.append(df_p['DEPTH'].iloc[0])
        pt_epidists.append(df_p['DIST'].iloc[0])
        pt_lats.append(df_p['LAT'].iloc[0])
        pt_lons.append(df_p['LON'].iloc[0])
        pt_shell_nos.append(df_p['SHELL#'].iloc[0])
        pt_block_nos.append(df_p['BLOCK#'].iloc[0])
        pt_boundary_types.append(df_p['BOUND'].iloc[0])

        cumulative_epidists.append(df_p['DIST'].iloc[0])
        cumulative_lengths.append(df_p['DIST'].iloc[0])
        cumulative_times.append(df_p['DIST'].iloc[0])
        
        
        # initialize the while loop
        r = len(df_all)
        idx = 1
        idx_final = r - 1
        
        # begin the main loop
        while True:
            # define the current pair of points
            # point 1
            p1_depth_resamp = pt_depths[-1]
            p1_dist_resamp = pt_epidists[-1]
            p1_lat_resamp = pt_lats[-1]
            p1_lon_resamp = pt_lons[-1]
            p1_shell_resamp = pt_shell_nos[-1]
            p1_block_resamp = pt_block_nos[-1]
            p1_bound_resamp = pt_boundary_types[-1]

            # point 2
            p2_depth_resamp = df_all['DEPTH'].iloc[idx]
            p2_dist_resamp = df_all['DIST'].iloc[idx]
            p2_lat_resamp = df_all['LAT'].iloc[idx]
            p2_lon_resamp = df_all['LON'].iloc[idx]
            p2_shell_resamp = df_all['SHELL#'].iloc[idx]
            p2_block_resamp = df_all['BLOCK#'].iloc[idx]
            p2_bound_resamp = df_all['BOUND'].iloc[idx]

            # calculate the path length between the pair and the depth of the midpoint
            path_len = mod_geo.path_length(p1_depth_resamp, p2_depth_resamp, p1_dist_resamp, p2_dist_resamp)
            az_resamp = mod_geo.azimuth(p1_lat_resamp, p1_lon_resamp, p2_lat_resamp, p2_lon_resamp)
            midpoint_depth_resamp = mod_geo.midpt_depth(p1_depth_resamp, p2_depth_resamp)
            delta_dist_resamp = p2_dist_resamp - p1_dist_resamp
            
            # if the current pair is too close, move on to the next point UNLESS IT'S A SPECIAL CASE
            if path_len < lower_lim:
                # unless it's the last point in the file. then just keep it and break.
                if idx == idx_final:
                    add_pt_info(path_len, midpoint_depth_resamp, az_resamp, p2_depth_resamp, p2_dist_resamp, delta_dist_resamp, p2_lat_resamp, p2_lon_resamp, p2_shell_resamp, p2_block_resamp)
                    break
                # unless it's a model discontinuity or a bottoming point or a mesh boundary. Then just keep it and move on.
                elif (p2_depth_resamp in discontinuities) or (p2_depth_resamp == max_depth) or (p2_bound_resamp != 'O'):
                    add_pt_info(path_len, midpoint_depth_resamp, az_resamp, p2_depth_resamp, p2_dist_resamp, delta_dist_resamp, p2_lat_resamp, p2_lon_resamp, p2_shell_resamp, p2_block_resamp)
                    idx += 1
                # if it's none of these special cases and the points are still too close, move on to test the next point.
                else:
                    idx += 1
            # if the current pair is within the desired length range, just keep p2 and move on.
            elif lower_lim <= path_len <= upper_lim:
                add_pt_info(path_len, midpoint_depth_resamp, az_resamp, p2_depth_resamp, p2_dist_resamp, delta_dist_resamp, p2_lat_resamp, p2_lon_resamp, p2_shell_resamp, p2_block_resamp)
                idx += 1
                # if the second point is the last point, break out of the loop
                if idx == r:
                    break
            # if the current pair is too far apart, calculate a new point
            elif path_len > upper_lim:
                new_pt_info = mod_boundary.resample(target_length, path_len, az_resamp, delta_dist_resamp, p1_depth_resamp, p1_dist_resamp, p1_lat_resamp, p1_lon_resamp, p2_depth_resamp)
                add_pt_info(new_pt_info[0], new_pt_info[1], az_resamp, new_pt_info[2], new_pt_info[3], new_pt_info[4], new_pt_info[5], new_pt_info[6], new_pt_info[7], new_pt_info[8])
        ## While loop ends.

        ## crop for start and end point values        
        start_depths = pt_depths[:-1]
        start_epidists = pt_epidists[:-1]
        start_lats = pt_lats[:-1]
        start_lons = pt_lons[:-1]
        start_shell_nos = pt_shell_nos[:-1]
        start_block_nos = pt_block_nos[:-1]
        start_bound_types = pt_boundary_types[:-1]
        
        end_depths = pt_depths[1:]
        end_epidists = pt_epidists[1:]
        end_lats = pt_lats[1:]
        end_lons = pt_lons[1:]
        end_shell_nos = pt_shell_nos[1:]
        end_block_nos = pt_block_nos[1:]
        end_bound_types = pt_boundary_types[1:]
        
        cumulative_epidists = cumulative_epidists[1:]
        cumulative_lengths = cumulative_lengths[1:]
        cumulative_times = cumulative_times[1:]
        
        ## other calculations
        p = [ray_param] * len(start_depths)
        percent_pathlens = [((val / cumulative_lengths[-1]) * 100) for val in segment_lengths]
        percent_times = [((val / cumulative_times[-1]) * 100) for val in segment_times]
        
        df_resamp = pd.DataFrame(data = {'RAY_PARAM': p, 'START_DEPTH': start_depths, 'START_DIST': start_epidists, 'START_LAT': start_lats, 'START_LON': start_lons, 'START_SHELL#': start_shell_nos, 'START_BLOCK#': start_block_nos, 'START_TYPE': start_bound_types, 'END_DEPTH': end_depths, 'END_DIST': end_epidists, 'END_LAT': end_lats, 'END_LON': end_lons, 'END_SHELL#': end_shell_nos, 'END_BLOCK#': end_block_nos, 'END_TYPE': end_bound_types, 'AZIMUTH': azimuths, 'SEG_DIST_DEG': segment_epidists, 'SEG_DIST_KM': segment_lengths, 'SEG_TIME': segment_times, 'SEG_SHELL#': segment_shell_nos, 'SEG_BLOCK#': segment_block_nos, 'MID_DEPTH': segment_mid_depths, mod_input.segment_property_header: segment_properties, 'LENGTH_PERCENT': percent_pathlens, 'TIME_PERCENT': percent_times, 'TIME': cumulative_times, 'DIST_DEG': cumulative_epidists, 'DIST_KM': cumulative_lengths})

        h = []
        cols = df_resamp.columns
        for col in cols:
            col = str(col)
            h.append(col)
        
        genf = f'%1.{cdp}f,'
        startf = f'%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%i,%i,%s,'
        endf = f'%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%i,%i,%s,'
        segf = f'%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%i,%i,%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,%1.{cdp}f,'
        cumulativef = f'%1.{cdp}f,%1.{cdp}f,%1.{cdp}f'
        
        format_string = genf+startf+endf+segf+cumulativef
        headers = ','.join(h)

        # save the dataframe:
        np.savetxt(f'./{mod_input.phases_directory}/{phase}{file_mod}/{mod_input.resampled_paths_directory}/{phase}_{path_id}_resampled_segments.csv', df_resamp, fmt = format_string, header = headers, delimiter = ',', comments = '')

        print(f'Path number: {itr} of {tot_paths}; path: {path}')

    except:
        print(f'Path number: {itr} of {tot_paths}; path: {path}; ***** BUG *****')
        # save the name of the path to a new file
        with open(f'./{mod_input.phases_directory}/{phase}{file_mod}/{mod_input.data_directory}/{phase}{file_mod}_resample_bugs.txt', 'a') as fout:
            fout.write(f'{path}\n')
