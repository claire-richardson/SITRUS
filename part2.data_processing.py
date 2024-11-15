import os
import time
import glob
import shuitl
import mod_geo
import fnmatch
import mod_input
import mod_track
import mod_database
import mod_boundary
import numpy as np
import pandas as pd
import mod_refmodels
import multiprocessing as mp
from obspy.taup import taup_geo
from obspy.taup import TauPyModel
from ellipticipy import ellipticity_correction

pid = os.getpid()
start_time = time.time()
start_subj = f'Process began (PID: {pid}); part2.data_processing.py'
start_text = f'Process {pid} began;\nDataset: {mod_input.dataset}'
try:
    mod_track.SendMsg(start_subj, start_text)
except:
    pass

phases_directory = mod_input.phases_directory
data_directory = mod_input.data_directory
paths_directory = mod_input.paths_directory
model = TauPyModel(model = mod_input.reference_model)
rdp = mod_input.rounded_decimal_places

try:
    os.mkdir(phases_directory)
except:
    pass

if mod_input.data_wave_type == 'S':
    df_crust = pd.read_csv(f'./{mod_input.tomography_model_directory}/crust/CRUST_1.0_vsh.csv')
elif mod_input.data_wavey_type == 'P':
    df_crust = pd.read_csv(f'./{mod_input.tomography_model_directory}/crust/CRUST_1.0_vp.csv')
df_data = pd.read_csv(mod_input.dataset)
all_phases = mod_input.all_phases

headers = mod_input.main_headers
for header in mod_input.raw_headers_to_keep:
    headers.append(header)

# make directories & subdirectories for each phase in master file; separate phase information
for p in all_phases:
    phase_name = p.split('_')[0]
    try:
        os.mkdir(f'./{phases_directory}/{p}')
        os.mkdir(f'./{phases_directory}/{p}/{data_directory}')
        os.mkdir(f'./{phases_directory}/{p}/{paths_directory}')
        # os.mkdir(f'./{phases_directory}/{p}/{orig_paths_directory}')
        # os.mkdir(f'./{phases_directory}/{p}/{resampled_paths_directory}')
        df_specific_phase = df_data[headers].copy()
        df_specific_phase = df_specific_phase.loc[df_specific_phase['PHASE'] == phase_name]
        df_specific_phase['DIST'] = 0.
        df_specific_phase['ELLIP_CORR'] = 0.
        df_specific_phase['ELLIP_DT'] = 0.
        df_specific_phase['CRUST_1.0_CORR'] = 0.
        df_specific_phase['CRUST_1.0_DT'] = 0.
        df_specific_phase['CRUST_1.0_ELLIP_DT'] = 0.
        df_specific_phase['PATH_ID'] = 0
        df_specific_phase.to_csv(f'./{phases_directory}/{p}/{data_directory}/{phase_name}_master_data.csv', index = False)
    except:
        pass


# define function for making raypaths and adding ellipticity corrections:
def make_raypaths(phase):
    phase_name = phase.split('_')[0]
    try:
        os.remove(f'./{data_directory}/{phase}/{data_directory}/{phase_name}_pt2_bugs_make_raypaths.txt')
    except:
        pass
    start_raypaths = time.time()
    
    with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_bugs_make_raypaths.txt', 'w') as fout:
        fout.write(f'EQ_DEPTH,EQ_LAT,EQ_LON,STA_LAT,STA_LON,BUG\n')
    
    df_data = pd.read_csv(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_master_data.csv')
    
    path_id = 0
    for path in range(len(df_data)):
        with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_log_make_raypaths.txt', 'a') as fout:
            fout.write(f'- working on path {path + 1} of {len(df_data)}\n')

        slat = df_data['STA_LAT'].iloc[path]
        slon = df_data['STA_LON'].iloc[path]
        elat = df_data['EQ_LAT'].iloc[path]
        elon = df_data['EQ_LON'].iloc[path]
        depth = df_data['EQ_DEP'].iloc[path]
        dt = df_data['DT'].iloc[path]
        
        path_az = mod_geo.azimuth(elat, elon, slat, slon)
        
        if 'm' in phase_name:
            phase2 = phase_name.replace('m', '')
            arrivals = model.get_ray_paths_geo(depth, elat, elon, slat, slon, [phase2])
            if not arrivals:
                pass
                save = False
                issue = 'major arc; no listed arrival'
            else:
                try:
                    for arrival in arrivals:
                        arr_dist = arrival.purist_distance
                        if arr_dist > 180.:
                            ellip_correction = ellipticity_correction(arrival, azimuth = path_az, source_latitude = slat)
                            pathfile = arrival.path
                            df_pathfile = pd.DataFrame(pathfile)
                            df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
                            save = True
                            break
                except:
                    save = False
                    issue = 'major arc; listed arrival/compute bug'
                    pass
    
        else:
            arrivals = model.get_ray_paths_geo(depth, elat, elon, slat, slon, [phase_name])
            if not arrivals:
                pass
                save = False
                issue = 'minor arc; no listed arrival'
            else:            
                arrival = arrivals[0]
                ellip_correction = ellipticity_correction(arrival, azimuth = path_az, source_latitude = slat)
                pathfile = arrival.path
                df_pathfile = pd.DataFrame(pathfile)
                df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
                save = True
    
        if save == True:
            out_path = f'./{phases_directory}/{phase}/{paths_directory}/orig_{phase_name}_{path_id}_{elat}_{elon}_{slat}_{slon}.csv'
            df_pathfile = df_pathfile.rename(columns = {'time': 'TIME', 'dist': 'DIST', 'depth': 'DEPTH', 'lat': 'LAT', 'lon': 'LON'})
            df_pathfile['SHELL#'] = 0
            df_pathfile['BLOCK#'] = 0
            df_pathfile['BOUND'] = 'hold'

            df_pathfile = np.array(df_pathfile)
            
            # find block#, shell#, and bound type for each point in the new pathfile
            for point in range(len(df_pathfile)):

                path_depth = df_pathfile[point, 3]
                path_lat = df_pathfile[point, 4]
                path_lon = df_pathfile[point, 5]
    
                if path_lon >= 180.:
                    path_lon -= 360.
    
                # add other identifying attributes
                shell = mod_database.find_shell_id(path_depth)
                block = mod_database.find_block_id(path_lat, path_lon)
                bound_type = mod_database.find_boundary_type(path_lat, path_lon, path_depth)

                df_pathfile[point, 6] = shell
                df_pathfile[point, 7] = block
                df_pathfile[point, 8] = bound_type            

            df_pathfile = pd.DataFrame(df_pathfile, columns = ['p', 'TIME', 'DIST', 'DEPTH', 'LAT', 'LON', 'SHELL#', 'BLOCK#', 'BOUND'])
            df_pathfile['SHELL#'] = df_pathfile['SHELL#'].astype(int)
            df_pathfile['BLOCK#'] = df_pathfile['BLOCK#'].astype(int)
            
            df_data.loc[path, 'DIST'] = df_pathfile['DIST'].iloc[-1]
            df_data.loc[path, 'ELLIP_CORR'] = round(ellip_correction, rdp)
            df_data.loc[path, 'ELLIP_DT'] = dt - (round(ellip_correction, rdp))
            df_data.loc[path, 'PATH_ID'] = int(path_id)

            df_pathfile.to_csv(out_path, index = False)
            path_id += 1
    
        elif save == False:
            df_data.loc[path, 'PATH_ID'] = np.nan
            with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_bugs_make_raypaths.txt', 'a') as fout:
                fout.write(f'{depth},{elat},{elon},{slat},{slon},{issue}\n')
                
    df_data = df_data.dropna(subset = ['PATH_ID']).reset_index(drop = True)
    df_data.to_csv(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_master_data.csv', index = False)
    
    with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_log_make_raypaths.txt', 'a') as fout:
        fout.write(f'FINISHED; runtime: {mod_track.runtime(time.time() - start_raypaths)}\n')

## MAKE RAYPATHS:
if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        process_idx += 1

        p = mp.Process(target = make_raypaths, args = (phase,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

#######################################################################
#######################################################################
##          DATA IS PARSED; ORIGINAL TAUP RAYPATHS ARE MADE          ##
#######################################################################
##                     START RESAMPLING RAYPATHS                     ##
#######################################################################
#######################################################################

cdp = mod_input.computed_decimal_places
rdp = mod_input.rounded_decimal_places
discontinuities = mod_input.discontinuities
depth_bounds = mod_input.shell_bounds
target_length = mod_input.target_path_length
target_length_tol = mod_input.target_path_length_tolerance
lower_lim = target_length - target_length_tol
upper_lim = target_length + target_length_tol

if mod_input.data_wave_type == 'S':
    df_crust = pd.read_csv(f'./{mod_input.tomography_model_directory}/crust/CRUST_1.0_vsh.csv')
    vel_label = 'vs'
if mod_input.data_wave_type == 'P':
    df_crust = pd.read_csv(f'./{mod_input.tomography_model_directory}/crust/CRUST_1.0_vp.csv')
    vel_label = 'vp'

crust_cols = list(df_crust)
crust_layers = []
for col in crust_cols:
    if '_thickness' in col:
        crust_layers.append('_'.join(col.split('_')[:-1]))

def convert_velocity(seismic_velocity, depth):
    '''
    function to convert velocity (`seismic_velocity`) from km/s to deg/s.
    function is total radius + depth (rather than minus) because CRUST1.0 gives depth as negative values.
    inputs:
        - `seismic_velocity`: seismic velocity in terms of km/s
        - `depth`: the depth at which the seismic velocity is computed.
    
    returns:
        - seismic velocity in terms of degrees/second, with the length of a degree dependent on the passed `depth`, which is used to determine the radius of the sphere at which the length of a degree is calculated.
    '''
    return seismic_velocity / ((2 * np.pi * (mod_input.total_radius + depth)) / 360)

prem_vel_upper_crust = convert_velocity(mod_refmodels.prem_vel(mod_input.data_wave_type, 0.), 0.)
prem_vel_lower_crust = convert_velocity(mod_refmodels.prem_vel(mod_input.data_wave_type, 15.), 15.)

def prem_crust_travel_time(ray_param, depth_max):
    if depth_max <= 15.:
        # find the time through the upper crust layer
        theta_uc = np.arcsin(ray_param * prem_vel_upper_crust)
        d_uc = depth_max / np.cos(theta_uc)
        t_uc = d_uc / mod_refmodels.prem_vel(mod_input.data_wave_type, 0.)
        prem_crust_time = t_uc

    elif depth_max > 15. and depth_max <= 24.4:
        # find the time through the upper crust layer
        theta_uc = np.arcsin(ray_param * prem_vel_upper_crust)
        d_uc = 15. / np.cos(theta_uc)
        t_uc = d_uc / mod_refmodels.prem_vel(mod_input.data_wave_type, 0.)
        
        # find the time through the lower crust layer
        theta_lc = np.arcsin(ray_param * prem_vel_lower_crust)
        d_lc = (depth_max - 15.) / np.cos(theta_lc)
        t_lc = d_lc / mod_refmodels.prem_vel(mod_input.data_wave_type, 15.)
        
        # add the times from the two layers together
        prem_crust_time = (t_uc + t_lc)

    elif depth_max > 24.4:
        # find the time through the upper crust layer
        theta_uc = np.arcsin(ray_param * prem_vel_upper_crust)
        d_uc = 15. / np.cos(theta_uc)
        t_uc = d_uc / mod_refmodels.prem_vel(mod_input.data_wave_type, 0.)
        
        # find the time through the lower crust layer
        theta_lc = np.arcsin(ray_param * prem_vel_lower_crust)
        d_lc = 9.4 / np.cos(theta_lc)
        t_lc = d_lc / mod_refmodels.prem_vel(mod_input.data_wave_type, 15.)
        
        # find the time through the layer that enters the mantle
        m_thickness = depth_max - 24.4
        m_mid_depth = (m_thickness / 2) + 24.4
        m_prem_vel = mod_refmodels.prem_vel(mod_input.data_wave_type, m_mid_depth)
        theta_m = np.arcsin(ray_param * convert_velocity(m_prem_vel, m_mid_depth))
        d_m = m_thickness / np.cos(theta_m)
        t_m = d_m / m_prem_vel
        
        # add the times from all of the layers together
        prem_crust_time = (t_uc + t_lc + t_m)
        
    return prem_crust_time

# boundary finding code/function for each parallel process:
def find_boundaries_resample(phase):
    phase_name = phase.split('_')[0]
    ## Define functions ##
    def add_pt_info(path_len, mid_depth, p_azimuth, p_depth, p_dist, p_seg_dist, p_lat, p_lon, p_shell, p_block):
        km_dist = path_len + cumulative_lengths[-1]
        physical_property = mod_refmodels.prem_vel(mod_input.data_wave_type, mid_depth)
        seg_time = path_len / physical_property
        time = seg_time + cumulative_times[-1]
        p_bound_type = mod_database.find_boundary_type(p_lat, p_lon, p_depth)
        mid_pt = mod_geo.pt_from_dist(p_azimuth, pt_lats[-1], pt_lons[-1], (p_seg_dist / 2))
        seg_shell = mod_database.find_shell_id(mid_depth)
        seg_block = mod_database.find_block_id(mid_pt[0], mid_pt[1])
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

    try:
        os.mkdir(f'./{mod_input.phases_directory}/{phase}/{mod_input.resampled_directory}')
    except:
        pass
    try:
        os.remove(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_boundary_finding.txt')
    except:
        pass
    
    start_boundary_finding_time = time.time()
    ## Start code ##
    df_phase_data = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_master_data.csv')
    paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.paths_directory}/orig_{phase_name}_*.csv')
    tot_paths = len(paths)
    itr = 0

    # begin master loop to loop through all raypaths
    for path in paths:
        itr += 1
        try:
            df_p = pd.read_csv(path)
            path_id = int(path.split('_')[-5])
            ray_param = df_p['p'][1]
            df_p = df_p.drop(columns = ['p', 'TIME'])
            max_depth = df_p['DEPTH'].max()
            with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_boundary_finding.txt', 'a') as fout:
                fout.write(f'- working on path {itr} (path id: {path_id}) of {len(paths)}\n')

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
                    pi_shell = mod_database.find_shell_id(pi_depth)
                    pi_block = mod_database.find_block_id(pi_lat, pi_lon)
                    pi_bound_type = mod_database.find_boundary_type(pi_lat, pi_lon, pi_depth)
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
            df_all = pd.concat([df_bounds, df_p], sort = False).sort_values(by = ['DIST']).reset_index(drop = True)

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
    
            cumulative_epidists.append(0.)
            cumulative_lengths.append(0.)
            cumulative_times.append(0.)

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
        
            df_resamp = pd.DataFrame(data = {'RAY_PARAM': p, 'START_DEPTH': start_depths, 'START_DIST': start_epidists, 'START_LAT': start_lats, 'START_LON': start_lons, 'START_SHELL#': start_shell_nos, 'START_BLOCK#': start_block_nos, 'START_TYPE': start_bound_types, 'END_DEPTH': end_depths, 'END_DIST': end_epidists, 'END_LAT': end_lats, 'END_LON': end_lons, 'END_SHELL#': end_shell_nos, 'END_BLOCK#': end_block_nos, 'END_TYPE': end_bound_types, 'AZIMUTH': azimuths, 'SEG_DIST_DEG': segment_epidists, 'SEG_DIST_KM': segment_lengths, 'SEG_TIME': segment_times, 'SEG_SHELL#': segment_shell_nos, 'SEG_BLOCK#': segment_block_nos, 'MID_DEPTH': segment_mid_depths, f'SEG_V{mod_input.data_wave_type}': segment_properties, 'LENGTH_PERCENT': percent_pathlens, 'TIME_PERCENT': percent_times, 'TIME': cumulative_times, 'DIST_DEG': cumulative_epidists, 'DIST_KM': cumulative_lengths})

            ##### now, tragically, after all that, it's time to add crustal corrections.
            df_resampled_path_crust = df_resamp.loc[df_resamp['END_DEPTH'] == 0.].copy().reset_index(drop = True)
            ray_param_deg = (ray_param / (180. / np.pi))
            
            prem_crustal_model_time = 0.
            crustal_model_time = 0.
            # find the total amount of time that the path leading up to the surface point(s) spent traveling through CRUST1.0.
            for resampled_seg in range(len(df_resampled_path_crust)):
                surface_lat = df_resampled_path_crust['END_LAT'].iloc[resampled_seg]
                surface_lon = df_resampled_path_crust['END_LON'].iloc[resampled_seg]
                surface_dist = df_resampled_path_crust['DIST_DEG'].iloc[resampled_seg]
                
                if surface_lat == 90.:
                    df_crust_block = df_crust.loc[(df_crust['lat_max'] == 90.) & (df_crust['lon_min'] <= surface_lon) & (df_crust['lon_max'] > surface_lon)].copy().reset_index(drop = True)
                else:
                    df_crust_block = df_crust.loc[(df_crust['lat_min'] <= surface_lat) & (df_crust['lat_max'] > surface_lat) & (df_crust['lon_min'] <= surface_lon) & (df_crust['lon_max'] > surface_lon)].copy().reset_index(drop = True)
                crust_block_depth = abs(df_crust_block['mantle_top'].iloc[0])

                # find the total travel time of prem through the crust block:
                prem_crustal_model_time += prem_crust_travel_time(ray_param_deg, crust_block_depth)
                # find the total path travel time through CRUST1.0
                for layer in crust_layers:
                    layer_thickness = df_crust_block[f'{layer}_thickness'].iloc[0]
                    layer_vel = df_crust_block[f'{layer}_{vel_label}'].iloc[0]
                    if layer_thickness != 0. and layer_vel != 0.:
                        layer_top = df_crust_block[f'{layer}_top'].iloc[0]
                        layer_theta = np.arcsin(ray_param_deg * convert_velocity(layer_vel, layer_top))
                        layer_d = layer_thickness / np.cos(layer_theta)
                        layer_t = layer_d / layer_vel

                        if surface_dist == df_resamp['DIST_DEG'].iloc[-1]:
                            crustal_model_time += layer_t
                        else:
                            crustal_model_time += (layer_t * 2)

            # one last check to see if the starting point is in the crust:
            path_start_depth = df_resamp['START_DEPTH'].iloc[0]
            path_start_lat = df_resamp['START_LAT'].iloc[0]
            path_start_lon = df_resamp['START_LON'].iloc[0]
            if path_start_lat == 90.:
                df_crust_block = df_crust.loc[(df_crust['lat_max'] == 90.) & (df_crust['lon_min'] <= path_start_lon) & (df_crust['lon_max'] > path_start_lon)].copy().reset_index(drop = True)
            else:
                df_crust_block = df_crust.loc[(df_crust['lat_min'] <= path_start_lat) & (df_crust['lat_max'] > path_start_lat) & (df_crust['lon_min'] <= path_start_lon) & (df_crust['lon_max'] > path_start_lon)].copy().reset_index(drop = True)
            crust_block_depth = abs(df_crust_block['mantle_top'].iloc[0])
            
            if path_start_depth >= crust_block_depth:
                pass
            else:
                # find the total travel time of prem through the crust block:
                prem_crustal_model_time += prem_crust_travel_time(ray_param_deg, crust_block_depth)
                # find the total path travel time through CRUST1.0
                for layer in crust_layers:
                    layer_thickness = df_crust_block[f'{layer}_thickness'].iloc[0]
                    layer_vel = df_crust_block[f'{layer}_{vel_label}'].iloc[0]
                    if layer_thickness != 0. and layer_vel != 0.:
                        layer_top = df_crust_block[f'{layer}_top'].iloc[0]
                        layer_bottom = df_crust_block[f'{layer}_top'].iloc[0] - layer_thickness
                        # if the EQ starts in this layer, start the calc. here:
                        if layer_top >= -path_start_depth > layer_bottom:
                            layer_theta = np.arcsin(ray_param_deg * convert_velocity(layer_vel, layer_top))
                            layer_d = (abs(layer_bottom) - path_start_depth) / np.cos(layer_theta)
                            layer_t = layer_d / layer_vel
                            crustal_model_time += layer_t
                        elif -path_start_depth > layer_top > layer_bottom:
                            layer_top = df_crust_block[f'{layer}_top'].iloc[0]
                            layer_theta = np.arcsin(ray_param_deg * convert_velocity(layer_vel, layer_top))
                            layer_d = layer_thickness / np.cos(layer_theta)
                            layer_t = layer_d / layer_vel
                            crustal_model_time += layer_t

            # then, FINALLY, add this information to the master_data.csv file.
            master_idx = df_phase_data.loc[df_phase_data['PATH_ID'] == path_id].index
            df_phase_data.loc[master_idx, 'CRUST_1.0_CORR'] = crustal_model_time - prem_crustal_model_time

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
            np.savetxt(f'./{mod_input.phases_directory}/{phase}/{mod_input.resampled_directory}/{phase_name}_{path_id}_resampled_segments.csv', df_resamp, fmt = format_string, header = headers, delimiter = ',', comments = '')
            # print(f'Path number: {itr} of {tot_paths}; path: {path}')
    
        except:
            # print(f'Path number: {itr} of {tot_paths}; path: {path}; ***** BUG *****')
            # save the name of the path to a new file
            with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_bugs_boundary_finding.txt', 'a') as fout:
                fout.write(f'{path}\n')
        
        # now that the path has been resampled, remove the original taup path file:
        os.remove(path)

    # do this at the very very end with the whole thing finished
    df_phase_data['CRUST_1.0_DT'] = df_phase_data['DT'] - df_phase_data['CRUST_1.0_CORR']
    df_phase_data['CRUST_1.0_ELLIP_DT'] = df_phase_data['DT'] - df_phase_data['CRUST_1.0_CORR'] - df_phase_data['ELLIP_CORR']
    df_phase_data.to_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_master_data.csv', index = False)
    shutil.rmtree(f'./{mod_input.phases_directory}/{phase}/{mod_input.paths_directory}')

    with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_boundary_finding.txt', 'a') as fout:
        fout.write(f'FINISHED; total time: {mod_track.runtime(time.time() - start_boundary_finding_time)}\n')

## Start parallel processes:
if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        process_idx += 1

        p = mp.Process(target = find_boundaries_resample, args = (phase,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()


########################################################################
########################################################################
##                    FINISHED RESAMPLING RAYPATHS                    ##
########################################################################
##                      START COMPUTING COVERAGE                      ##
########################################################################
########################################################################

df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
df_all = pd.DataFrame()
all_shells = df_shells['SHELL#']
all_blocks = df_blocks['BLOCK#']
total_phases = len(mod_input.all_phases)
azimuthal_sector_extent = 180. / mod_input.azimuthal_sectors

df_sectors = pd.DataFrame(columns = ['sector', 'min_extent', 'max_extent'])
sector_nos = list(range(1, mod_input.azimuthal_sectors + 1))
df_sectors['sector'] = sector_nos
df_sectors['min_extent'] = (df_sectors['sector'] - 1) * azimuthal_sector_extent
df_sectors['max_extent'] = df_sectors['sector'] * azimuthal_sector_extent
df_sectors = np.array(df_sectors.T)

def azimuthal_sector(az):
    '''
    define a function to determin the azimuthal sector of an input azimuth
    '''
    if az == 360.:
        az = 0.
    if az >= 180.:
        az -= 180.
    return int(np.where((df_sectors[1] <= az) & (df_sectors[2] > az))[0][0] + 1)

column_names = ['SHELL#', 'BLOCK#', 'TOTAL_PATHS', 'TOTAL_SECTORS']
sh_col_all = []
bl_col_all = []
for sh_all in all_shells:
    for bl_all in all_blocks:
        sh_col_all.append(sh_all)
        bl_col_all.append(bl_all)
df_all['SHELL#'] = sh_col_all
df_all['BLOCK#'] = bl_col_all
df_all['TOTAL_PATHS'] = 0
df_all['TOTAL_SECTORS'] = 0

for az_sector in df_sectors[0]:
    sector_header = f'SECTOR_{int(az_sector)}'
    column_names.append(sector_header)
    df_all[sector_header] = 0

try:
    os.mkdir(f'./coverage')
except:
    pass

def find_phase_coverage(phase):
    phase_name = phase.split('_')[0]
    try:
        os.remove(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_coverage.txt')
    except:
        pass
        
    start_coverage_time = time.time()
    paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.resampled_directory}/{phase_name}_*.csv')
    df_phase = pd.DataFrame()
    
    # make the master dataframe for all of the sectors defined in mod_input. this will be the empty dataframe for each individual phase.
    sh_col = []
    bl_col = []
    for sh in all_shells:
        for bl in all_blocks:
            sh_col.append(sh)
            bl_col.append(bl)
    df_phase['SHELL#'] = sh_col
    df_phase['BLOCK#'] = bl_col
    df_phase['TOTAL_PATHS'] = 0
    df_phase['TOTAL_SECTORS'] = 0
    for az_sector in df_sectors[0]:
        df_phase[f'SECTOR_{int(az_sector)}'] = 0
    df_phase = np.array(df_phase)
    
    total_paths = len(paths)
    itr = 0
    
    # for each path file, find all of the unique block & shell id combinations that the ray passes through    
    for path in paths:
        itr += 1
        with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_coverage.txt', 'a') as fout:
            fout.write(f'- working on path {itr} of {len(paths)}\n')
    
        df_p = pd.read_csv(path)
        df_slice = df_p[['SEG_SHELL#', 'SEG_BLOCK#']]
        df_count = df_slice.groupby(['SEG_SHELL#','SEG_BLOCK#']).size().reset_index().rename(columns = {0: 'TOTAL_PATHS'})
        
        df_p = np.array(df_p)
        df_slice = np.array(df_slice)        
        df_count = np.array(df_count.T)
    
        # find the mean azimuth for each of those unique block & shell id combinations, and its corresponding azimuthal sector
        for line in range(len(df_count[0])):
            segment_s_no = df_count[0, line]
            segment_b_no = df_count[1, line]
            element_idx = np.where((df_phase.T[0] == segment_s_no) & (df_phase.T[1] == segment_b_no))[0]
            segment_a_indices = np.where((df_p.T[19] == segment_s_no) & (df_p.T[20] == segment_b_no))[0]
            segment_a = np.mean(df_p[segment_a_indices, 15])
            sector = azimuthal_sector(segment_a)
    
            
            # add one to the count of total paths per block/shell combination for the current combo in the phase specific dataframe
            df_phase[element_idx, 2] += 1
            
            # add one to the total sector count IF it's not already accounted for
            individual_phase_sector_count = df_phase[element_idx, int(sector + 3)][0]
    
            if individual_phase_sector_count == 0.:
                df_phase[element_idx, 3] += 1
            
            # add one to the count of the specific segment azimuthal sector for the current combo in the phase specific dataframe
            df_phase[element_idx, int(sector + 3)] += 1
    
    df_phase = pd.DataFrame(df_phase, columns = column_names)
    df_phase.to_csv(f'./coverage/{phase}_coverage.csv', index = False)
    
    with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_coverage.txt', 'a') as fout:
        fout.write(f'FINISHED; total time: {mod_track.runtime(time.time() - start_coverage_time)}\n')

## Start parallel processes:
coverage_start = time.time()

if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        phase_name = phase.split('_')[0]
        with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_coverage.txt', 'a') as fout:
            fout.write(f'START PARSING INDIVIDUAL PHASE COVERAGE\n')
        process_idx += 1

        p = mp.Process(target = find_phase_coverage, args = (phase,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

    coverage_time = time.time() - coverage_start
    with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase_name}_pt2_log_coverage.txt', 'a') as fout:
        fout.write(f'FINISHED PARSING INDIVIDUAL PHASE COVERAGE; runtime: {mod_track.runtime(coverage_time)}\n')

# loop through each newly made individual coverage file to compile the sums into the master/total coverage file.
for phase in mod_input.all_phases:
    df_phase_coverage = pd.read_csv(f'./coverage/{phase}_coverage.csv')
    df_all['TOTAL_PATHS'] = df_all['TOTAL_PATHS'] + df_phase_coverage['TOTAL_PATHS']

    for sector in sector_nos:
        df_all[f'SECTOR_{sector}'] = df_all[f'SECTOR_{sector}'] + df_phase_coverage[f'SECTOR_{sector}']

df_all = np.array(df_all)
# finally, loop through each line in df_all to get the final total sector counts
for line in range(len(df_all)):
    total_shell_no = df_all[line, 0]
    total_block_no = df_all[line, 1]
    sector_count = 0

    for sector in sector_nos:
        if df_all[line, int(sector + 3)] != 0.:
            sector_count += 1

    df_all[line, 3] = sector_count

df_all = pd.DataFrame(df_all, columns = column_names)
df_all.to_csv(f'./coverage/total_coverage_{mod_input.dataset}', index = False)

end_subj = f'Process ended (PID: {pid}); part2.data_processing.py'
end_text = f'Process {pid} complete;\nRuntime: {mod_track.runtime(time.time() - start_time)}'
try:
    mod_track.SendMsg(end_subj, end_text)
except:
    pass


