import os
import time
import mod_geo
import fnmatch
import mod_input
import mod_database
import numpy as np
import pandas as pd
import mod_refmodels
import multiprocessing as mp
from obspy.taup import taup_geo
from obspy.taup import TauPyModel
from ellipticipy import ellipticity_correction

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
        os.remove(f'./{data_directory}/{phase}/{data_directory}/{phase_name}_pt2_pathfile_bugs.txt')
    except:
        pass
    start_raypaths = time.time()
    
    with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_pathfile_bugs.txt', 'w') as fout:
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
            with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_pathfile_bugs.txt', 'a') as fout:
                fout.write(f'{depth},{elat},{elon},{slat},{slon},{issue}\n')
                
    df_data = df_data.dropna(subset = ['PATH_ID']).reset_index(drop = True)
    df_data.to_csv(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_master_data.csv', index = False)
    
    with open(f'./{phases_directory}/{phase}/{data_directory}/{phase_name}_pt2_log_make_raypaths.txt', 'a') as fout:
        fout.write(f'FINISHED; total time: {(time.time() - start_raypaths) / 60} minutes; {((time.time() - start_raypaths) / 60) / 60} hours\n')

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
