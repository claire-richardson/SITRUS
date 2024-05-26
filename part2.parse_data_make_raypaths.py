import os
import mod_geo
import fnmatch
import mod_input
import mod_pandas
import numpy as np
import pandas as pd
from obspy.taup import taup_geo
from obspy.taup import TauPyModel
from ellipticipy import ellipticity_correction

# directories:
phases_directory = mod_input.phases_directory
data_directory = mod_input.data_directory
orig_paths_directory = mod_input.orig_paths_directory
resampled_paths_directory = mod_input.resampled_paths_directory

## PARSE DATASET:
try:
    os.mkdir(phases_directory)
except:
    pass

df_data = pd.read_csv(mod_input.dataset)
df_phases = df_data['PHASE'].unique()
all_phases = list(df_phases)

headers = mod_input.main_headers
for header in mod_input.raw_headers_to_keep:
    headers.append(header)

# make directories & subdirectories for each phase in master file; separate phase information
for p in all_phases:
    try:
        os.mkdir(f'./{phases_directory}/{p}')
        os.mkdir(f'./{phases_directory}/{p}/{data_directory}')
        os.mkdir(f'./{phases_directory}/{p}/{orig_paths_directory}')
        os.mkdir(f'./{phases_directory}/{p}/{resampled_paths_directory}')     
        df_specific_phase = df_data[headers].copy()
        df_specific_phase['DIST'] = 0.
        df_specific_phase['ELLIP_CORR'] = 0.
        df_specific_phase['ELLIP_DT'] = 0.
        df_specific_phase['CRUST_1.0_CORR'] = 0.
        df_specific_phase['CRUST_1.0_DT'] = 0.
        df_specific_phase['CRUST_1.0_ELLIP_DT'] = 0.
        df_specific_phase['PATH_ID'] = 0
        df_specific_phase.to_csv(f'./{phases_directory}/{p}/{data_directory}/{p}_master_data.csv', index = False)
    except:
        pass

## MAKE RAYPATHS:
# define variables and files to be used:
model = TauPyModel(model = mod_input.reference_model)

print(f'START MAKING RAYPATHS')
make_raypaths_start = time.time()

if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        process_idx += 1
        print(f'  - starting process {process_idx} of {len(mod_input.all_phases)}')

        p = mp.Process(target = make_raypaths, args = (, ,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

make_raypaths_time = time.time() - make_raypaths_start
print(f'FINISHED MAKING RAYPATHS; runtime: {backmapping_time / 60} minutes / {(backmapping_time / 60) / 60} hours')











def make_raypaths(phase):
    # parse phase name for multi-bounce and major arc indicators so it can be passed to TauP
    major = False
    mult = None
    for letter in phase:
        try:
            mult = int(letter)
        except:
            pass
        if letter == 'm':
            major = True
    # if the phase is multi-bounce and major arc:
    if mult != None and major == True:
        phase_name = (phase.split(str(mult))[0] * mult) + 'm'
    # if the phase is not multi-bounce but is major arc
    elif mult != None and major == False:
        phase_name = (phase.split(str(mult))[0] * mult)
    else:
        phase_name = phase
    with open(f'./{data_directory}/{phase}/{data_directory}/{phase}_pathfile_bugs.txt', 'w') as fout:
        fout.write(f'EQ_DEPTH,EQ_LAT,EQ_LON,STA_LAT,STA_LON,BUG\n')
    
    df_data = pd.read_csv(f'./{phases_directory}/{phase}/{data_directory}/{phase}_master_data.csv')
    path_id = 0
    for path in range(len(df_data)):
        path_2 = path + 1
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
                    path_nos = []
                    path_dists = []
    
                    for a in range(len(arrivals)):
                        path_nos.append(a)
                        arrival = arrivals[a]
                        pathfile = arrival.path[0:]
                        df_pathfile = pd.DataFrame(np.array(pathfile))
                        path_dist = df_pathfile['dist'].iloc[-1] * (180. / np.pi)
                        path_dists.append(path_dist)
    
                    major_arcs = []
                    for di in path_dists:
                        if di > 180.:
                            major_arcs.append(di)
    
                    idx = path_dists.index(major_arcs[0])
    
                    arrival = arrivals[idx]
                    ellip_correction = ellipticity_correction(arrival, azimuth = path_az, source_latitude = slat)
                    pathfile = arrival.path[0:]
                    df_pathfile = pd.DataFrame(np.array(pathfile))
                    df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
                    save = True
                
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
                pathfile = arrival.path[0:]
                df_pathfile = pd.DataFrame(np.array(pathfile))
                df_pathfile['dist'] = df_pathfile['dist'] * (180. / np.pi)
                save = True

        if save == True:
            out_path = f'./{phases_directory}/{phase}/{orig_paths_directory}/{phase}_{path_id}_{elat}_{elon}_{slat}_{slon}.csv'
            df_pathfile = df_pathfile.rename(columns = {'time': 'TIME', 'dist': 'DIST', 'depth': 'DEPTH', 'lat': 'LAT', 'lon': 'LON'})
            df_pathfile['SHELL#'] = 0
            df_pathfile['BLOCK#'] = 0
            df_pathfile['BOUND'] = 'hold'
            
            # find block#, shell#, and bound type for each point in the new pathfile
            for point in range(len(df_pathfile)):
                path_depth = df_pathfile['depth'].iloc[point]
                path_lat = df_pathfile['lat'].iloc[point]
                path_lon = df_pathfile['lon'].iloc[point]
    
                if path_lon >= 180.:
                    path_lon -= 360.

                shell = mod_pandas.find_shell_id(path_depth)
                block = mod_pandas.find_block_id(path_lat, path_lon)
                bound_type = mod_pandas.find_boundary_type(path_lat, path_lon, path_depth)

                df_pathfile.loc[point, 'SHELL#'] = shell
                df_pathfile.loc[point, 'BLOCK#'] = block
                df_pathfile.loc[point, 'BOUND'] = bound_type

            df_data.loc[path, 'DIST'] = df_pathfile['DIST'].iloc[-1]
            df_data.loc[path, 'ELLIP_CORR'] = round(ellip_correction, rdp)
            df_data.loc[path, 'ELLIP_DT'] = dt - (round(ellip_correction, rdp))
            df_data.loc[path, 'CRUST_1.0_CORR'] = 0.
            df_data.loc[path, 'CRUST_1.0_DT'] = 0.
            df_data.loc[path, 'CRUST_1.0_ELLIP_DT'] = 0.
            df_data.loc[path, 'PATH_ID'] = int(path_id)

            df_pathfile.to_csv(out_path, index = False)
            path_id += 1

        elif save == False:
            df_data.loc[path, 'PATH_ID'] = np.nan
            with open(f'{data_path}{phase}{outfile_mod}_pathfile_bugs.txt', 'a') as fout:
                fout.write(f'{depth},{elat},{elon},{slat},{slon},{issue}\n')

    df_data = df_data.dropna(subset = ['PATH_ID']).reset_index(drop = True)
    df_data.to_csv(f'./{phases_directory}/{phase}/{data_directory}/{phase}_master_data.csv', index = False)







