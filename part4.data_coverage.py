import os
import glob
import time
import mod_input
import numpy as np
import pandas as pd
import mod_database
import multiprocessing as mp

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
    try:
        os.remove(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_pt4_log_coverage.txt')
    except:
        pass
        
    start_coverage_time = time.time()
    paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.resampled_directory}/{phase}_*.csv')
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
        with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_pt4_log_coverage.txt', 'a') as fout:
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

    lats = np.arange(mod_input.start_lat, mod_input.final_lat, mod_input.reference_lat)
    lons = np.arange(mod_input.start_lon, mod_input.final_lon, mod_input.reference_lon)
    try:
        os.mkdir(f'./coverage/{phase}_plot_files')
    except:
        pass
        
    for shell_to_plot in all_shells:
        avg_phase_paths_file = f'./coverage/{phase}_plot_files/gridded_paths_shell_{shell_to_plot}_{mod_input.reference_lat}deg_lat_by_{mod_input.reference_lon}deg_lon.csv'
        avg_phase_sectors_file = f'./coverage/{phase}_plot_files/gridded_sectors_shell_{shell_to_plot}_{mod_input.reference_lat}deg_lat_by_{mod_input.reference_lon}deg_lon.csv'
    
        df_phase_shell = df_phase.loc[df_phase['SHELL#'] == shell_to_plot].copy()
    
        top_depth = mod_database.get_shell_info(shell_to_plot)[1]
        bottom_depth = mod_database.get_shell_info(shell_to_plot)[3]
    
        df_phase_paths_grid = pd.DataFrame(data = {'LON': lons})
        df_phase_sectors_grid = pd.DataFrame(data = {'LON': lons})
        
        # loop through all latitudes in the model space
        for lat in lats:
            # make empty lists to fill in the values for the current latitude band
            lat_phase_paths = []
            lat_phase_sectors = []
        
            # loop through all longitudes in the model space
            for lon in lons:
                # find the block that the current lat/lon pair falls in
                block = mod_database.find_block_id(lat, lon)
    
                phase_paths = int(df_phase_shell.loc[df_phase_shell['BLOCK#'] == block]['TOTAL_PATHS'])
                phase_sectors = int(df_phase_shell.loc[df_phase_shell['BLOCK#'] == block]['TOTAL_SECTORS'])
                
                lat_phase_paths.append(phase_paths)
                lat_phase_sectors.append(phase_sectors)
                
            df_phase_paths_grid[f'{lat}'] = lat_phase_paths
            df_phase_sectors_grid[f'{lat}'] = lat_phase_sectors
    
        df_phase_paths_grid.to_csv(avg_phase_paths_file, index = False)
        df_phase_sectors_grid.to_csv(avg_phase_sectors_file, index = False)
    
    with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_pt4_log_coverage.txt', 'a') as fout:
        fout.write(f'FINISHED; total time: {(time.time() - start_coverage_time) / 60} minutes; {((time.time() - start_coverage_time) / 60) / 60} hours\n')


## Start parallel processes:
coverage_start = time.time()

if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_pt4_log_coverage.txt', 'a') as fout:
            fout.write(f'START PARSING INDIVIDUAL PHASE COVERAGE\n')
        process_idx += 1
        print(f'  - starting process {process_idx} of {len(mod_input.all_phases)}')

        p = mp.Process(target = find_phase_coverage, args = (phase,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

    coverage_time = time.time() - coverage_start
    with open(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_pt4_log_coverage.txt', 'a') as fout:
        fout.write(f'FINISHED PARSING INDIVIDUAL PHASE COVERAGE; runtime: {coverage_time / 60} minutes / {(coverage_time / 60) / 60} hours\n\n')

# loop through each newly made individual coverage file to compile the sums into the master/total coverage file.
for phase in mod_input.all_phases:
    df_phase_coverage = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_coverage.csv')
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
df_all.to_csv(f'./coverage/total_coverage.csv', index = False)

