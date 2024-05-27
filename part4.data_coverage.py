import glob
import time
import mod_input
import numpy as np
import pandas as pd
import multiprocessing as mp

df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
df_all = pd.DataFrame()
all_shells = df_shells['SHELL#']
all_blocks = df_blocks['BLOCK#']
total_phases = len(mod_input.all_phases)
azimuthal_sector_extent = 180. / mod_input.azimuthal_sectors

df_sectors = pd.DataFrame(columns = ['sector', 'min_extent', 'max_extent'])
df_sectors['sector'] = list(range(1, mod_input.azimuthal_sectors + 1))
df_sectors['min_extent'] = (df_sectors['sector'] - 1) * azimuthal_sector_extent
df_sectors['max_extent'] = df_sectors['sector'] * azimuthal_sector_extent

def azimuthal_sector(az):
    '''
    define a function to determin the azimuthal sector of an input azimuth
    '''
    if az == 360.:
        az = 0.
    if az >= 180.:
        az -= 180.
    for line in range(len(df_sectors)):
        s = df_sectors['sector'].iloc[line]
        n = df_sectors['min_extent'].iloc[line]
        x = df_sectors['max_extent'].iloc[line]
        if az >= n and az < x:
            return int(s)

# make the master dataframe for all of the sectors defined in mod_input. this will be the empty dataframe for ALL phase information
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
for az_sector in df_sectors['sector']:
    df_all[f'SECTOR_{az_sector}'] = 0


def find_phase_coverage(phase):
    try:
        os.remove(f'./{phases_directory}/{phase}/{data_directory}/{phase}_pt4_coverage_log.txt')
    except:
        pass
    start_coverage_time = time.time()
    paths = glob.glob(f'./{mod_input.phases_directory}/{phase}/{mod_input.paths_directory}/{phase}_*.csv')
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
    for az_sector in sector_nos:
        df_phase[f'SECTOR_{az_sector}'] = 0
    
    total_paths = len(paths)
    itr = 0
    percent_completes = []
    
    # for each path file, find all of the unique block & shell id combinations that the ray passes through    
    for path in paths:
        itr += 1
        with open(f'./{phases_directory}/{phase}/{data_directory}/{phase}_pt4_coverage_log.txt', 'a') as fout:
            fout.write(f'- working on path {itr} of {len(paths)}\n')
        df_p = pd.read_csv(path)
        df_slice = df_p[['SEG_SHELL#', 'SEG_BLOCK#']]
        df_count = df_slice.groupby(['SHELL#','BLOCK#']).size().reset_index().rename(columns = {0: 'TOTAL_PATHS'})
        
        # find the mean azimuth for each of those unique block & shell id combinations, and its corresponding azimuthal sector
        for line in range(len(df_count)):
            segment_s_no = df_count['SHELL#'].iloc[line]
            segment_b_no = df_count['BLOCK#'].iloc[line]
            element_idx = df_phase.loc[(df_phase['SHELL#'] == segment_s_no) & (df_phase['BLOCK#'] == segment_b_no)].index
            df_az = df_p.loc[(df_p['SEG_SHELL#'] == segment_s_no) & (df_p['SEG_BLOCK#'] == segment_b_no)]
            segment_a = df_az['AZIMUTH'].mean()
            sector = azimuthal_sector(segment_a)
            
            # add one to the count of total paths per block/shell combination for the current combo in the phase specific dataframe
            df_phase.loc[element_idx, 'TOTAL_PATHS'] += 1

            # add one to the total sector count IF it's not already accounted for
            individual_phase_sector_count = int(df_phase.loc[(df_phase['SHELL#'] == segment_s_no) & (df_phase['BLOCK#'] == segment_b_no)][f'SECTOR_{sector}'])

            if individual_phase_sector_count == 0:
                df_phase.loc[element_idx, 'TOTAL_SECTORS'] += 1

            # add one to the count of the specific segment azimuthal sector for the current combo in the phase specific dataframe
            df_phase.loc[element_idx, f'SECTOR_{sector}'] += 1

    df_phase.to_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_coverage.csv', index = False)
    with open(f'./{phases_directory}/{phase}/{data_directory}/{phase}_pt4_coverage_log.txt', 'a') as fout:
        fout.write(f'FINISHED; total time: {(time.time() - start_coverage_time) / 60 minutes; {((time.time() - start_coverage_time) / 60) / 60} hours\n')


## Start parallel processes:
print(f'START PARSING COVERAGE')
coverage_start = time.time()

if __name__ == '__main__':
    process_list = []
    process_idx = 0
    for phase in mod_input.all_phases:
        process_idx += 1
        print(f'  - starting process {process_idx} of {len(mod_input.all_phases)}')

        p = mp.Process(target = find_phase_coverage, args = (phase,))
        p.start()
        process_list.append(p)

    for process in process_list:
        process.join()

coverage_time = time.time() - coverage_start
print(f'FINISHED PARSING COVERAGE; runtime: {coverage_time / 60} minutes / {(coverage_time / 60) / 60} hours')


# loop through each newly made individual coverage file to compile the sums into the master/total coverage file.
for phase in mod_input.all_phases:
    df_phase_coverage = pd.read_csv(f'./{mod_input.phases_directory}/{phase}/{mod_input.data_directory}/{phase}_coverage.csv')
    df_all['TOTAL_PATHS'] = df_all['TOTAL_PATHS'] + df_phase_coverage['TOTAL_PATHS']

    for sector in sector_nos:
        df_all[f'SECTOR_{sector}'] = df_all[f'SECTOR_{sector}'] + df_phase_coverage[f'SECTOR_{sector}']


shells = []
# finally, loop through each line in df_all to get the final total sector counts
for line in range(len(df_all)):
    total_shell_no = df_all['SHELL#'].iloc[line]
    total_block_no = df_all['BLOCK#'].iloc[line]
    total_element_idx = df_all.loc[(df_all['SHELL#'] == total_shell_no) & (df_all['BLOCK#'] == total_block_no)].index
    sector_count = 0
    if total_shell_no not in shells:
        shells.append(total_shell_no)
    for sector in sector_nos:
        if int(df_all.loc[total_element_idx, f'SECTOR_{sector}']) != 0:
            sector_count += 1

    df_all.loc[total_element_idx, 'TOTAL_SECTORS'] = sector_count

df_all.to_csv(f'./{mod_input.phases_directory}/total_coverage.csv', index = False)





