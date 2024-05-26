import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mod_input
import glob

sector_nos = []
min_extents = []
max_extents = []
sector_val = 0.

header_list = ['SHELL#', 'BLOCK#', 'TOTAL_PATHS', 'TOTAL_SECTORS']
for i in range(mod_input.azimuthal_sectors):
    sector_no = int(i + 1)
    sector_nos.append(sector_no)
    min_extents.append(sector_val)
    sector_val += mod_input.azimuthal_sector_extent
    max_extents.append(sector_val)
    header_name = f'SECTOR_{sector_no}'
    header_list.append(header_name)
    df_sectors = pd.DataFrame(data = {'sector': sector_nos, 'min_extent': min_extents, 'max_extent': max_extents})    
header_str = ','.join(header_list)

def azimuthal_sector(az):
    '''
    define a function to determin the azimuthal sector of an input azimuth
    '''
    if az >= 180.:
        az -= 180.
    for line in range(len(df_sectors)):
        s = df_sectors['sector'].iloc[line]
        n = df_sectors['min_extent'].iloc[line]
        x = df_sectors['max_extent'].iloc[line]
        if az >= n and az < x:
            return int(s)

df_shells = pd.read_csv(mod_input.shell_file)
df_blocks = pd.read_csv(mod_input.block_file)
df_all = pd.DataFrame()
all_shells = df_shells['SHELL#']
all_blocks = df_blocks['BLOCK#']
total_phases = len(mod_input.all_phases)

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
for az_sector in sector_nos:
    df_all[f'SECTOR_{az_sector}'] = 0

phases_completed = []
# find azimuthal coverage for all phases, by phase.
for phase in mod_input.all_phases:
    phases_completed.append(phase)
    print(f'working on {phase} coverage')
    paths = glob.glob(f'./{mod_input.phases_directory}/{phase}{mod_input.dataset_name_mod}/{mod_input.resampled_paths_directory}/{phase}_*.csv')
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
        percent_complete = int((itr / total_paths) * 100)
        if percent_complete not in percent_completes:
            percent_completes.append(percent_complete)
            print(f'-  computing phase specific coverage; (phase {len(phases_completed)} of {total_phases}): {percent_complete}% complete')
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

    np.savetxt(f'./{mod_input.phases_directory}/{phase}{mod_input.dataset_name_mod}/{mod_input.data_directory}/coverage_{phase}.csv', df_phase, fmt = '%i', delimiter = ',', header = header_str, comments = '')


# loop through each newly made individual coverage file to compile the sums into the master/total coverage file.
for phase in mod_input.all_phases:
    print(f'compiling files, phase: {phase}')
    df_phase_coverage = pd.read_csv(f'./{mod_input.phases_directory}/{phase}{mod_input.dataset_name_mod}/{mod_input.data_directory}/coverage_{phase}.csv')
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
        print(f'working on shell {total_shell_no}')
        shells.append(total_shell_no)
    for sector in sector_nos:
        if int(df_all.loc[total_element_idx, f'SECTOR_{sector}']) != 0:
            sector_count += 1

    df_all.loc[total_element_idx, 'TOTAL_SECTORS'] = sector_count

np.savetxt(f'coverage_all{mod_input.dataset_name_mod}.csv', df_all, fmt = '%i', delimiter = ',', header = header_str, comments = '')

