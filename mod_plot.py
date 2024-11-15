## Functions to plot info from approximate equal area grid (e.g., SITRUS part1.py/Stage 1 style
import mod_input
import numpy as np
import pandas as pd
import mod_database

all_shells = pd.read_csv(mod_input.shell_file)['SHELL#']
all_blocks = pd.read_csv(mod_input.block_file)['BLOCK#']

def shell_grid_register(df_to_register, column_name, plot_grid_lat_increment, plot_grid_lon_increment):
    '''
    function to convert a model on an approximate equal area grid (according to `all_blocks` and `all_shells`) to a plottable array
    '''
    ar_to_register = np.array(df_to_register[['SHELL#', 'BLOCK#', column_name]])
    reg_lats = np.arange(mod_input.start_lat, mod_input.final_lat, plot_grid_lat_increment)
    reg_lons = np.arange(mod_input.start_lon, mod_input.final_lon, plot_grid_lon_increment)
    ar = np.zeros([len(reg_lons), len(reg_lats)])
    lat_idx = 0
    for lat in reg_lats:
        lon_idx = 0
        for lon in reg_lons:
            block = mod_database.find_block_id(lat, lon)
            value = ar_to_register[np.where(ar_to_register.T[1] == block)[0][0], 2]
            ar[lon_idx, lat_idx] = value
            lon_idx += 1
        lat_idx += 1
    ar = ar.T
    return ar


def shell_coverage_mesh(phases, shell_to_plot, column_to_plot, plot_grid_lat_increment, plot_grid_lon_increment):
    reg_lats = np.arange(mod_input.start_lat, mod_input.final_lat, plot_grid_lat_increment)
    reg_lons = np.arange(mod_input.start_lon, mod_input.final_lon, plot_grid_lon_increment)
    ar = np.zeros([len(reg_lons), len(reg_lats)])

    ar_all = np.zeros([len(all_blocks), 3 + mod_input.azimuthal_sectors])
    ar_all[:,0] = all_blocks

    for phase in phases:
        ar_phase_coverage = np.loadtxt(f'./coverage/{phase}_coverage.csv', skiprows = 1, delimiter = ',')
        ar_phase_coverage = ar_phase_coverage[np.where(ar_phase_coverage.T[0] == shell_to_plot)[0]][:, 1:]
        ar_all[:, 1] = ar_all[:, 1] + ar_phase_coverage[:, 1]

        sector_ids = list(range(3, len(ar_phase_coverage.T)))
        for sector in sector_ids:
            ar_all[:, sector] = ar_all[:, sector] + ar_phase_coverage[:, sector]

    for block in range(len(ar_all)):
        sector_array = ar_all[block, 3:]
        sector_count = np.count_nonzero(sector_array)
        ar_all[block, 2] = sector_count

    if column_to_plot == 'TOTAL_PATHS':
        col_idx_to_plot = 1
    elif column_to_plot == 'TOTAL_SECTORS':
        col_idx_to_plot = 2
    elif 'SECTOR_' in column_to_plot:
        col_idx_to_plot = 2 + int(column_to_plot.split('_')[-1])

    lat_idx = 0
    for lat in reg_lats:
        lat_vals = []
        lon_idx = 0
        for lon in reg_lons:
            block = mod_database.find_block_id(lat, lon)
            value = ar_all[np.where(ar_all.T[0] == block)[0][0], col_idx_to_plot]
            ar[lon_idx, lat_idx] = value
            lon_idx += 1
        lat_idx += 1
    return ar.T

















