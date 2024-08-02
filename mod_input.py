########################
#
# MASTER INPUT PARAMETER FILE FOR SITRUS
#
# READ ME:
# 
# This is the master input parameter file that controls the output for all components
# of the SITRUS software, including the file structure, mesh definition, data and
# model preprocessing, and model updates themselves. While all variables in this
# module are user-modifiable, many should not be changed to ensure the functionality
# of the software. These are labeled in the section "NON-MODIFIABLE VARIABLES". If
# you do modify these variables, proceed at your own risk: you must ensure
# compatibility with all other SITRUS components, including all five modules and
# all eight scripts.
#
# Please refer to the *.readme.md file(s) for a list and description of the
# parameters that need to be set for a given SITRUS script.
#
########################
########################
# MODIFIABLE VARIABLES #
########################

# PROGRESS TRACKING # 
user_email = 'crricha5@asu.edu'

# PHASE AND DATA NAMING #
#['S_G', 'ScS_G', 'SS_G', 'SSS_G', 'SSSS_G', 'ScSScS_G', 'SSSSS_G', 'sS_G', 'sSS_G', 'sScSScS_G', 'sScS_G', 'sSSS_G', 'ScSScSScS_G', 'sSSSS_G', 'sScSScSScS_G']
#['SSS_H_vs', 'S_H_vs', 'ScSScS_H_vs', 'ScSScSScSScS_H_vs', 'ScS_H_vs', 'ScSScSScS_H_vs', 'SSSSm_H_vs', 'SSm_H_vs', 'Sdiff_H_vs', 'SSSm_H_vs', 'SSSSS_H_vs', 'ScSScSScSScSScSm_H_vs', 'SSSS_H_vs', 'SS_H_vs', 'SSSSSm_H_vs', 'ScSScSScSm_H_vs', 'ScSScSScSScSScS_H_vs', 'SSSSSSm_H_vs', 'ScSScSScSScSm_H_vs']
#['SSSm_H_mb', 'ScSScSScS_H_mb', 'ScSScSScSScS_H_mb', 'SSSSm_H_mb', 'SSSSSm_H_mb', 'ScSScSScSScSScS_H_mb', 'SSSSSSm_H_mb', 'ScSScSScSScSm_H_mb']
#['ScSScS_H', 'ScS_H', 'Sdiff_H', 'SS_H', 'S_H', 'SSS_H']
all_phases = ['S_R', 'SS_R', 'SSS_R', 'SSSm_R', 'ScSScS_R', 'ScSScSScS_R', 'ScS_R', 'SSm_R', 'SSSSm_R', 'sSS_R', 'sScSScS_R', 'sS_R', 'sSSSSm_R', 'sSSSm_R', 'sScSScSScS_R', 'sScS_R', 'sSSS_R', 'sSSm_R']
#['S_G', 'ScS_G', 'SS_G', 'SSS_G', 'SSSS_G', 'ScSScS_G', 'SSSSS_G', 'sS_G', 'sSS_G', 'sScSScS_G', 'sScS_G', 'sSSS_G', 'ScSScSScS_G', 'sSSSS_G', 'sScSScSScS_G']


# INPUT DATA FILES #
dataset = 'data_ritsema_clean.csv'
raw_headers_to_keep = ['EQ_DATE', 'STA_NAME', 'COMPREHENSIVE_WEIGHT', 'HIPR_BIN', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT']
#['EQ_DATE', 'STNM', 'STA_ELV', 'COMPREHENSIVE_WEIGHT', 'HIPR_BIN', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT']
#['EQ_DATE', 'STNM', 'COMPREHENSIVE_WEIGHT', 'HIPR_BIN', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT']
data_wave_type = 'S' # either 'S' or 'P'


# REFERENCE MODEL #
reference_model = 'prem' # prem, (ak135, iasp91 eventually?)
total_radius = 6371. # km
core_mantle_boundary = 2891. # km


# 3D MESH DEFINITION #
shell_bounds = [0., 24.4, 80., 160., 220., 310., 400., 490., 580., 670., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2891.] # km
discontinuities = [0., 15., 24.4, 220., 400., 600., 670., 771., 2741., 2891.] # km; should include discontinuities of the reference model
cardinal_azimuths = [0., 90., 180., 270., 360.] # degrees
make_near_neighbors_files = False
max_near_neighbors_radius = 15. # degrees
reference_lat = 2 # degrees
reference_lon = 2 # degrees
start_lat = -90  # degrees; starting latitude for your coordinate system
final_lat = 90  # degrees; ending latitude for your coordinate system
start_lon = -180  # degrees; starting longitude for your coordinate system
final_lon = 180 # degrees; ending longitude for your coordinate system


# PATH RESAMPLING PARAMETERS #
target_path_length = 80 # km
target_path_length_tolerance = 5 # km


# COVERAGE PARAMETERS #
azimuthal_sectors = 6


# MODEL PREPROCESSING
convert_nc_to_csv = False
zero_shells = [1]


## INPUT MODEL PARAMETERES #
input_model = 'S40RTS_vsh' # name of input model file
delimiter = ',' # delimiter in converted model.csv file
lat_header = 'latitude' # header in converted model.csv file
lon_header = 'longitude' # header in converted model.csv file
depth_header = 'depth' # header in converted model.csv file
property_header = 'vsh' # header in converted model.csv file
voigt = [False] #[False] [True, 'vpv', 'vph']
convert_vel_to_perturb = True


## MODEL UPDATE PARAMETERS
dataset_description = ['Grand+Ritsema+Hongyux3 (all)', 'Grand+Ritsema+Hongyux3 (all)'] # DONT USE COMMAS!
update_phases = [['ScSScS_H', 'ScS_H', 'Sdiff_H', 'SS_H', 'S_H', 'SSS_H', 'SSSm_H_mb', 'ScSScSScS_H_mb', 'ScSScSScSScS_H_mb', 'SSSSm_H_mb', 'SSSSSm_H_mb', 'ScSScSScSScSScS_H_mb', 'SSSSSSm_H_mb', 'ScSScSScSScSm_H_mb', 'SSS_H_vs', 'S_H_vs', 'ScSScS_H_vs', 'ScSScSScSScS_H_vs', 'ScS_H_vs', 'ScSScSScS_H_vs', 'SSSSm_H_vs', 'SSm_H_vs', 'Sdiff_H_vs', 'SSSm_H_vs', 'SSSSS_H_vs', 'ScSScSScSScSScSm_H_vs', 'SSSS_H_vs', 'SS_H_vs', 'SSSSSm_H_vs', 'ScSScSScSm_H_vs', 'ScSScSScSScSScS_H_vs', 'SSSSSSm_H_vs', 'ScSScSScSScSm_H_vs', 'S_R', 'SS_R', 'SSS_R', 'SSSm_R', 'ScSScS_R', 'ScSScSScS_R', 'ScS_R', 'SSm_R', 'SSSSm_R', 'sSS_R', 'sScSScS_R', 'sS_R', 'sSSSSm_R', 'sSSSm_R', 'sScSScSScS_R', 'sScS_R', 'sSSS_R', 'sSSm_R', 'S_G', 'ScS_G', 'SS_G', 'SSS_G', 'SSSS_G', 'ScSScS_G', 'SSSSS_G', 'sS_G', 'sSS_G', 'sScSScS_G', 'sScS_G', 'sSSS_G', 'ScSScSScS_G', 'sSSSS_G', 'sScSScSScS_G'], ['ScSScS_H', 'ScS_H', 'Sdiff_H', 'SS_H', 'S_H', 'SSS_H', 'SSSm_H_mb', 'ScSScSScS_H_mb', 'ScSScSScSScS_H_mb', 'SSSSm_H_mb', 'SSSSSm_H_mb', 'ScSScSScSScSScS_H_mb', 'SSSSSSm_H_mb', 'ScSScSScSScSm_H_mb', 'SSS_H_vs', 'S_H_vs', 'ScSScS_H_vs', 'ScSScSScSScS_H_vs', 'ScS_H_vs', 'ScSScSScS_H_vs', 'SSSSm_H_vs', 'SSm_H_vs', 'Sdiff_H_vs', 'SSSm_H_vs', 'SSSSS_H_vs', 'ScSScSScSScSScSm_H_vs', 'SSSS_H_vs', 'SS_H_vs', 'SSSSSm_H_vs', 'ScSScSScSm_H_vs', 'ScSScSScSScSScS_H_vs', 'SSSSSSm_H_vs', 'ScSScSScSScSm_H_vs', 'S_R', 'SS_R', 'SSS_R', 'SSSm_R', 'ScSScS_R', 'ScSScSScS_R', 'ScS_R', 'SSm_R', 'SSSSm_R', 'sSS_R', 'sScSScS_R', 'sS_R', 'sSSSSm_R', 'sSSSm_R', 'sScSScSScS_R', 'sScS_R', 'sSSS_R', 'sSSm_R', 'S_G', 'ScS_G', 'SS_G', 'SSS_G', 'SSSS_G', 'ScSScS_G', 'SSSSS_G', 'sS_G', 'sSS_G', 'sScSScS_G', 'sScS_G', 'sSSS_G', 'ScSScSScS_G', 'sSSSS_G', 'sScSScSScS_G']] # list of lists of phases to include in the model update. include name modifier. the number and order of lists must match the number and order of layers.
type_of_phase_subselection = ['proportion', 'proportion'] #['proportion'] # 'proportion' or 'number'
subselection_of_phase_data_to_use = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
residual_header = ['CRUST_ELLIP_DT', 'CRUST_ELLIP_DT'] # column name of the residual that you want to use in the data file.
residual_limits = [[True, -20., 20.], [True, -20., 20.]] #[[False]] or [[True, lower_lim, upper_lim]] argument for whether or not there are limits imposed on residual values. # list of lists should be the same length as the number of layers
layer_top_shells = [2, 2] # list of the top-most shell(s) in a given layer(s)
layer_base_shells = [31, 31] # list of the bottom-most shell(s) in a given layer(s)
freeze_previous_layers = [True, True] # [None, True, False]
remove_residual_mean = [True, True]
starting_RMS_model_to_use = input_model
iteration_to_stop_RMS_weighting = [2, 2] # the iteration for each layer of `n` layers on which to stop weighting backmapped perturbations by RMS weighting. a list with `n` elements, either integers or 'None'.
cutoff_type = ['total iterations', 'total iterations'] # 'reduction' or 'total iterations'
cutoff = [5, 5] # if 'reduction', type == float; if 'total iterations', type == int
HPC_cores = 31


# SMOOTHING PARAMETERS #
smoothing_radii = [[8.], [3., 4., 5.]] # list of lists of smoothing radii for each main layer in the current update
total_required_paths = [0, 20] # this is the minimum number of paths that are required in the radius for smoothing
total_required_azimuths = [0, 2] # this is the minimum number of azimuths that are required in each azimuthal sector for smoothing
gaussian_cutoff_weight = [0.5, 0.5] # Value of the Gaussian at the smoothing cuttoff
azimuthal_weighting = [True, False] # True or False to include azimuthal weighting when computing the smoothed azimuthal pertubation mean.
path_length_weighting = [False, False]
special_weights = [['COMPREHENSIVE_WEIGHT', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT'], ['COMPREHENSIVE_WEIGHT', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT']]


############################
# NON-MODIFIABLE VARIABLES #
############################
# FILE STRUCTURE #
shell_file = 'shell_dimensions.csv'
block_file = 'block_dimensions.csv'
phases_directory = 'phases'
data_directory = 'phase_data'
backmapped_paths_directory = 'backmapped_path_files'
paths_directory = 'raypath_files'
resampled_directory = 'resampled_path_files'
tomography_model_directory = 'models'
near_neighbors_directory = 'near_neighbors'
block_centric_directory = 'model_block_information'
main_headers = ['PHASE', 'EQ_LAT', 'EQ_LON', 'STA_LAT', 'STA_LON', 'EQ_DEP', 'DT']


# rounding values
computed_decimal_places = 10
rounded_decimal_places = 5




