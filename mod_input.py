##########
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
##########

########################
# MODIFIABLE VARIABLES #
########################

# PHASE AND DATA NAMING #
phase = 'SS'
all_phases = ['S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff'] # ['SSSm', 'ScSScSScS', 'ScSScSScSScS', 'SSSSm', 'SSSSSm', 'ScSScSScSScSScS', 'SSSSSSm', 'ScSScSScSScSm'] # ['S3', 'ScS', 'ScSScS', 'Sdiff', 'S', 'SSm', 'ScS4', 'ScS3', 'S4m', 'S3m', 'ScS5m', 'S5m', 'ScS3m', 'ScS5', 'S6m', 'ScS4m', 'SS', 'S4', 'S5'] # 
# dataset_name_mod = '' # '_mb' # '_vs' # 
all_phases_full = ['S3_vs', 'ScS_vs', 'ScSScS_vs', 'Sdiff_vs', 'S_vs', 'SSm_vs', 'ScS4_vs', 'ScS3_vs', 'S4m_vs', 'S3m_vs', 'ScS5m_vs', 'S5m_vs', 'ScS3m_vs', 'ScS5_vs', 'S6m_vs', 'ScS4m_vs', 'SS_vs', 'S4_vs', 'S5_vs', 'S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff', 'SSSm_mb', 'ScSScSScS_mb', 'ScSScSScSScS_mb', 'SSSSm_mb', 'SSSSSm_mb', 'ScSScSScSScSScS_mb', 'SSSSSSm_mb', 'ScSScSScSScSm_mb']
hipr_phases = ['S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff']


# INPUT DATA FILES #
dataset = 'eventinfo.clean.virtual_stack.all.csv' # 'eventinfo.multi_phase_single_pick.csv' #'eventinfo.comprehensive.6phase.Nov15.2019.csv'
raw_headers_to_keep = ['COMPREHENSIVE_WEIGHT']
data_wave_type = 'S' # either 'S' or 'P'

# REFERENCE MODEL #
reference_model = 'prem' # prem, (ak135, iasp91 eventually?)
total_radius = 6371. # km
core_mantle_boundary = 2891. # km


# 3D COORDINATE SYSTEM #
shell_bounds = [0., 24.4, 80., 160., 220., 310., 400., 490., 580., 670., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2891.] # km
discontinuities = [0., 15., 24.4, 220., 400., 600., 670., 771., 2741., 2891.] # km; should include discontinuities of the reference model
cardinal_azimuths = [0., 90., 180., 270., 360.] # degrees
make_near_neighbors_files = True
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


# INPUT MODEL PARAMETERES #
# crustal_model = 'CRUST_1.0' # name of crustal model to correct for
# crustal_model_increment = 1
# crustal_model_layers = ['water', 'ice', 'upper_sediments', 'middle_sediments', 'lower_sediments', 'upper_crust', 'middle_crust', 'lower_crust', 'mantle'] # ['ice', 'water', 'soft_sediments', 'hard_sediments', 'upper_crust', 'middle_crust', 'lower_crust', 'mantle'] # 
# crustal_correction_data_label = 'CRUST_1.0_CORR'

input_model = 'PREM' # name of input model file
residual_header = 'SYNTH_DT_S40RTS_dvs_1.0_5.0' #'CRUST_1.0_ELLIP_DT' #'DT' #'DT_CORR' # column name of the residual that you want to use in the data file.
wave_type = 'S' # 'S' or 'P' or 'SYNTH'
delimiter = '|' # delimiter in converted model.csv file
lat_header = 'latitude' # header in converted model.csv file
lon_header = 'longitude' # header in converted model.csv file
depth_header = 'depth' # header in converted model.csv file
property_header = 'vs' # header in converted model.csv file
update_reference = True # True if the input model needs to be re-referenced to the reference model used for the update, e.g., PREM
voigt = [False] #[False] [True, 'vpv', 'vph']


# OUTPUT FILE PARAMETERS #
# physical_property = 'vs'
# segment_property_header = f'SEG_V{data_wave_type}'
perturbation_header = 'dVs_%'
out_property_header = 'dvs (%)'
RMS = 'RMS_dVs'


# MODEL UPDATE LAYER STRIPPING PARAMETERS #
update_phases = [['S3_vs', 'ScS_vs', 'ScSScS_vs', 'Sdiff_vs', 'S_vs', 'SSm_vs', 'ScS4_vs', 'ScS3_vs', 'S4m_vs', 'S3m_vs', 'ScS5m_vs', 'S5m_vs', 'ScS3m_vs', 'ScS5_vs', 'S6m_vs', 'ScS4m_vs', 'SS_vs', 'S4_vs', 'S5_vs', 'S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff', 'SSSm_mb', 'ScSScSScS_mb', 'ScSScSScSScS_mb', 'SSSSm_mb', 'SSSSSm_mb', 'ScSScSScSScSScS_mb', 'SSSSSSm_mb', 'ScSScSScSScSm_mb']] # [['ScS_vs', 'ScSScS_vs', 'S_vs', 'SS_vs', 'S', 'SS', 'SSS', 'ScS', 'ScSScS']] #
# list of lists of phases to include in the model update. include name modifier. the number and order of lists must match the number and order of layers.
residual_limits = [[False]] #[[False]] or [[True, lower_lim, upper_lim]] argument for whether or not there are limits imposed on residual values. # list of lists should be the same length as the number of layers
layer_stripping_top_shells = [2] # list of the top-most shell(s) in a given layer(s)
layer_stripping_base_shells = [31] # list of the bottom-most shell(s) in a given layer(s)
freeze_previous_layers = [True] # [None, True, False]
iteration_to_stop_RMS_weighting = [2] # the iteration for each layer of `n` layers on which to stop weighting backmapped perturbations by RMS weighting. a list with `n` elements, either integers or 'None'.
variance_cutoff_type = ['total iterations'] # 'reduction' or 'total iterations'
variance_cutoff = [10] # if 'reduction', type == float; if 'total iterations', type == int
proportion_of_dataset_to_use = [1] #.5
multiprocessing_path_increments = [9000] # 9000 for total dataset, 7500 for limited dataset


# COVERAGE PARAMETERS #
azimuthal_sectors = 6


# SMOOTHING PARAMETERS #
smoothing_radii = [[3.0, 6.0, 9.0]] # list of lists of smoothing radii for each main layer in the current update
total_required_paths = [20] # this is the minimum number of paths that are required in the radius for smoothing
total_required_azimuths = [2] # this is the minimum number of azimuths that are required in each azimuthal sector for smoothing
peak_value = 1. # the value of the Gaussian at the center of the smoothing circle
peak_center = 0. # the center of the Gaussian fn on the x-axis. The radius at which values will get the highest weight
cutoff_weight = [0.5] # Value of the Gaussian at the smoothing cuttoff
azimuthal_weighting = [False] # True or False to include azimuthal weighting when computing the smoothed azimuthal pertubation mean.
hipr_weighting = [True] # True or False to include HIPR (JGR 2021) weights when computing the final weighted average perturbation.


############################
# NON-MODIFIABLE VARIABLES #
############################
# FILE STRUCTURE #
shell_file = 'shell_dimensions.csv'
block_file = 'block_dimensions.csv'
phases_directory = 'phases'
data_directory = 'phase_data'
# orig_paths_directory = 'original_path_files'
# resampled_paths_directory = 'resampled_path_files'
backmapped_paths_directory = 'backmapped_path_files'
# tomography_paths_directory = 'model_path_files'
paths_directory = 'raypath_files'
tomography_model_directory = 'models'
near_neighbors_directory = 'near_neighbors'
block_centric_directory = 'model_block_information'
main_headers = ['PHASE', 'EQ_LAT', 'EQ_LON', 'STA_LAT', 'STA_LON', 'EQ_DEP', 'DT']

# rounding values
computed_decimal_places = 10
rounded_decimal_places = 5





# dubious v
# start_depth = shell_bounds[0]
# final_depth = shell_bounds[-1]
# end_depth = shell_bounds[-2]
# depth_mins = shell_bounds[:-1]
# depth_maxs = shell_bounds[1:]
# shell_numbers = list(range(1, len(shell_bounds)))

# end_lat = final_lat - reference_lat  # second to last latitude value
# total_lat = final_lat - start_lat
# latitudes = list(range(start_lat, final_lat + reference_lat, reference_lat))

# end_lon = final_lon - reference_lon  # second to last longitude value
# total_lon = final_lon - start_lon

# lat_steps = int(total_lat / reference_lat)

azimuthal_sector_extent = 180. / azimuthal_sectors





