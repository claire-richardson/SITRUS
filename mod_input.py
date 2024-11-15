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
user_email = 'crricha5@asu.edu' # full email address for process start/end notifications

# PHASE AND DATA NAMING #
all_phases = ['S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff'] # list of phases in the travel time dataset (example phases from Lai et al., 2019 and Lai & Garnero 2020)

# INPUT DATA FILES #
dataset = 'traveltimes.csv' # travel time dataset filename
raw_headers_to_keep = ['COMPREHENSIVE_WEIGHT', 'HIPR_BIN', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT'] # list of headers in the dataset file to keep. Comprehensive weight is from Lai et al., 2019 and Lai & Garnero 2020; HIPR weights are from Hansen et al., 2021.
data_wave_type = 'S' # either 'S' or 'P', wave type of the model to be updated and of the measurements in the dataset. S refers to SH-polarity S-waves


# REFERENCE MODEL #
reference_model = 'prem' # currently prem is the only supported model (Dziewonski & Anderson, 1981)
total_radius = 6371. # km; total radius of the reference model
core_mantle_boundary = 2891. # km; core-mantle boundary depth of the reference model


# 3D MESH DEFINITION #
shell_bounds = [0., 24.4, 80., 160., 220., 310., 400., 490., 580., 670., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2891.] # km
discontinuities = [0., 15., 24.4, 220., 400., 600., 670., 771., 2741., 2891.] # km; should include discontinuities of the reference model. these should also be included in `shell_bounds`.
make_near_neighbors_files = True
max_near_neighbors_radius = 15. # degrees; radius at which to calculate all near neighbors; default: 15
reference_lat = 2 # degrees; side length of 2D approximate equal area reference block for lateral component of the model space
reference_lon = 2 # degrees; side length of 2D approximate equal area reference block for lateral component of the model space
start_lat = -90  # degrees; starting latitude for your coordinate system (currently only global models are supported)
final_lat = 90  # degrees; ending latitude for your coordinate system (currently only global models are supported)
start_lon = -180  # degrees; starting longitude for your coordinate system (currently only global models are supported)
final_lon = 180 # degrees; ending longitude for your coordinate system (currently only global models are supported)


# PATH RESAMPLING PARAMETERS #
target_path_length = 80 # km
target_path_length_tolerance = 5 # km


# COVERAGE PARAMETERS #
azimuthal_sectors = 6


# MODEL PREPROCESSING
convert_nc_to_csv = False
zero_shells = [1]


## INPUT MODEL PARAMETERES #
input_model = 'S40RTS' # name of input model file
delimiter = ',' # delimiter in converted model.csv file
lat_header = 'latitude' # header in converted model.csv file
lon_header = 'longitude' # header in converted model.csv file
depth_header = 'depth' # header in converted model.csv file
property_header = 'vsh' # header in converted model.csv file
voigt = [False] #[False] [True, 'vpv', 'vph'] 
convert_vel_to_perturb = True


## MODEL UPDATE PARAMETERS
dataset_description = ['sitrus model update'] # DONT USE COMMAS!, short description of the model update
update_phases = ['S', 'SS', 'SSS', 'ScS', 'ScSScS', 'Sdiff'] # the same set or subset of `all_phases` to be used in the model update
type_of_phase_subselection = ['proportion', 'proportion'] #['proportion'] # 'proportion' or 'number'
subselection_of_phase_data_to_use = [1, 1, 1, 1, 1, 1]
residual_header = ['CRUST_ELLIP_DT'] # column name of the residual that you want to use in the data file. set to CRUST_ELLIP_DT to use the SITRUS corrected residuals with CRUST1.0 (Laske et al., 2013) and EllipticiPy (Russell et al., 2022).
residual_limits = [False] #[[False]] or [[True, lower_lim, upper_lim]] argument for whether or not there are limits imposed on residual values. # list of lists should be the same length as the number of layers
layer_top_shells = [2] # list of the top-most shell(s) in a given layer(s)
layer_base_shells = [31] # list of the bottom-most shell(s) in a given layer(s)
freeze_previous_layers = [True] # None, True, or False; indicate whether to freeze the previous update layer
remove_residual_mean = [True] # indicate whether or not to remove the residual mean from each event as a proxy for event relocation
starting_RMS_model_to_use = input_model # default is to use the same RMS profile as the starting model
iteration_to_stop_RMS_weighting = [1] # the iteration for each layer of `n` layers on which to stop weighting backmapped perturbations by RMS weighting. a list with `n` elements, either integers or 'None'.
cutoff_type = ['total iterations'] # 'reduction' or 'total iterations'
cutoff = [5] # if 'reduction', type == float; if 'total iterations', type == int
HPC_cores = 31 # should equal the number of depth shells (one less than the length of the `depth_bounds` list)


# SMOOTHING PARAMETERS #
smoothing_radii = [3., 4., 5.] # list of lists of smoothing radii for each main layer in the current update
total_required_paths = [20] # this is the minimum number of paths that are required in the radius for smoothing
total_required_azimuths = [2] # this is the minimum number of azimuths that are required in each azimuthal sector for smoothing
gaussian_cutoff_weight = [0.5] # Value of the Gaussian at the smoothing cuttoff
azimuthal_weighting = [True] # True or False to include azimuthal weighting when computing the smoothed azimuthal pertubation mean.
special_weights = ['COMPREHENSIVE_WEIGHT', 'HIPR_STA_WEIGHT', 'HIPR_EQ_WEIGHT']


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
main_headers = ['EQ_DATE', 'PHASE', 'EQ_LAT', 'EQ_LON', 'STA_LAT', 'STA_LON', 'EQ_DEP', 'DT']
cardinal_azimuths = [0., 90., 180., 270., 360.] # degrees

# rounding values
computed_decimal_places = 10
rounded_decimal_places = 5




