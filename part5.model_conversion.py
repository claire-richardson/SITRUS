import sys
import os
import glob
import shutil
import getopt
import mod_input
from datetime import datetime, timezone
from netCDF4 import Dataset


def cd(path):
    os.chdir(os.path.expanduser(path))

models_name = mod_input.tomography_model_directory
input_csv = str(mod_input.input_model)+'.csv'
input_nc = str(mod_input.input_model)+'.nc'


if os.path.exists(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}') == True:
    print('path exists')
else:
    # make the .csv
    print('path does not exist')
    os.mkdir(f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/')
    cd(f'./{mod_input.tomography_model_directory}/')
    
    SCRIPT = os.path.basename(sys.argv[0])
    VERSION = 'V.2020.273'

    DEBUG = False

    GEOCSV_VERSION = 'GeoCSV2.0'

    # root _directory
    ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # optional directory to store netCDF and GeoCSV files
    DATA_DIR = os.path.join(ROOT_DIR, models_name)

    
    LON_VARIABLE = 'longitude'  # must match the netCDF file's longitude variable
    LAT_VARIABLE = 'latitude'  # must match the netCDF file's longitude variable
    DEPTH_VARIABLE = 'depth'  # must match the netCDF file's depth variable
    VALID_MODES = {'depth': 'km', 'single': ''}
    DELIMITER = '|'
    NETCDF_FILE_NAME = input_nc #None
    BASE_NAME = None
    VIEW_HEADER = False

    # output mode (depth | lat | lon) as individual files based on depth, lat, lon
    # or (single) as a single file
    OUTPUT_MODE = 'single'

    def get_variable_attributes(model_data, header, variable, variable_name):
        """add  variable attributes to the header

        Keyword arguments:
        model_data: Dataset instance of the model_file
        header: GeoCSV header variable for the model
        variable : variable to add
        variable_name: variable name used to represent this variable

        Return values:
        the geoCSV header for the model
        """
        header.append('# {}_column: {}\n'.format(variable, variable_name))
        for attr, value in vars(model_data.variables[variable]).items():
            if '_range' in attr:
                header.append('# {}_{}: {},{}\n'.format(variable, attr, value[0], value[1]))
            else:
                header.append('# {}_{}: {}\n'.format(variable, attr, value))
        return header


    def get_model_header(model_file, model_data):
        """create GeoCSV header for the model

        Keyword arguments:
        model_file: the netCDF model file name
        model_data: Dataset instance of the model_file

        Return values:
        the geoCSV header for the model
        """
        header = list()
        # GeoCSV header
        header.append('# dataset: {}\n'.format(GEOCSV_VERSION))
        header.append('# created: {} UTC ({})\n'.format(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), SCRIPT))
        header.append('# netCDF_file: {}\n'.format(os.path.basename(model_file)))
        header.append('# delimiter: {}\n'.format(DELIMITER))

        # global attributes
        history_done = False
        history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")} Converted to GeoCSV by {SCRIPT} ,' \
            f'{VERSION} ' \
            f'from {NETCDF_FILE_NAME}'
        for attr, value in vars(model_data).items():
            if isinstance(value, str):
                value = value.replace('\n', '; ')
            if attr.lower() == 'history':
                value = f'{history}; {value}'
                history_done = True
            header.append(f'# global_{attr}: {value}\n')

        if not history_done:
            header.append(f'# global_history: {history}\n')

        # variable s
        header = get_variable_attributes(model_data, header, LAT_VARIABLE, 'latitude')#, LAT_VARIABLE)
        header = get_variable_attributes(model_data, header, LON_VARIABLE, 'longitude')#, LON_VARIABLE)
        header = get_variable_attributes(model_data, header, DEPTH_VARIABLE, 'depth')#, DEPTH_VARIABLE)

        return ''.join(header)


    def get_var_header(model_data, var):
        """create GeoCSV header for a variable

        Keyword arguments:
        model_data: Dataset instance of the model_file
        var: the netCDF model file variable

        Return values:
        the geoCSV header for the variable
        """
        header = list()
        header = get_variable_attributes(model_data, header, var, var)
        return ''.join(header)

    def check_netcdf_file():
        """check the input netCDF model file and make sure it exist and extract info

        Return values:
        this_file:  file name with no extension
        """
        # check the model file and extract necessary information
        # must be in the argument list
        if NETCDF_FILE_NAME is None:
            print('[ERR] the netCDF model file name is required', flush=True)
            # usage_csv()
            sys.exit(1)

        # user may provide full path
        elif os.path.isfile(NETCDF_FILE_NAME):
            model_file_name = NETCDF_FILE_NAME
            base_file_name, ext = os.path.splitext(model_file_name)

        # user may place it under the data directory
        elif os.path.isfile(os.path.join(DATA_DIR, NETCDF_FILE_NAME)):
            model_file_name = os.path.join(DATA_DIR, NETCDF_FILE_NAME)
            base_file_name, ext = os.path.splitext(model_file_name)

        # could not find the file
        else:
            print('[ERR] could not find the netCDF model file {}'.format(NETCDF_FILE_NAME), flush=True)
            # usage_csv()
            sys.exit(1)

        return model_file_name, base_file_name

    def make_model_geocsv():
        """create GeoCSV file from a netCDF model file

        Keyword arguments:
        model_file: the netCDF model file name
        """
        model_file, base_file_name = check_netcdf_file()

        data_header = list()
        model_data = Dataset(model_file)

        try:
            # conversion to string is done to preserve precision
            lat = list()
            lon = list()
            depth = list()
            for this_value in model_data.variables[LAT_VARIABLE][:]:
                lat.append("{}".format(str(this_value)))
            for this_value in model_data.variables[LON_VARIABLE][:]:
                lon.append("{}".format(str(this_value)))
            for this_value in model_data.variables[DEPTH_VARIABLE][:]:
                depth.append("{}".format(str(this_value)))
        except Exception as ex:
            print('\n[ERR] the expected variables ({}, {}, {}) not in the variable list: {}\n'.format(
                LAT_VARIABLE, LON_VARIABLE, DEPTH_VARIABLE, str(list(model_data.variables.keys()))))
            sys.exit(1)

        emcin = {}

        # make sure this is a 3D netCDF file
        var_3d = list()
        for var in model_data.variables.keys():
            if len(model_data.variables[var].shape) == 3:
                var_3d.append(var)
        if len(var_3d) <= 0:
            print('\n[ERR] not a 3D netCDF file\n\n', flush=True)
            sys.exit(1)

        output_data = list()
        for k, this_depth in enumerate(depth):
            if OUTPUT_MODE == 'single' and k == 0:
                data_header = list()
                output_file = '{}.csv'.format(base_file_name)
                fp = open(output_file, 'w')
                # print('[INFO] Output file: {}'.format(output_file), flush=True)
                fp.write(get_model_header(model_file, model_data))
                data_header.append('{}{}{}{}{}'.format(LAT_VARIABLE, DELIMITER, LON_VARIABLE,
                                                       DELIMITER, DEPTH_VARIABLE))

            if DEBUG:
                print('[INFO] Processing depth: {}'.format(this_depth), flush=True)

            index = [-1, -1, -1]
            for i, this_lat in enumerate(lat):
                for j, this_lon in enumerate(lon):
                    if OUTPUT_MODE == 'single':
                        output_data.append('{}{}{}{}{}'.format(str(this_lat), DELIMITER, str(this_lon), DELIMITER,
                                                               str(this_depth)))
                    for var in model_data.variables.keys():
                        depth_index = None
                        lat_index = None
                        lon_index = None
                        if var.encode('ascii', 'ignore').decode("utf-8") not in [LAT_VARIABLE, LON_VARIABLE,
                                                                                 DEPTH_VARIABLE]:
                            if ((OUTPUT_MODE == 'single' and (not i and not j and not k)) or
                                    (OUTPUT_MODE == 'depth' and (not i and not j))):
                                fp.write(get_var_header(model_data, var))
                                data_header.append('{}{}'.format(DELIMITER, var))
                            # find the variable ordering
                            if lat_index is None:
                                for l in range(len(model_data.variables[var].dimensions)):
                                    if model_data.variables[var].dimensions[l].encode('ascii', 'ignore').decode(
                                           "utf-8") == DEPTH_VARIABLE:
                                        depth_index = l
                                    elif model_data.variables[var].dimensions[l].encode('ascii', 'ignore').decode(
                                            "utf-8") == LON_VARIABLE:
                                        lon_index = l
                                    else:
                                        lat_index = l

                                if var not in emcin.keys():
                                    try:
                                        emcin[var] = model_data.variables[var][:]
                                    except Exception as err:
                                        print('\n[ERR] problem reading variable "{}"'.format(var))
                                        print('{0}\n'.format(err))
                                        sys.exit(2)

                            index[depth_index] = k
                            index[lat_index] = i
                            index[lon_index] = j
                            # nan values, otherwise we write string to preserve the precision
                            if str(emcin[var][index[0]][index[1]][index[2]]) == '--':
                                output_data.append('{}{}'.format(DELIMITER,
                                                                 float(emcin[var][index[0]][index[1]][index[2]])))
                            else:
                                # conversion to string is done to preserve precision
                                output_data.append('{}{}'.format(DELIMITER,
                                                                 str(emcin[var][index[0]][index[1]][index[2]])))
                    output_data.append('\n')

        if OUTPUT_MODE == 'single':
            fp.write('{}\n'.format(''.join(data_header)))
            fp.write(''.join(output_data))
            fp.close()
    make_model_geocsv()

    cd('../')
    shutil.move(f'./{mod_input.tomography_model_directory}/{input_csv}', f'./{mod_input.tomography_model_directory}/{mod_input.data_wave_type}/{mod_input.input_model}_update/{input_csv}')
    print('done')
