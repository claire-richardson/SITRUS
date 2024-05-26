README file for: `.py`


PURPOSE OF SCRIPT:
`pt3.src.make_plot_ready_files.py` makes files that can be read and plotted by the matplotlib utility Cartopy. This is achieved by further interpolating the grid registered files that are outputs of pt2 to a standardized user-defined grid.


OUTPUT/SAVED FILES:
    1. interpolated plot-ready files for each shell in `shell_dimensions.csv`: `./models/{wave_type}_models/{model_name}_{wave_type}/original_model_plot_files/{model_name}_shell_{shell_no}_*.csv`


BEFORE RUNNING:
    1. define necessary variables in `mod_input.py`.


INPUT VARIABLES TO DEFINE IN `mod_input.py`:
    1. `input_model`: the name of the model to be re-interpolated. this should exactly match the name of original converted CSV file (e.g., `SPiRaL` for `SPiRaL.csv`).
    2. `wave_type`: the wave type(s) (S and/or P) of the model to be re-interpolated.
    3. `reference_lat`: the spacing between latitude lines for the re-interpolated model, which will ultimately be the plotted model. value must be a factor of 180.
    4. `reference_lon`: the spacing between longitude lines for the re-interpolated model, which will ultimately be the plotted model. value must be a factor of 360.
    