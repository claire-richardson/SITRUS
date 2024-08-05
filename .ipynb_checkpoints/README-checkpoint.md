# SITRUS (Seismic Iterative TomogRaphy Update Scheme)
SITRUS is an end-to-end software package designed to update tomography models by mapping predicted travel time residual misfits back into the starting model.


**Key things to do before you begin:**
- set up the SITRUS conda environment using `sitrus.yml`. from the shell, run `conda env create -f sitrus.yml`. this environment must be active (`conda activate sitrus`) to successfully execute SITRUS scripts.      

- read all `*readme.md` files. these will tell you the purpose of each script, the input parameters to set in the input file `mod_input.py`, and all of the input and output items.

- ensure your raw dataset is a CSV file formatted according to instructions in `_part1.readme.md`
  