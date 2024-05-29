## Purpose of script:
`part1.grid_definition.py` makes the global 3D block mesh to be used for SITRUS model updates. It also computes files for a static `near_neighbors` directory to be referenced in other parts of the software.


## Output/saved files:
1. `shell_dimensions.csv`: a file containing the dimensions of the depth shells as defined in `mod_input.py`. The depth shells are constant thickness, and each shell is given a Shell ID.
2. `block_dimensions.csv`: a file containing the dimensions of the 2D spherical squares on a given depth surface as defined in `mod_input.py`. SITRUS computes approximate equal area spherical squares for an arbitrary depth surface based on target latitude/longitude extents. Each block is given a Block ID. To make the 3D mesh, the same set of 2D blocks in `block_dimensions.csv` are projected radially between each depth surface in `shell_dimensions.csv`. Each unique model element/3D block in the mesh therefore requires both a Shell ID and a Block ID to be identified. Using this convention, all model elements in each depth shell share latitude/longitude boundaries.
3. `/near_neighbors` directory: a directory to hold all near neighbors files for each 2D block in `block_dimensions.csv`. This part of the script is *optional*, but should be run anytime a new mesh is generated. It takes much more time to run than the block and shell dimension computations.
4. `/near_neighbors/block_*_neighbors.csv`: CSV files for each 2D block in `block_dimensions.csv`. Each file contains:
    - the Block IDs of each of that block's near neighbors (`NEIGHBOR`)
    - the epicentral distance between the center of the primary block and each of its neighbors in degrees (`RADIUS_DEG`)


## Before running:
1. define necessary variables in `mod_input.py`.


## Input variables to define in `mod_input.py`:
1. `shell_bounds`: list of boundaries that define constant thickness depth shells. This should include the major discontinuities of the reference model (default: PREM). It should include all boundaries, including topmost and bottommost values. **the first shell must be consistent with the crust of the reference model.**
2. `total_radius`: total radius of the Earth according to the reference model (default: PREM)
3. `reference_lat`: the latitudinal extent of a reference 2D block at the equator of a depth surface from which to calculate the approximate equal area of all 2D blocks in `block_dimensions.csv` (default: 2; _must be an integer factor of the total latitude of the model space_)
4. `reference_lon`: the longitudinal extent of a reference 2D block at the equator of a depth surface from which to calculate the approximate equal area of all 2D blocks in `block_dimensions.csv` (default: 2; _must be an integer factor of the total longitude of the model space_)
5. `start_lat`: the latitude from which to start defining the 2D block dimensions (default: -90)
6. `final_lat`: the latitude at which to stop defining the 2D block dimensions (default: 90)
7. `start_lon`: the longitude from which to start defining the 2D block dimensions (default: -180)
8. `final_lon`: the latitude at which to stop defining the 2D block dimensions (default: 180)
9. `make_near_neighbors_files`: boolean indicating whether to make near neighbors files (default: True; **warning**: this will increase the walltime to ~several hours)
10. `max_near_neighbors_radius`: the total radius at which a neighboring block is considered a "near neighbor" (default: 15). Choose this number carefully, as it will be the maximum radius to which you can smooth during the adaptive smoothing process in a model update.


## Necessary compute resources:
1. HPC cores: 1
2. Time allocation: With near neighbors-generally a few hours, depending on block size and total number. Without near neighbors: >1 minute.