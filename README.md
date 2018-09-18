# loop_helix_loop_reshaping
Reshape a part of the protein into a loop-helix-loop unit

## Installation
I recommand to create a python virtual environment before installing the package.
```
virtualenv -p /usr/bin/python3 --system-site-packages venv
```
Activate the virtual environment and install the package by
```
pip install -e .
```

## Running an example
This is an example that remodels two loop-helix-loop units in the `test_inputs/2lv8_cleaned.pdb` structure. The insertion points for the LHL units are defined in the file `test_inputs/2lv8_insertion_points.json`. 

### Select linker loops
The first step is selecting linker loops that do not clash with the scaffold.
```
./run_jobs.py test_select_linkers job_scripts/select_linkers_example.py
```
This command will create the directory `data/test_select_linkers` and dump 16 json libraries for selected linker loops at different positions.

### Screen loop helix loop units
The second step is screening LHL units that at each insertion points individually. At each insertion point, pairs of linker loops will be inserted and a helix will be introduced to bridge the gap. The built LHL unit will be filtered by the quality of the gap closure, the helix-scaffold contact and the number of buried unsatisfied H-bonds.
```
./run_jobs.py test_screen_single_insertion_loop_helix_loop_units job_scripts/screen_single_insertion_loop_helix_loop_units_example.py
```
The command will create the directory `data/test_screen_single_insertion_loop_helix_loop_units` and dump 2 json libraries for LHL units at the two insertion points.
