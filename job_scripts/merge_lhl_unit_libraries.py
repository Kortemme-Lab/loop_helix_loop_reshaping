#!/usr/bin/env python3
'''Merge the loop helix loop unit libraries. Libraries for the same
insertion point will be merged together.
Usage:
    ./merge_lhl_unit_libraries.py data_dir
'''

import os
import sys

import json


if __name__ == '__main__':
    data_dir = sys.argv[1]
    
    # Load the lhl_units 
    
    lhl_units = {}

    for lf in os.listdir(data_dir):
        if lf.startswith('selected_lhl_units_'):
            i = int(lf.split('_')[3])

            if not (i in lhl_units.keys()):
                lhl_units[i] = []

            with open(os.path.join(data_dir, lf), 'r') as f:
                lhl_units[i] += json.load(f)
   
    # Remove the existing files
    
    for lf in os.listdir(data_dir):
        if lf.startswith('selected_lhl_units'):
            os.remove(os.path.join(data_dir, lf))

    # Dump the merged libraries

    for i in lhl_units.keys():
        print('There are {0} LHL units for insertion point {1}'.format(len(lhl_units[i]), i))
        
        with open(os.path.join(data_dir, 'selected_lhl_units_{0}_0.json'.format(i)), 'w') as f:
            json.dump(lhl_units[i], f)
