#!/usr/bin/env python2.4

"""
Convert one mapping from columns to atoms (AAA 4) and a second mapping from
atoms to a smaller set (4 2) and write a mapping from columns to the 
smaller set (AAA 2).
"""

import sys

atom_mapping = dict()
for line in open( sys.argv[2] ):
    fields = line.split()
    atom_mapping[ fields[0] ] = fields[1]
    
for line in open( sys.argv[1] ):
    fields = line.split()
    print fields[0], atom_mapping[ fields[1] ]