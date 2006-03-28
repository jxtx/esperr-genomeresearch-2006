#!/usr/bin/env python

import sys

adists = open( sys.argv[1] )
clusters = open( sys.argv[2] )

desired_k = int( sys.argv[3] )

for line in clusters:
    fields = line.split()
    if int( fields[0] ) == desired_k:
        cl = fields[2:]
        break

cols = []
for line in adists: 
    if line.startswith( "#" ):
        continue
    cols.append( line.split()[0] )

assert len( cl ) == len( cols )

for a, b in zip( cols, cl ):
    print a, b

