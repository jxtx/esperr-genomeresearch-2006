#!/usr/bin/env python2.4

from __future__ import division

import sys

labels = []

total_in = 0
num_in = 0

total_out = 0
num_out = 0

for line in open( sys.argv[1] ):
    labels.append( int( line ) )

for i, line in enumerate( sys.stdin ):
    fields = line.split()
    if fields[0].startswith( "V" ):
        fields = fields[1:]
    for j, val in enumerate( map( float, fields[(i+1):] ) ):
        if labels[i] == labels[j]:
            num_in += val
            total_in += 1
        else:
            num_out += val
            total_out += 1
    sys.stdout.flush()

print
print "average same class:       ", (num_in/total_in)
print "average different classes:", (num_out/total_out)
