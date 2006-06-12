#!/usr/bin/env python2.4

from __future__ import division

import sys

labels = []
for line in open( sys.argv[1] ):
    labels.append( line.strip() )

for i, line in enumerate( open( sys.argv[2] ) ):
    fields = line.split()
    if fields[0].startswith( "V" ):
        fields = fields[1:]
    for j, val in enumerate( map( float, fields[(i+1):] ) ):
        if labels[i] == labels[j]:
            group = labels[i]
            print "IN_" + group, val
        else:
            print "OUT_" + labels[i], val
            print "OUT_" + labels[j], val
            print "OUT", val

