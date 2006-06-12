#!/usr/bin/env python2.4

"""
Read 
    1) a mapping from columns to values
    2) a set of ancestral distributions for columns in that mapping
    3) a second set of ancestral distributions for columns (maybe) not in the
       first mapping
       
For each column not in the mapping the nearest neighbor that IS in the mapping
is identified and a new mapping is printed in which such columns are mapped to
the same symbol as said nearest neighbor.

"""

import sys
from math import *

from kdtree import *

inf = float( "inf" )

k = kdtree()
points = []
points_to_labels = dict()

# Read the mapping
mapping = {}
for line in open( sys.argv[1] ):
    fields = line.split()
    mapping[ fields[0] ] = fields[1]

# Read the ancestral distributions associated with that mapping
for line in open( sys.argv[2] ):
    if line.startswith( '#' ):
        continue
    fields = line.split()
    label = fields[0]
    point = tuple( map( float, fields[2:] ) )
    assert point not in points_to_labels, "Multiple columns in input map to the same point in space"
    # assert label in mapping, "Column in ancestral distributions but not in mapping"
    if label not in mapping:
        # Ignore points in the first adists but not in the mapping
        continue
    points_to_labels[point] = label
    k.addPoint( point )
    points.append( point )
   
# Read the second ancestral distributions and write a mapping
for line in open( sys.argv[3] ):
    if line.startswith( '#' ):
        continue
    fields = line.split()
    label = fields[0]
    point = tuple( map( float, fields[2:] ) )
    # <SLOW>
    ## best_dist = inf
    ## best_point = None
    ## for p in points:
    ##     d = dist( point, p )
    ##     if d < best_dist:
    ##         best_dist = d
    ##         best_point = p
    # </SLOW>
    nearest = k.nearestNeighbor( point, inf )[1] 
    ## assert nearest == best_point
    nearest = points_to_labels[ nearest ]
    print label, mapping[ nearest ]
