#!/usr/bin/env python2.4

import sys
#import kdtree
from math import *
from Numeric import *

# def dist(p,q):
#     """Squared distance between p and q."""
#     d = 0
#     for i in range(len(p)):
#         d += (p[i]-q[i])**2
#     return sqrt( d )

inf = float( "inf" )

# k = kdtree.kdtree()
points = []
labels = []

wildcard_points = list()

mapping = {}
for line in open( sys.argv[1] ):
    fields = line.split()
    mapping[ fields[0] ] = fields[1]

for line in sys.stdin:
    if line.startswith( '#' ):
        continue
    fields = line.split()
    label = fields[0]
    point = tuple( map( float, fields[1:] ) )
    if '*' in label:
        wildcard_points.append( ( label, point ) )
    else:
        #assert point not in points_to_labels, "point already see! point: %r, old label: %r, new label: %r" % ( point, label, points_to_labels[point] )
        labels.append( label )
        points.append( point )

# Matrix of ceratin points
X = array( points )

# Now find nearest for each wildcard point

for label, point in wildcard_points:
#     best_dist = inf
#     best_point = None
#     for p in points:
#         d = dist( point, p )
#         if d < best_dist:
#             best_dist = d
#             best_point = p
#    nearest = points_to_labels[ k.nearestNeighbor( point, inf )[1] ]
    dists = sqrt( sum( ( X - array(point) )**2, 1 ) )
    closest = argmin( dists )
    print label, labels[ closest ], mapping[ labels[ closest ] ]

for key in sorted( mapping ):
    print key, mapping[key]
