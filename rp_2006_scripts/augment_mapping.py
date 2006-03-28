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

# ---- KDTree taken from python cookbook (Eppstein) -------------------------

def dist2(p,q):
    """Squared distance between p and q."""
    d = 0
    for i in range(len(p)):
        d += (p[i]-q[i])**2
    return d

class kdtree:
    def __init__(self,dim=2,index=0):
        self.dim = dim
        self.index = index
        self.split = None

    def addPoint(self,p):
        """Include another point in the kD-tree."""
        if self.split is None:
            self.split = p
            self.left = kdtree(self.dim, (self.index + 1) % self.dim)
            self.right = kdtree(self.dim, (self.index + 1) % self.dim)
        elif self.split[self.index] < p[self.index]:
            self.left.addPoint(p)
        else:
            self.right.addPoint(p)

    def nearestNeighbor(self,q,maxdist2):
        """Find pair (d,p) where p is nearest neighbor and d is squared
        distance to p. Returned distance must be within maxdist2; if
        not, no point itself is returned.
        """
        solution = (maxdist2+1,None)
        if self.split is not None:
            solution = min(solution, (dist2(self.split,q),self.split))
            d2split = (self.split[self.index] - q[self.index])**2
            if self.split[self.index] < q[self.index]:
                solution = min(solution,
                    self.left.nearestNeighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.right.nearestNeighbor(q,solution[0]))
            else:
                solution = min(solution,
                    self.right.nearestNeighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.left.nearestNeighbor(q,solution[0]))
        return solution

# ---------------------------------------------------------------------------

inf = float( "inf" )

k = kdtree()
points = []
points_to_labels = dict()
wildcard_points = list()

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
    assert label in mapping, "Column in ancestral distributions but not in mapping"
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
    nearest = points_to_labels[ k.nearestNeighbor( point, inf )[1] ]
    print label, mapping[ nearest ]
