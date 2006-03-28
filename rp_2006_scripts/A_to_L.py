#!/usr/bin/env python2.4

"""
"Standardize/scale an afinity matix"
"""

import sys
from numpy import *

def read_data( f, row_skip=0, col_skip=1 ):
    """Read a file containing column names, counts, and ancestral distributions"""
    input = iter( f )
    for i in range( row_skip ):
        input.readline()
    return array( [ map( float, line.split()[col_skip:] ) for line in input ] )

def save_data( D, fname ):
    f = open( fname, "w" )
    for row in D:
        for item in row:
            print >>f, item,
        print >> f
    f.close()

def main():
    # Read affinity / similarity / kernel matrix
    print "Reading"
    A = read_data( open( sys.argv[1] ) )
    N = A.shape[0]
    assert A.shape[1] == N, "Input affinity matrix must be square"
    # Make self similarities zero (per Ng, Jordan, Weiss)
    for i in range( N ):
        A[i,i] = 0
    # Diagonal matrix containing row sums of A
    D = identity( N ) * sum( A, 1 )
    D = D.astype( Float )
    # D = D ^ (-1/2)
    for i in range( len(D) ): 
        if D[i,i] > 0:
            D[i,i] = 1/sqrt(D[i,i])
    # L = D^(-1/2) * A * D^(-1/2)
    print "Scaling"
    L = dot( D, dot( A, D ) )
    # Save it
    save_data( L, sys.argv[2] )
    
    
if __name__ == "__main__": main()
