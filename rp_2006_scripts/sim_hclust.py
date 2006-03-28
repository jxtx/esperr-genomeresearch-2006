#!/usr/bin/env python2.4

import sys

from numpy import *

from cookbook.progress_bar import ProgressBar

INF = float( "inf" )

def read_data( f, row_skip=0, col_skip=0 ):
    """Read a file containing column names, counts, and ancestral distributions"""
    input = iter( f )
    for i in range( row_skip ):
        input.readline()
    return array( [ map( float, line.split()[col_skip:] ) for line in input ] )

def setup_nn( D ):
    N = D.shape[0]
    nearest = zeros( N )
    nearest_sim = zeros( N, Float )
    for i in range( N ):
        nearest[i] = argmax( D[i] )
        nearest_sim[i] = max( D[i] )
        assert nearest_sim[i] >= 0
    return nearest, nearest_sim
        
def update_nn( D, i, j, nearest, nearest_sim ):
    nearest[i] = argmax( D[i] )
    nearest_sim[i] = max( D[i] )
    for row in range( len(D) ):
        if nearest[row] == i or nearest[row] == j:
            nearest[row] = argmax( D[row] )
            nearest_sim[row] = max( D[row] )
    
def max_indexes( D, nearest, nearest_sim ):
    i = argmax( nearest_sim )
    j = nearest[i]
    assert i != j, "Nearest neighbor should not be self"
    return i, j
    
# def max_indexes( D, good_indexes ):
#     """Probably slow, there should be an easy Numeric way to do this"""
#     max_value = 0
#     max_i = None
#     max_j = None
#     good = nonzero( good_indexes )
#     for i in good:
#         if i == len(D) - 1: continue
#         t = where( good_indexes[i+1:], D[i,i+1:], -1 )
#         max_t = max( t )
#         if max_t > max_value:
#             max_value = max_t
#             max_i = i
#             max_j = argmax( t ) + (i+1)
#     return max_i, max_j

def main():
    # Read affinity / similarity / kernel matrix
    print >> sys.stderr, "Reading"
    L = read_data( sys.stdin )
    N = len( L )
    # Make self similarities small
    for i in range( N ): L[i,i] = -1
    
    nearest, nearest_sim = setup_nn( L )
    
    # ---- Cluster ----------------------------------------------------
        
    print >> sys.stderr, "Clustering"
    bar = ProgressBar( 0, N, 80 )
    
    
    n_clusters = N
    cluster_sizes = ones( N )
    good_indexes = ones( N )
    stage = 0
    cluster_stages = dict()
    while n_clusters > 1:
#       if n_clusters < 10:
#             print >> sys.stderr, "----------"
#             print >> sys.stderr, compress( good_indexes, compress( good_indexes, L, 1 ), 0 )
#             print >> sys.stderr, "----------"
#             print nonzero( good_indexes )
#             print "----------"
        # Find closest pair
        i, j = max_indexes( L, nearest, nearest_sim )
        similarity = L[i,j]
        print n_clusters - 1,
        # Print merge info
        if i in cluster_stages: 
            print cluster_stages[i] + 1,
        else: 
            print - (i+1),
        if j in cluster_stages: 
            print cluster_stages[j] + 1,
        else: 
            print - (j+1),
        print similarity,
        print
        # Update L
        L[i] = ( L[i] * cluster_sizes[i] + L[j] * cluster_sizes[j] ) / ( cluster_sizes[i] + cluster_sizes[j] )
        L[:,i] = L[i]
        # Keep self similarities small
        L[i,i] = -1
        cluster_sizes[i] += cluster_sizes[j]
        cluster_stages[i] = stage
        good_indexes[j] = 0   
        # Just to be safe
        L[j] = -1; L[:,j] = -1
        # Update NN
        nearest[j] == -1
        nearest_sim[j] == -1
        update_nn( L, i, j, nearest, nearest_sim )
        # Counters and status
        stage += 1
        n_clusters -= 1
        sys.stdout.flush()
        bar.update_and_print( N-n_clusters, sys.stderr )
    
    
if __name__ == "__main__": main()
