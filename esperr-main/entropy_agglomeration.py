#!/usr/bin/env python2.4

"""
Agglomeration procedure that maximizes entropy at each stage while considering
only agglomerations that join "neighboring" groups in some space.
"""

from __future__ import division
import sys
from numpy import *
# from rpy import *

def read_data( f ):
    """Read a file containing column names, counts, and ancestral distributions"""
    names = []
    counts = []
    rows = []    
    for line in f:
        if line.startswith( "#" ): continue
        fields = line.split()
        names.append( fields[0] )
        counts.append( int( fields[1] ) )
        rows.append( map( float, fields[2:] ) )
    return names, array( counts, dtype="i" ), array( rows, dtype="d" )
    
def euclidean_distance( X ):
    """
    Compute euclidean distance matrix between all rows of 'X'
    """
    n, m = X.shape
    D = zeros( (n,n), dtype="d" )
    for i in range( n ):
        D[i] = sqrt( sum( ( X - X[i] )**2, 1 ) )
        D[i,i] = float( "inf" )
    return D 
    
def update_euclidean_distance( X, D, i ):
    """
    Update distance matrix following a change affecting only X[i]
    """
    D[i] = sqrt( sum( ( X - X[i] )**2, 1 ) )
    D[:,i] = D[i] # Symmetry
    D[i,i] = float( "inf" )

def row_and_col_min_indexes( X ):
    """
    Return the indexes of the minimum values in each row and column
    """
    col_mins = argmin( X, 0 )
    row_mins = argmin( X, 1 )
    return row_mins, col_mins
           
def calc_entropy( p_c ):
    return - sum( p_c * log( p_c ) )
    #rval = 0
    #for p in p_c:
    #    rval -= p_c * log( p_c )
    #return rval
                
def update_entropy( h_prev, p_c, i, j ):
    return h_prev + p_c[i] * log( p_c[i] ) + p_c[j] * log( p_c[j] ) \
        - ( p_c[i] + p_c[j] ) * log( p_c[i] + p_c[j] ) 
      
def knn_agglomeration( X, counts, p_c, threshold, names ):
    """
    Group all entries in X that occur < threshold times with the nearest entry
    that occurs >= threshold times
    """
    X_to_clusters = arange( len( X ) )
    is_over_thresh = ( counts >= threshold )
    # The new cluster centroids
    new_X = compress( is_over_thresh, X, 0 )
    # Temporarily store the old positions for distance calculation (makes
    # this algorithm deterministic, otherwise depends on input order)
    temp_X = compress( is_over_thresh, X, 0 )
    new_counts = compress( is_over_thresh, counts, 0 )
    new_p_c = compress( is_over_thresh, p_c, 0 )
    new_index = 0
    for i, val in enumerate( is_over_thresh ):
        if val:
            X_to_clusters[i] = new_index
            new_index += 1
        else:
            dists = sqrt( sum( ( temp_X - X[i] )**2, 1 ) )
            closest = argmin( dists )
            sum_p = new_p_c[closest] + p_c[i]
            new_X[closest] = new_X[closest] * ( new_p_c[closest] / sum_p ) + X[i] * ( p_c[i] / sum_p )
            new_counts[closest] += counts[i]
            new_p_c[closest] += p_c[i]
            X_to_clusters[i] = closest
    return new_X, new_counts, new_p_c, X_to_clusters
    
def find_best_pair( X, D, counts, p_c, entropy, n_clusters ):
    """
    Find the pair of "contiguous" clusters which when combined result in
    the largest entropy
    """
    best_i = best_j = best_entropy = None
    n, m = X.shape
    row_min_indexes, col_min_indexes = row_and_col_min_indexes( D )
    
    # if plot:
    #     r.plot( X[:,0], counts, xlab="", ylab="", type="n" )
    #     if last_i is not None:
    #        r.abline( v=oldi_pos, lty=2 )
    #        r.abline( v=oldj_pos, lty=2 )
    #        r.abline( v=X[last_i,0] )
    
    for i in range( n ):
        j = row_min_indexes[i]
        if i < j: # and col_min_indexes[j] == i:
            
            # if plot:
            #     r.points( X[:,0], counts, pch=19, cex=0.5 )
            #     r.points( X[i,0], counts[i], cex=4 )
            #     r.points( X[j,0], counts[j], cex=4 )
            
            h = update_entropy( entropy, p_c, i, j )
            if h > best_entropy:
                best_i, best_j, best_entropy = i, j, h
            
            # print "Trying", i, j, "entropy:", h  / log(n_clusters)
        
            ##  Erase circles from plot
            # if plot:
            #     r.points( X[i,0], counts[i], cex=4, col="white" )
            #     r.points( X[j,0], counts[j], cex=4, col="white" )

    return best_i, best_j, best_entropy

                
def main():
    
    plot = False
   
    if len( sys.argv ) > 3:
        min_occurences = int( sys.argv[3] )
    else:
        min_occurences = 10
    
    # Read initial data points
    print "Reading data"
    names, counts, X = read_data( open( sys.argv[1] ) )
    print "Loaded %d data points" % len( X )
    
    # Number of clusters
    n_clusters = len( X )
    # Total number of observations
    total_mass = sum( counts )
    # Probability mass of each cluster
    p_c = counts / total_mass
    
    # If neccesary, do an initial knn agglomeration of infrequently
    # occuring points
    if min_occurences > 0 and not alltrue( counts >= min_occurences ):
        X, counts, p_c, X_to_clusters = knn_agglomeration( X, counts, p_c, min_occurences, names )
        n_clusters = len( X )
        print "Nearest neighbor agglomeration reduced data set to %d points" % len( X )
    else:
        X_to_clusters = arange( len( X ) )
    
    # Initial entropy
    entropy = - sum( p_c * log( p_c )  )
    
    # if plot:
    #     r.print_( "starting plot" ) 
    #     r.pdf( "movie.pdf" )
    
    # Distance matrix 
    print "Building initial distance matrix"
    D = euclidean_distance( X )
    
    out = open( sys.argv[2], "w" )
  
    last_i = None
    
    while n_clusters > 5:

        n, m = X.shape

        best_i, best_j, best_entropy = find_best_pair( X, D, counts, p_c, entropy, n_clusters )
        print "best entropy: ", best_entropy / log(n_clusters)
        
        i, j = best_i, best_j
        print "Merging", i, j, "clusters:", n_clusters-1
        
        entropy = update_entropy( entropy, p_c, i, j )
        
        # # Save for plot
        # if plot:
        #     oldi_pos, oldi_count = X[i,0], counts[i]
        #     oldj_pos, oldj_count = X[j,0], counts[j]
        
        # Merge the best clusters
        sum_p = p_c[i] + p_c[j]
        X[i] = X[i] * ( p_c[i] / sum_p ) + X[j] * ( p_c[j] / sum_p )
        p_c[i] += p_c[j]
        counts[i] += counts[j]
    
        # Compress
        X[j] = X[n_clusters-1]; X = X[:(n-1)]
        counts[j] = counts[n_clusters-1]; counts = counts[:(n-1)]
        p_c[j] = p_c[n_clusters-1]; p_c = p_c[:(n-1)]
        for row in range(len(X_to_clusters)): 
            if X_to_clusters[row] == j:
                X_to_clusters[row] = i 
            elif X_to_clusters[row] == (n-1):
                X_to_clusters[row] = j
    
        # Update D
        D = D[:(n-1),:(n-1)]
        update_euclidean_distance( X, D, i )
        if j < (n-1): update_euclidean_distance( X, D, j )
    
        # raw_input( "all done (whack a key)" )
        
        n_clusters -= 1 
        if n_clusters < 100:        
            print >>out, n_clusters, entropy / log(n_clusters), " ".join( map( str, X_to_clusters ) )
        last_i = i

    # if plot:
    #     r.dev_off()

if __name__ == "__main__": main()
