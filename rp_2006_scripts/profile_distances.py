#!/usr/bin/env python2.4

from __future__ import division

import psyco_full
import sys

from numpy import *
from random import *

from itertools import *
from math import *

from sets import Set

from cookbook.progress_bar import ProgressBar

# Word size
N = 5

# Treshold for a word to be considered over represented
T = .90

def main():
    # Read all the seqs
    seqs = [ tuple( map( int, line.split() ) ) for line in sys.stdin ]
    # Initialize progress bar
    bar = ProgressBar( 0, len( seqs ), 80 )
    # 'Matrix' of standardized values
    M = []
    # Open output
    out = open( sys.argv[1], "w" )
    # Iterate over all samples
    for count, ints in enumerate( seqs ):
        m = dict()
        # Standardize each N-mer that was observed
        R = 100
        # Observed 
        observed = kmer_counts( ints, N )
        allowed_keys = observed.keys()
        sample_counts = dict( ( key, zeros( R ) ) for key in allowed_keys )
        radix = max(ints) + 1
        l = len( ints )       
        tx = learn( ints, radix )
        for r in range( R ):
            seq = draw( tx, 0, l )
            for i in range( l - N ):
                key = seq[i:i+N]
                if key in allowed_keys:
                    sample_counts[key][r] += 1
        for key in allowed_keys:
            # Standardize by quantile style
            sample_counts[key].sort()
            q = int( len( sample_counts[key] ) * .75 )
            p75 = sample_counts[key][q]
            m[key] = observed[ key ] / p75
            # Threshold style
            #pi = sum( observed[key] > sample_counts[key] ) / R
            #if pi > T: m[key] = 1
            # Just standardize by length style
            #m[key] = observed[ key ] / l
            
        # Store the standardized values (sparsely)
        M.append( m )
        # Progress bar update
        bar.update_and_print( count )
    print
    # Compute 'distance'
    for i in range( len( seqs ) ):
        out.write( "V%d" % i )
        for j in range( len( seqs ) ):
            in_both = 0
            for key in M[i]:
                if key in M[j]:
                    in_both += ( M[i][key] * M[j][key] )
            # in_both = sum( [ key in M[j] for key in M[i] ] )
            out.write( "\t%.9f" % in_both )
        out.write( "\n" )    
        bar.update_and_print( i )

def learn( ints, radix ):
    tx = zeros( ( radix, radix ), Float )
    for i in range( 1, len( ints ) ):
        tx[ ints[i-1], ints[i] ] += 1
    # Convert each row to cumulative probabilities to allow sampling
    for row in tx:
        if sum( row ) > 0:
            row[:] = cumsum( row / sum( row ) )
    return tx

def draw( tx, prev, n ):
    rval = []
    for i in range( n ):
        prev = sample_transition( tx[prev] )
        rval.append( prev )
    return tuple( rval )

def sample_transition( tx_row ):
    r = random()
    for i, val in enumerate( tx_row ):
        if r < val: return i
    return len( tx_row ) - 1

def kmer_counts( ints, k  ):
    f = dict()
    for i in range( len(ints) - k ):
        key = ints[i:i+k]
        try: f[key] += 1
        except: f[key] = 1
    return f

if __name__ == "__main__": main()
