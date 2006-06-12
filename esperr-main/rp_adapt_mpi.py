#!/usr/bin/env python

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'. Uses MPI via pypar. 

usage: %prog pos_data neg_data out_dir [options]
   -f, --format=NAME:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence
"""

import align.maf
import cookbook.doc_optparse
import os.path
import sys
import traceback

from cookbook.progress_bar import *

import rp.cv
import rp.io
import rp.mapping
import rp.models.averaging

# Startup pypar and get some info about what node we are
import pypar 
nodes = pypar.size() 
node_id = pypar.rank() 
print "I am node %d of %d" % ( node_id, nodes )

# Things that shouldn't be hardcoded'
stop_size = 5
fold = 10
passes = 10

model = rp.models.averaging

def run( pos_file, neg_file, out_dir, format, align_count, mapping ):

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    # Read integer sequences
    pos_strings = list( rp.io.get_reader( pos_file, format, None ) )
    neg_strings = list( rp.io.get_reader( neg_file, format, None ) )

    symbol_count = mapping.get_out_size()

    # Collapse
    while symbol_count > stop_size:

        # Sync nodes on each pass, may not be required
        pypar.barrier()

        if node_id == 0:
            print "Collapsing from:", symbol_count

        pairs = all_pairs( symbol_count )

        # Decide which subrange of all pairs this node will handle
        lo, hi = pypar.balance( len( pairs ), nodes, node_id )

        # Find best collapsed mapping in interval
        best_i, best_j, best_merit = None, None, 0
        for i, j in pairs[lo:hi]:
            merit = calc_merit( pos_strings, neg_strings, mapping.collapse( i, j )  ) 
            if merit > best_merit:
                best_i, best_j, best_merit = i, j, merit
            
        # Aggregate results
        if node_id != 0:
            # Send best i, j, merit to the master
            pypar.send( ( best_i, best_j, merit ), 0 )
        else:
            # I am the master, get results from all other nodes and determine
            # which had the best merit
            for other_node_id in range( 1, nodes ):
                i, j, merit = pypar.receive( other_node_id )
                if merit > best_merit:
                    best_i, best_j, best_merit = i, j, merit
  
        # Collapse the two symbols that resulted in the best merit   
        mapping = mapping.collapse( best_i, best_j )
        symbol_count -= 1
        
        # Ensure only the master writes files
        if node_id == 0:
        
            # Append merit to merit output        
            print >>merit_out, symbol_count, best_merit
        
            print "\nBest Merit %d." % best_merit,
        
            # Write best mapping to a file
            mapping_out = open( os.path.join( out_dir, "%03d.mapping" % symbol_count ), 'w' )
            for i, symbol in enumerate( mapping.get_table() ): 
                print >>mapping_out, str.join( '', rp.mapping.DNA.reverse_map( i, align_count ) ), symbol
            mapping_out.close()

def all_pairs( n ):
    rval = []
    for i in range( 0, n ):
        for j in range( i + 1, n ):
            rval.append( ( i, j ) )
    return rval   
    
def max_order( radix ):
    """Determine max order based on size of alphabet"""
    if radix <= 4: return 4
    elif radix <= 7: return 3
    elif radix <= 14: return 2
    else: return 1     

def calc_merit( pos_strings, neg_strings, mapping ):
    # Apply mapping to strings
    pos_strings = [ mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ mapping.translate( s ) for s in neg_strings ]
    # Cross validate using those strings
    radix = mapping.get_out_size()
    order = max_order( radix )
    model_factory = lambda d0, d1: model.train( order, radix, d0, d1 )
    cv_engine = rp.cv.CV( model_factory, pos_strings, neg_strings, fold=fold, passes=passes )
    cv_engine.run()
    # Merit is TP + TN
    return cv_engine.cls1.pos + cv_engine.cls2.neg

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    #try:
    pos_fname, neg_fname, out_dir = args
    align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
    #except:
    #    cookbook.doc_optparse.exit()

    try:
        run( open( pos_fname ), open( neg_fname ), out_dir, options.format, align_count, mapping )
    finally:
        pypar.finalize()


if __name__ == "__main__": main()
