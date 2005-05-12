#!/usr/bin/env python

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'.

usage: %prog pos_data neg_data out_dir [options]
   -f, --format=NAME:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence
"""

from __future__ import division

import align.maf
import cookbook.doc_optparse
import os.path
import random
import sys
import traceback

from cookbook.progress_bar import *
from rp import cv, io
from itertools import *

#import rp.models.averaging as model
#import rp.models.standard as model
#import rp.models.tree as model
import rp.models.tree_pruned_1 as model

import rp.mapping

stop_size = 5
fold = 5 
passes = 5 

samp_size_collapse = 30
samp_size_expand = 10

def run( pos_file, neg_file, out_dir, format, align_count, initial_mapping ):

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    # Read integer sequences
    pos_strings = list( io.get_reader( pos_file, format, None ) )
    neg_strings = list( io.get_reader( neg_file, format, None ) )

    # Apply initial mapping immediately, to get the 'atoms' we will then collapse
    pos_strings = [ initial_mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ initial_mapping.translate( s ) for s in neg_strings ]

    # Handling bad columns in the training data is not obvious, so don't do it for now
    for string in chain( pos_strings, neg_strings ):
        assert -1 not in string, "Cannot have invalid columns (map to -1) in training data"

    # And create the mapping we will collapse 
    mapping = rp.mapping.identity_mapping( initial_mapping.get_out_size() )

    best_merit_overall = 0
    best_mapping_overall = None
    index_best_merit_overall = 0
    out_counter = 0

    step_counter = 0

    last_force_counter = 0

    # Collapse
    while 1:

        symbol_count = mapping.get_out_size()

        best_merit = 0
        best_mapping = None

        # First try a bunch of collappses
        if symbol_count > stop_size:
            pairs = all_pairs( symbol_count )
            if len( pairs ) > samp_size_collapse: pairs = random.sample( pairs, samp_size_collapse ) 
            for i, j in pairs:
                new_mapping = mapping.collapse( i, j )
                merit = calc_merit( pos_strings, neg_strings, new_mapping )
                if merit > best_merit:
                    best_merit = merit
                    best_mapping = new_mapping

        # Also try expansions (some of them)
        elements = range( mapping.get_out_size() )
        if len( elements ) > samp_size_expand: elements = random.sample( elements, samp_size_expand )
        for i in elements:
            new_mapping = mapping.expand( i )
            if new_mapping.get_out_size() == symbol_count: continue
            merit = calc_merit( pos_strings, neg_strings, new_mapping )
            if merit > best_merit:
                best_merit = merit
                best_mapping = new_mapping

        mapping = best_mapping

        if best_merit >= best_merit_overall:
            best_merit_overall = best_merit
            best_mapping_overall = best_mapping
            best_merit_overall_index = step_counter
            last_force_counter = step_counter
            # Append merit to merit output        
            print >>merit_out, symbol_count, best_merit
            # Write best mapping to a file
            mapping_out = open( os.path.join( out_dir, "%03d.mapping" % out_counter ), 'w' )
            for i, symbol in enumerate( initial_mapping.get_table() ): 
                # Apply the 'second' mapping to the atom symbol
                if symbol >= 0: symbol = mapping[ symbol ]
                print >>mapping_out, str.join( '', rp.mapping.DNA.reverse_map( i, align_count ) ), symbol
            mapping_out.close()
            out_counter += 1

        print >>sys.stderr, "%06d, New best merit: %2.2f%%, size: %d, overall best: %2.2f%% at %06d" \
            % ( step_counter, best_merit * 100, mapping.get_out_size(), best_merit_overall * 100, best_merit_overall_index  )

        # If we have gone 100 steps without improving over the best, restart from best
        if step_counter > best_merit_overall_index + 50:
            print >>sys.stderr, "Restarting from best mapping"
            mapping = best_mapping_overall
            best_merit_overall_index = step_counter
            last_force_counter = step_counter

        if step_counter > last_force_counter + 20:
            last_force_counter = step_counter
            print >>sys.stderr, "Forcing expansions"
            for i in range( 5 ):
                symbol_count = mapping.get_out_size()
                best_merit = 0
                best_mapping = None
                for i in random.sample( range( mapping.get_in_size() ), samp_size_expand ):
                    new_mapping = mapping.expand( i )
                    if new_mapping.get_out_size() == symbol_count: continue
                    merit = calc_merit( pos_strings, neg_strings, new_mapping )
                    if merit > best_merit:
                        best_merit = merit
                        best_mapping = new_mapping
                mapping = best_mapping
        step_counter += 1

def all_pairs( n ):
    rval = []
    for i in range( 0, n ):
        for j in range( i + 1, n ):
            rval.append( ( i, j ) )
    return rval

def calc_merit( pos_strings, neg_strings, mapping ):
    # Apply mapping to strings
    pos_strings = [ mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ mapping.translate( s ) for s in neg_strings ]
    # Cross validate using those strings
    radix = mapping.get_out_size()
    #order = max_order( radix )
    order = 2
    model_factory = lambda d0, d1: model.train( order, radix, d0, d1, D=0.01, N=5 )
    cv_engine = cv.CV( model_factory, pos_strings, neg_strings, fold=fold, passes=passes )
    cv_engine.run()
    # Merit is TP + TN
    return ( cv_engine.cls1.pos / ( len( pos_strings ) * passes ) + cv_engine.cls2.neg / ( len( neg_strings ) * passes ) ) / 2

def max_order( radix ):
    """Determine max order based on size of alphabet"""
    if radix <= 4: return 4
    elif radix <= 7: return 3
    elif radix <= 14: return 2
    else: return 1

def not_max_order( radix ):
    return 3

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        pos_fname, neg_fname, out_dir = args
        align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
    except:
        cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), out_dir, options.format, align_count, mapping )


if __name__ == "__main__": main()
