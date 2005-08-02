#!/usr/bin/env python

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'.

usage: %prog pos_data neg_data out_dir [options]
   -f, --format=NAME:  Format of input data. 'ints' by default, or 'maf'
   -a, --atoms=FILE:   A mapping specifying the largest set of symbols (these never get broken)
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence
   -M, --model=NAME:   Name of model to user
   -o, --order=D:      Order (max order?) or model
"""

from __future__ import division

import align.maf
import cookbook.doc_optparse
import os.path
import random
import sys
import time
import traceback

from Numeric import *
from cookbook.progress_bar import *
from itertools import *

import rp.cv
import rp.io
import rp.models
import rp.mapping


stop_size = 5
fold = 5
passes = 5

samp_size_collapse = 30
samp_size_expand = 10

def run( pos_file, neg_file, out_dir, format, align_count, atom_mapping, mapping, modname, modorder ):

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    # Read integer sequences
    print >>sys.stderr, "Loading training data"
    pos_strings = list( rp.io.get_reader( pos_file, format, None ) )
    neg_strings = list( rp.io.get_reader( neg_file, format, None ) )

    # Apply initial mapping immediately, to get the 'atoms' we will then collapse
    print >>sys.stderr, "Applying initial mapping"
    pos_strings = [ atom_mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ atom_mapping.translate( s ) for s in neg_strings ]

    # Count how many times each atom appears in the training data
    atom_counts = zeros( atom_mapping.get_out_size() )
    for string in chain( pos_strings, neg_strings ):
        for val in string:
            atom_counts[ val ] += 1

    # Valid candiates for expansion must occur more than 10 times in the training data
    can_expand = compress( atom_counts > 10, arange( len( atom_counts ) ) )

    # Handling bad columns in the training data is not obvious, so don't do it for now
    # for string in chain( pos_strings, neg_strings ):
    #    assert -1 not in string, "Cannot have invalid columns (map to -1) in training data"

    best_merit_overall = 0
    best_mapping_overall = None
    index_best_merit_overall = 0
    out_counter = 0

    step_counter = 0

    # 
    last_force_counter = 0

    print >>sys.stderr, "Searching"

    # Collapse
    while 1:

        symbol_count = mapping.get_out_size()


        best_merit = 0
        best_mapping = None

        clock = time.clock()
        cv_runs = 0

        # First try a bunch of collapses
        if symbol_count > stop_size:
            pairs = all_pairs( symbol_count )
            if len( pairs ) > samp_size_collapse: pairs = random.sample( pairs, samp_size_collapse ) 
            for i, j in pairs:
                new_mapping = mapping.collapse( i, j )
                merit = calc_merit( pos_strings, neg_strings, new_mapping, modname, modorder )
                cv_runs += 1
                if merit > best_merit:
                    best_merit = merit
                    best_mapping = new_mapping

        # Also try a bunch of expansions
        elements = random.sample( can_expand, samp_size_expand )
        for i in elements:
            new_mapping = mapping.expand( i )
            if new_mapping.get_out_size() == symbol_count: continue
            merit = calc_merit( pos_strings, neg_strings, new_mapping, modname, modorder )
            cv_runs += 1
            if merit > best_merit:
                best_merit = merit
                best_mapping = new_mapping

        clock = time.clock() - clock

        mapping = best_mapping

        # Append merit to merit output        
        print >>merit_out, step_counter, symbol_count, best_merit
        merit_out.flush()

        if best_merit >= best_merit_overall:
            best_merit_overall = best_merit
            best_mapping_overall = best_mapping
            # So we know what step the best mapping was encountered at
            best_merit_overall_index = step_counter
            restart_counter = step_counter
            # Reset the counter we use to force expansions
            last_force_counter = step_counter
            # Write best mapping to a file
            mapping_out = open( os.path.join( out_dir, "%03d.mapping" % out_counter ), 'w' )
            for i, symbol in enumerate( atom_mapping.get_table() ): 
                # Apply the 'second' mapping to the atom symbol
                if symbol >= 0: symbol = mapping[ symbol ]
                print >>mapping_out, str.join( '', rp.mapping.DNA.reverse_map( i, align_count ) ), symbol
            mapping_out.close()
            out_counter += 1

        print >>sys.stderr, "%06d, New best merit: %2.2f%%, size: %d, overall best: %2.2f%% at %06d, cvs per sec: %f" \
            % ( step_counter, best_merit * 100, mapping.get_out_size(), best_merit_overall * 100, best_merit_overall_index, clock/cv_runs  )

        # If we have gone 50 steps without improving over the best, restart from best
        if step_counter > restart_counter + 50:
            print >>sys.stderr, "Restarting from best mapping"
            print >>merit_out, step_counter, "RESTART"
            mapping = best_mapping_overall
            restart_counter = step_counter
            # Immediately force expansions after restart
            last_force_counter = 0

        if step_counter > last_force_counter + 20:
            last_force_counter = step_counter
            print >>sys.stderr, "Forcing expansions"
            print >>merit_out, step_counter, "FORCED EXPANSIONS"
            for i in range( 5 ):
                symbol_count = mapping.get_out_size()
                best_merit = 0
                best_mapping = None
                for i in random.sample( can_expand, samp_size_expand ):
                    new_mapping = mapping.expand( i )
                    if new_mapping.get_out_size() == symbol_count: continue
                    merit = calc_merit( pos_strings, neg_strings, new_mapping, modname, modorder )
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

def calc_merit( pos_strings, neg_strings, mapping, modname, modorder ):
    # Apply mapping to strings
    pos_strings = [ mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ mapping.translate( s ) for s in neg_strings ]
    # Cross validate using those strings
    radix = mapping.get_out_size()
    model_factory = lambda d0, d1: rp.models.train( modname, modorder, radix, d0, d1 )
    cv_engine = rp.cv.CV( model_factory, pos_strings, neg_strings, fold=fold, passes=passes )
    cv_engine.run()
    # Merit is TP + TN
    return ( cv_engine.cls1.pos / ( len( pos_strings ) * passes ) + cv_engine.cls2.neg / ( len( neg_strings ) * passes ) ) / 2

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    if 1:
        pos_fname, neg_fname, out_dir = args
        align_count, atom_mapping = rp.mapping.alignment_mapping_from_file( file( options.atoms ) )
        if options.mapping: mapping = rp.mapping.second_mapping_from_file( file( options.mapping ), atom_mapping )
        else: mapping = rp.mapping.identity_mapping( atom_mapping.get_out_size() )
        if options.model: modname = options.model
        else: modname = default_modname
        if options.order: modorder = int( options.order )
        else: modorder = default_modorder
        
    #except:
    #    cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), out_dir, options.format, align_count, atom_mapping, mapping, modname, modorder )


if __name__ == "__main__": main()
