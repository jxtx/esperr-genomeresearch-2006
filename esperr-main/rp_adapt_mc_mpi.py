#!/usr/bin/env python2.4

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'.

usage: %prog training_set_fnames... [options]
   -d, --out=DIR:       Output directory
   -f, --format=NAME:   Format of input data. 'ints' by default, or 'maf'
   -a, --atoms=FILE:    A mapping specifying the largest set of symbols (these never get broken)
   -m, --mapping=FILE:  A mapping (alphabet reduction) to apply to each sequence
   -M, --model=NAME:    Name of model to user
   -o, --order=D:       Order (max order?) or model
   -F, --fold=N:        Fold for cross validation
   -P, --passes=N:      Passes for cross validation
   -l, --loo:           Use leave-one-out cross validation rather than folds
   -p, --mpi:           Expect to be run under mpi (requires pypar!)
"""

from __future__ import division

import pkg_resources
pkg_resources.require( "bx-python" )

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

mpi = False
pypar = None
node_id = -1
nodes = 0

def message( *args ):
    """
    Write a message to stderr (but only on the master node if we are using pypar)
    """
    global mpi, node_id
    if not mpi or node_id == 0:
        sys.stderr.write( ' '.join( map( str, args ) ) )
        sys.stderr.write( '\n' )
        sys.stderr.flush()

stop_size = 5

fold = 5
passes = 5
loo = False

min_cols = 50

def run( ts_fnames, out_dir, format, align_count, atom_mapping, mapping, modname, modorder  ):

    samp_size_collapse = 30
    samp_size_expand = 10
    
    if mpi:
        global pypar, node_id, nodes
        # Startup pypar and get some info about what node we are
        pypar = __import__( 'pypar' )
        nodes = pypar.size() 
        node_id = pypar.rank() 
        print "I am node %d of %d" % ( node_id, nodes )
        # Modify these, they get split over nodes
        samp_size_collapse = samp_size_collapse // nodes
        samp_size_expand = samp_size_expand // nodes

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    # Read integer sequences
    message( "Loading training data" )
    
    training_sets = []
    for fname in ts_fnames:
        strings = []
        skipped = 0
        for s in rp.io.get_reader( open( fname ), format, None ):
            # Apply initial mapping
            s = atom_mapping.translate( s )
            # Ensure required columns
            if sum( s != -1 ) < min_cols:
                skipped += 1
                continue
            # Add to set
            strings.append( s )
        # Informational
        message( "Loaded training data from '%s', found %d usable strings and %d skipped" \
            % ( fname, len( strings ), skipped ) )
        training_sets.append( strings )

    # Count how many times each atom appears in the training data, valid 
    # candiates for expansion must occur more than 10 times in the training 
    # data.
    message( "Finding expandable atoms" )
    atom_counts = zeros( atom_mapping.get_out_size() )
    for string in chain( * training_sets ):
        for val in string:
            atom_counts[ val ] += 1
    can_expand = compress( atom_counts > 10, arange( len( atom_counts ) ) )

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    best_merit_overall = 0
    best_mapping_overall = None
    index_best_merit_overall = 0
    out_counter = 0

    step_counter = 0
    last_force_counter = 0
    
    message(  "Searching" )

    # Collapse
    while 1:

        clock = time.clock()
        cv_runs = 0

        if mpi:
            # Sync up nodes at start of each pass
            pypar.barrier()

        symbol_count = mapping.get_out_size()

        best_i = None
        best_j = None
        best_merit = 0
        best_mapping = None

        # First try a bunch of collapses
        if symbol_count > stop_size:
            # Select some random pairs from the region owned by this node
            pairs = all_pairs( symbol_count )
            if mpi:
                lo, hi = pypar.balance( len( pairs ), nodes, node_id )
                pairs = pairs[lo:hi]
            if len( pairs ) > samp_size_collapse: 
                pairs = random.sample( pairs, samp_size_collapse )
            # Try collapsing each pair 
            for i, j in pairs:
                new_mapping = mapping.collapse( i, j )
                merit = calc_merit( training_sets, new_mapping, modname, modorder  )
                cv_runs += 1
                if merit > best_merit:
                    best_i, best_j = i, j
                    best_merit = merit
                    best_mapping = new_mapping

        # Also try a bunch of expansions
        if mpi:
            lo, hi = pypar.balance( len( can_expand ), nodes, node_id )
            elements = random.sample( can_expand[lo:hi], samp_size_expand )
        else:
            elements = random.sample( can_expand, samp_size_expand )
        for i in elements:
            new_mapping = mapping.expand( i )
            if new_mapping.get_out_size() == symbol_count: 
                continue
            merit = calc_merit( training_sets, new_mapping, modname, modorder  )
            cv_runs += 1
            if merit > best_merit:
                best_i, best_j = i, None
                best_merit = merit
                best_mapping = new_mapping

        clock = time.clock() - clock

        if mpi:
            best_i, best_j, best_merit, cv_runs = sync_nodes( best_i, best_j, best_merit, cv_runs )  
            # Collapse or expand (if j is None) to get the overall best mapping
            if best_j is None:
                best_mapping = mapping.expand( best_i )
            else:
                best_mapping = mapping.collapse( best_i, best_j )

        mapping = best_mapping

        # Append merit to merit output
        if not mpi or node_id == 0:
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
            if not mpi or node_id == 0:
                write_mapping( mapping, os.path.join( out_dir, "%03d.mapping" % out_counter ) )
            out_counter += 1

        message( "%06d, New best merit: %2.2f%%, size: %d, overall best: %2.2f%% at %06d, cvs/sec: %f" \
                  % ( step_counter, best_merit * 100, mapping.get_out_size(), best_merit_overall * 100, best_merit_overall_index, cv_runs/clock  ) )

        # If we have gone 50 steps without improving over the best, restart from best
        if step_counter > restart_counter + 50:
            message( "Restarting from best mapping" )
            if not mpi or node_id == 0:
                print >>merit_out, step_counter, "RESTART"
            mapping = best_mapping_overall
            restart_counter = step_counter
            # Immediately force expansions after restart
            last_force_counter = 0

        if step_counter > last_force_counter + 20:
            last_force_counter = step_counter
            message( "Forcing expansions" )
            if not mpi or node_id == 0:
                print >>merit_out, step_counter, "FORCED EXPANSIONS"
            if mpi:
                lo, hi = pypar.balance( len( can_expand ), nodes, node_id )
                my_can_expand = can_expand[lo:hi]
            else:
                my_can_expand = can_expand
            for i in range( 5 ):
                symbol_count = mapping.get_out_size()
                best_merit = 0
                best_i = None
                best_mapping = None
                for i in random.sample( my_can_expand, samp_size_expand ):
                    new_mapping = mapping.expand( i )
                    if new_mapping.get_out_size() == symbol_count: 
                        continue
                    merit = calc_merit( training_sets, new_mapping, modname, modorder  )
                    if merit > best_merit:
                        best_i = i
                        best_merit = merit
                        best_mapping = new_mapping
                if mpi:
                    best_i, best_j, best_merit, cv_runs = sync_nodes( best_i, None, best_merit, 0 )
                    assert best_j == None
                    best_mapping = mapping.expand( best_i )
                if best_mapping:
                    mapping = best_mapping
                
        step_counter += 1

def sync_nodes( best_i, best_j, best_merit, cv_runs ):
    # Aggregate results from all nodes
    if node_id != 0:
        # Send best i, j, merit to the master
        pypar.send( ( best_i, best_j, best_merit, cv_runs ), 0 )
        # Get back the overall best i, j, merit from the master
        best_i, best_j, best_merit, cv_runs  = pypar.receive( 0 )
    else:
        # I am the master, get results from all other nodes and determine
        # which had the best merit
        for other_node_id in range( 1, nodes ):
            i, j, merit, runs = pypar.receive( other_node_id )
            if merit > best_merit:
                best_i, best_j, best_merit = i, j, merit
            cv_runs += runs
        # Send back the overall bests
        for other_node_id in range( 1, nodes ):
            pypar.send( ( best_i, best_j, best_merit, cv_runs ), other_node_id )
    return best_i, best_j, best_merit, cv_runs

def write_mapping( mapping, fname ):
    """
    Writes mapping from atom symbols to collapsed symbols (NOT from original
    columns anymore!)
    """
    mapping_out = open( fname, 'w' )
    for i in range( mapping.get_in_size() ):
        print >>mapping_out, i, mapping[ i ]
    mapping_out.close()

def all_pairs( n ):
    rval = []
    for i in range( 0, n ):
        for j in range( i + 1, n ):
            rval.append( ( i, j ) )
    return rval

def calc_merit( training_sets, mapping, modname, modorder ):
    # Apply mapping to strings
    training_sets = [ [ mapping.translate( s ) for s in strings ] for strings in training_sets ]
    # Cross validate using those strings
    radix = mapping.get_out_size()
    
    if len( training_sets ) == 2:
        pos_strings, neg_strings = training_sets
        model_factory = lambda d0, d1: rp.models.train( modname, modorder, radix, d0, d1 )
        cv_engine = rp.cv.CV( model_factory, pos_strings, neg_strings, fold=fold, passes=passes )
        cv_engine.run()
        # Merit is TP + TN
        return ( cv_engine.cls1.pos / ( len( pos_strings ) * passes ) + cv_engine.cls2.neg / ( len( neg_strings ) * passes ) ) / 2
    elif len( training_sets ) > 2:
        model_factory = lambda d: rp.models.prob_train( modname, modorder, radix, d )    
        cv_engine = rp.cv.MultiCV( model_factory, training_sets, fold=fold, passes=passes )
        cv_engine.run()
        ## print >> sys.stderr, cv_engine.get_summary()
        ## print >> sys.stderr, cv_engine.get_success_rate()
        return cv_engine.get_success_rate()       
    else:
        raise Exception( "No support for '%d' training sets" % len( training_sets ) )

def main():

    global mpi, fold, passes, loo

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    if 1:
        ts_fnames = args
        out_dir = options.out
        align_count, atom_mapping = rp.mapping.alignment_mapping_from_file( file( options.atoms ) )
        if options.mapping: mapping = rp.mapping.second_mapping_from_file( file( options.mapping ), atom_mapping )
        else: mapping = rp.mapping.identity_mapping( atom_mapping.get_out_size() )
        modname = options.model
        modorder = int( options.order )
        if options.fold: fold = int( options.fold )
        if options.passes: fold = int( options.passes)
        loo = bool( options.loo )
        mpi = bool( options.mpi )
        
    #except:
    #    cookbook.doc_optparse.exit()

    run( ts_fnames, out_dir, options.format, align_count, atom_mapping, mapping, modname, modorder )


if __name__ == "__main__": 
    main()
