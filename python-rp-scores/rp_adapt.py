#!/usr/bin/env python

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'.

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
from rp import cv, io

import rp.models.averaging

import rp.mapping

STOP_SIZE = 22

fold = 10
passes = 10

# FIXME: select order dynamically
order = 1
model = rp.models.averaging

def run( pos_file, neg_file, out_dir, format, mapping ):

    # Open merit output
    merit_out = open( os.path.join( out_dir, 'merits.txt' ), 'w' )

    # Read integer sequences
    pos_strings = list( io.get_reader( pos_file, format, None ) )
    neg_strings = list( io.get_reader( neg_file, format, None ) )

    symbol_count = mapping.get_out_size()

    # Collapse
    while symbol_count > STOP_SIZE:

        print >> sys.stderr, "Collapsing from:", symbol_count

        best_mapping = None
        best_merit = 0

        pb = ProgressBar( 0, symbol_count * ( symbol_count - 1 ) / 2, 78 )

        count = 0

        for i in range( 0, symbol_count ):
            for j in range( i + 1, symbol_count ):

                collapsed = mapping.collapse( i, j )
                
                merit = calc_merit( pos_strings, neg_strings, collapsed ) 
                
                if merit > best_merit:
                    best_merit = merit
                    best_mapping = collapsed
                
                count += 1
                pb.update_and_print( count, sys.stderr )
                
        symbol_count -= 1
        
        # Append merit to merit output        
        print >>merit_out, symbol_count, best_merit
        
        print >>sys.stderr, "\nBest Merit %d." % best_merit,
        
        # Write best mapping to a file
        mapping_out = open( os.path.join( out_dir, "%03d.mapping" % symbol_count ), 'w' )
        for symbol in best_mapping.get_table(): print >>mapping_out, symbol
        mapping_out.close()

def calc_merit( pos_strings, neg_strings, mapping ):
    # Apply mapping to strings
    pos_strings = [ mapping.translate( s ) for s in pos_strings ]
    neg_strings = [ mapping.translate( s ) for s in neg_strings ]
    # Cross validate using those strings
    model_factory = lambda d0, d1: model.train( order, mapping.get_out_size(), d0, d1 )
    cv_engine = cv.CV( model_factory, pos_strings, neg_strings, fold=fold, passes=passes )
    cv_engine.run()
    # Merit is TP + TN
    return cv_engine.cls1.pos + cv_engine.cls2.neg

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        pos_fname, neg_fname, out_dir = args
        mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
    except:
        cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), out_dir, options.format, mapping )


if __name__ == "__main__": main()
