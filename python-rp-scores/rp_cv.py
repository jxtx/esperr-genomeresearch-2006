#!/usr/bin/env python

"""
Use cross validation to evaluate a model for some training data.

usage: %prog pos_data neg_data [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -r, --radix=N:      Radix (optional)
   -o, --orders=N,...:   Orders to cross validate over
"""

import align.maf
import array
import cookbook.doc_optparse
import sys
import traceback

import rp.cv
import rp.io
import rp.mapping
import rp.standard_model

def run( pos_file, neg_file, format, mapping, radix, orders ):

    # Read integer sequences
    pos_strings = list( rp.io.get_reader( pos_file, format, mapping ) )
    neg_strings = list( rp.io.get_reader( neg_file, format, mapping ) )

    # Determine radix
    if not radix:
        if mapping: radix = mapping.get_out_size()
        else: radix = max( map( max, pos_strings ) + map( max, neg_strings ) ) + 1

    # Cross validate for various orders
    for order in orders:
        model_factory = lambda d0, d1: rp.standard_model.train( order, radix, d0, d1 )
        cv_engine = rp.cv.CV( model_factory, pos_strings, neg_strings )
        cv_engine.run()

        print "Order:", order
        print "          +   ~+   ~-    -"
        print "Reg:  ", cv_engine.cls1
        print "AR:   ", cv_engine.cls2
        print "Rate: ", cv_engine.get_success_rate()

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        pos_fname, neg_fname = args
        orders = map( int, getattr( options, 'orders' ).split( ',' ) )
        radix = getattr( options, 'radix', None )
        if options.mapping:
            mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
    except:
        cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), options.format, mapping, radix, orders )


if __name__ == "__main__": main()
