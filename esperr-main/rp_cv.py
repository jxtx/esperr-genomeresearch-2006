#!/usr/bin/env python

"""
Use cross validation to evaluate a model for some training data.

usage: %prog pos_data neg_data [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -r, --radix=N:      Radix (optional)
   -o, --orders=N,...: Orders to cross validate over
   -F, --fold=N:       Fold (default 5)
   -M, --model=name:   Name of model to train (default 'standard')
   -l, --loo:          Use leave-one-out cross validation (fold is ignored in this case) 
"""

import array
import cookbook.doc_optparse
import sys
import traceback
import time

import rp.cv
import rp.io
import rp.mapping
import rp.models

default_fold = 5

def run( pos_file, neg_file, format, mapping, radix, orders, modname, fold, loo ):

    # Split up 

    # Read integer sequences
    pos_strings = list( rp.io.get_reader( pos_file, format, mapping ) )
    neg_strings = list( rp.io.get_reader( neg_file, format, mapping ) )

    # Determine radix
    if not radix:
        if mapping: radix = mapping.get_out_size()
        else: radix = max( map( max, pos_strings ) + map( max, neg_strings ) ) + 1

    print "Order     TP  ~TP  ~FN   FN   FP  ~FP  ~TN   TN       %    time"

    # Cross validate for various orders
    for order in orders:
        model_factory = lambda d0, d1: rp.models.train( modname, order, radix, d0, d1 )
        if loo: passes = 1
        else: passes = 5 
        cv_engine = rp.cv.CV( model_factory, pos_strings, neg_strings, fold, passes, loo )
        start_time = time.time()
        cv_engine.run()
        seconds = time.time() - start_time

        print "%5d  " % order,
        print cv_engine.cls1, cv_engine.cls2,
        print "  %2.2f    %2.2f" % ( cv_engine.get_success_rate()*100, seconds )
        
def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    #try:
    if 1:
        pos_fname, neg_fname = args
        orders = map( int, getattr( options, 'orders' ).split( ',' ) )
        radix = getattr( options, 'radix', None )
        if radix: radix = int( radix )
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        if options.fold:
            fold = int( options.fold )
        else:
            fold = default_fold
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
        loo = bool( options.loo )    
    #except:
    #    cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), options.format, mapping, radix, orders, modname, fold, loo )


if __name__ == "__main__": main()
