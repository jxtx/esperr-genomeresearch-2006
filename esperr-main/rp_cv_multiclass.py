#!/usr/bin/env python2.4

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
   -p, --print_results: Print Results instead of summary
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

def run( ts_fnames, format, mapping, radix, orders, modname, fold, loo, print_results ):

    # Read integer sequences
    radix = 0
    training_sets = []
    for f in ts_fnames:
        strings = []
        for s in rp.io.get_reader( open( f ), format, mapping ):
            if sum( s != -1 ) >= 50:
                strings.append( s )
                radix = max( radix, max( s ) + 1 )
        training_sets.append( strings )
                
    # Determine radix
    if mapping: 
        radix = mapping.get_out_size()

    # Cross validate for various orders
    for order in orders:
        model_factory = lambda d0: rp.models.prob_train( modname, order, radix, d0 )
        if loo: passes = 1
        else: passes = 5 
        cv_engine = rp.cv.MultiCV( model_factory, training_sets, fold, passes, loo, keep_results=print_results )
        start_time = time.time()
        cv_engine.run()
        seconds = time.time() - start_time

        if print_results:
            for r, c in zip( cv_engine.scores, cv_engine.classes ):
                print "\t".join( map( str, list( r ) + list( c ) ) )
        else:
            print cv_engine.get_success_rate()*100, cv_engine.get_summary(), seconds
        
def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    #try:
    if 1:
        ts_fnames = args
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
        print_results = bool( options.print_results )
    #except:
    #    cookbook.doc_optparse.exit()

    run( ts_fnames, options.format, mapping, radix, orders, modname, fold, loo, print_results )


if __name__ == "__main__": main()
