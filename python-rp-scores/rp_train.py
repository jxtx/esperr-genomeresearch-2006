#!/usr/bin/env python2.3

"""
Train a discriminating model from two data sets and store it.

usage: %prog pos_data neg_data out [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -r, --radix=N:      Radix
   -o, --order=N:      Order
"""

import align.maf
import alphabet
import array
import cookbook.doc_optparse
import sys
import traceback

from rp import io, standard_model

def main():
    
    # Parse command line
    try:
        options, args = cookbook.doc_optparse.parse( __doc__ )
        pos_fname, neg_fname, out_fname = args
	order = int( getattr( options, 'order' ) )
	radix = getattr( options, 'radix', None )
        if options.mapping: 
            mapping = alphabet.Mapping( file( options.mapping ) )
        else: 
            mapping = None
    except:
        cookbook.doc_optparse.exit()

    # Read integer sequences
    pos_strings = list( io.get_reader( open( pos_fname ), options.format, mapping ) )
    neg_strings = list( io.get_reader( open( neg_fname ), options.format, mapping ) )

    # Determine radix
    if not radix:
	if mapping: radix = mapping.symbol_count
        else: radix = max( map( max, pos_strings ) + map( max, neg_strings ) ) + 1

    # Build model
    model = standard_model.train( order, radix, pos_strings, neg_strings )

    # Write it to a file
    out = open( out_fname, "w" )
    model.to_file( out )
    out.close()

if __name__ == "__main__": main()
