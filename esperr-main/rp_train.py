#!/usr/bin/env python

"""
Train a discriminating model from two data sets and store it.

usage: %prog pos_data neg_data out [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -r, --radix=N:      Radix
   -o, --order=N:      Order
   -M, --model=name:   Name of model to train (default 'standard')
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import array
import cookbook.doc_optparse
import sys
import traceback

import rp.io
import rp.mapping
import rp.models

def run( pos_file, neg_file, out_file, format, mapping, radix, order, modname ):

    # Read integer sequences
    pos_strings = list( rp.io.get_reader( pos_file, format, mapping ) )
    neg_strings = list( rp.io.get_reader( neg_file, format, mapping ) )

    # Determine radix
    if not radix:
        if mapping: radix = mapping.get_out_size()
        else: radix = max( map( max, pos_strings ) + map( max, neg_strings ) ) + 1
               
    # Build model
    print "about to train"
    model = rp.models.train( modname, order, radix, pos_strings, neg_strings )
    print "trained"

    # Write to out file
    print "about to write"
    model.to_file( out_file )
    print "written"

def main():

    # Parse command line
    try:
        options, args = cookbook.doc_optparse.parse( __doc__ )
        pos_fname, neg_fname, out_fname = args
        order = int( getattr( options, 'order' ) )
        radix = getattr( options, 'radix' )
        format = getattr( options, 'format' )
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
            print "Align count:", align_count, "Mapping: ", mapping
        else:
            mapping = None
    except:
        cookbook.doc_optparse.exit()

    out = open( out_fname, "w" )
    run( open( pos_fname ), open( neg_fname ), out, format, mapping, radix, order, modname )
    out.close()

if __name__ == "__main__": main()
