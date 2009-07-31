#!/usr/bin/env python2.3

"""
Score a set of sequences using a model

usage: %prog data score_matrix out [options]
   -f, --format=FILE:   Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE:  A mapping (alphabet reduction) to apply to each sequence (optional)
   -M, --model=name:    Name of model to train (default 'standard')
   -g, --mingood=FLOAT: Minimum fraction of good columns required to produce a score
   -G, --mincols=INT:   Minimum number of scorable positions 
"""

from __future__ import division

import pkg_resources
pkg_resources.require( "bx-python" )

import bx.cookbook.doc_optparse
import sys
import traceback

from numpy import *

import rp.io
import rp.mapping
import rp.models

nan = float( 'nan' )

def run( data_file, model_file, out_file, format, mapping, modname, mingood, mincols ):

    # Read model
    model = rp.models.get( modname ).from_file( model_file )
    radix = model.get_radix()
    order = model.get_order()

    # Read integer sequences
    strings = rp.io.get_reader( data_file, format, mapping )

    # Score each
    for string in strings:
        if mingood is None and mincols is None:
            score = model.score( string )
        else:
            scores = array( [ nan ] * len( string ), dtype="f" )
            model.score_positions( string, scores )
            goodwords = equal(scores,scores)
            ngood = sum( goodwords )
            putmask( scores, not_equal(scores,scores) , 0 )
            if mingood and ( ngood / ( len( string ) - order ) >= mingood ): 
                score = sum( scores ) / ngood
            elif mincols and ( ngood >= mincols ):
                score = sum( scores ) / ngood
            else:
                score = None
        print >>out_file, score or "NA"

def main():

    # Parse command line
    try:
        options, args = bx.cookbook.doc_optparse.parse( __doc__ )
        data_fname, model_fname, out_fname = args
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
        mingood = getattr( options, 'mingood' )
        if mingood: mingood = float( mingood ) 
        mincols = getattr( options, 'mincols' )
        if mincols: mincols = int( mincols )
    except:
        bx.cookbook.doc_optparse.exception()

    out = open( out_fname, "w" )
    run( open( data_fname ), open( model_fname ), out, options.format, mapping, modname, mingood, mincols )
    out.close()

if __name__ == "__main__": main()
