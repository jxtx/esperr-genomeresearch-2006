#!/usr/bin/env python2.4

"""
Score a set of sequences using a model

usage: %prog data score_matrix out [options]
   -i, --input=FILE:       Input filename (otherwise stdin)
   -d, --output=FILE:      Output filename (otherwise stdout)
   -f, --format=FILE:   Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE:  A mapping (alphabet reduction) to apply to each sequence (optional)
   -M, --model=name:    Name of model to train (default 'standard')
   -G, --mincols=INT:   Minimum number of scorable positions 
   -c, --classnames=A,B: Names to use for the classes in the output. 
"""

from __future__ import division

import pkg_resources
pkg_resources.require( "bx-python" )

import cookbook.doc_optparse
import sys
import traceback

from Numeric import *

import rp.io
import rp.mapping
import rp.models

nan = float( 'nan' )

def run( data_file, out_file, model_fnames, format, mapping, modname, mincols, classnames ):

    max_order = 0
    models = []
    for f in model_fnames:
        model = rp.models.get( modname ).prob_from_file( f )
        radix = model.get_radix()
        order = model.get_order()
        max_order = max( max_order, order )
        models.append( model )

    # Score each
    for string in rp.io.get_reader( data_file, format, mapping ):
        if len( string ) < max_order or sum( string != -1 ) < mincols:
            print >> out_file, "\t".join( ( [ "NA" ] * len( models ) ) + [ "NO DATA" ] )
        else:
            probs = [ model.score( string ) for model in models ]
            best = classnames[ argmax( probs ) ]
            print >> out_file, "\t".join( map( str, list( exp( probs ) ) + [ best ] ) )

def main():

    # Parse command line
    try:
        options, args = cookbook.doc_optparse.parse( __doc__ )
        model_fnames = args
        if options.input:
            in_file = open( options.input )
        else:
            in_file = sys.stdin
        if options.output:
            out_file = open( options.output, "w" )
        else:
            out_file = sys.stdout
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
        mincols = getattr( options, 'mincols' )
        if mincols: mincols = int( mincols )
        if options.classnames:
            classnames = options.classnames.split( "," )
            assert len( classnames ) == len( model_fnames )
        else:
            classnames = map( str, range( len( model_fnames ) ) )
    except:
        cookbook.doc_optparse.exit()

    run( in_file, out_file, model_fnames, options.format, mapping, modname, mincols, classnames )
    out_file.close()

if __name__ == "__main__": main()
