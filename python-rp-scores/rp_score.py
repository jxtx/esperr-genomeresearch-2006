#!/usr/bin/env python2.3

"""
Score a set of sequences using a model

usage: %prog data score_matrix out [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -w, --window=N:     Size of window to scroll over sequence
   -s, --shift:        Amount to shift window
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
        data_fname, model_fname, out_fname = args
        if options.mapping: 
            mapping = alphabet.Mapping( file( options.mapping ) )
        else: 
            mapping = None
    except:
        cookbook.doc_optparse.exit()

    # Read model
    model = standard_model.from_file( open( model_fname ) )
    order = model.get_order()
    radix = model.get_radix()

    # Read integer sequences
    strings = io.get_reader( open( data_fname ), options.format, mapping )

    # Open output file
    out = open( out_fname, "w" )

    # Score each
    for string in strings:
	score = model.score( string, 0, len( string ) )
	print >>out, score
    
    # Close output file
    out.close()

if __name__ == "__main__": main()
