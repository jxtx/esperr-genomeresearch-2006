#!/usr/bin/env python2.3

"""
Score a set of alignments (MAF format) using a model

usage: %prog data score_matrix out [options]
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -w, --window=N:     Size of window to scroll over sequence (default 100)
   -s, --shift=N:      Amount to shift window (deafult 5)
"""

import align.maf
import alphabet
import array
import cookbook.doc_optparse
import seq_numarray
import sys
import traceback

from rp import io, standard_model

def run( data_file, model_file, out_file, mapping, window, shift ):

    # Read model
    model = standard_model.from_file( model_file )
    order = model.get_order()
    radix = model.get_radix()

    # Open maf file
    mafs = align.maf.Reader( data_file )

    # Score each alignment
    for maf in mafs:
	ints = seq_numarray.DNA.translate_alignment( [ c.text for c in maf.components ] )
	if mapping: ints = mapping.translate( ints )
	score_windows( maf, array.array( 'i', list( ints ) ), model, out, window, shift )

def score_windows( maf, string, model, out, window, shift ):

    abs_pos = maf.components[0].start
    text = maf.components[0].text
    for i, c in enumerate( text ):
	if i + window >= len( text ): break
	if c != '-': abs_pos += 1
	if abs_pos % shift == 0:
	    print >>out, abs_pos, model.score( string, i, window )

def main():
    
    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    try:
        data_fname, model_fname, out_fname = args
	window = int( getattr( options, 'window', 100 ) )
	shift = int( getattr( options, 'shift', 5 ) )
        if options.mapping: 
            mapping = alphabet.Mapping( file( options.mapping ) )
        else: 
            mapping = None
    except:
        cookbook.doc_optparse.exit()
    
    out = open( out_fname, "w" )
    run( open( data_fname ), open( model_fname ), out, mapping, window, shift )
    out.close()

if __name__ == "__main__": main()
