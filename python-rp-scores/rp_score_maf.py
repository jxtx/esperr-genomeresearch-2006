#!/usr/bin/env python2.3

"""
Score a set of alignments (MAF format) using a model

usage: %prog data score_matrix out [options]
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -w, --window=N:     Size of window to scroll over sequence (default 100)
   -s, --shift=N:      Amount to shift window (deafult 5)
"""

import align.maf
import array
import cookbook.doc_optparse
import seq_numarray
import sys
import traceback

import rp.io 
import rp.mapping
import rp.models.standard

def run( data_file, model_file, out_file, mapping, window, shift ):

    # Read model
    model = rp.models.standard.from_file( model_file )
    radix = model.get_radix()

    # Open maf file
    mafs = align.maf.Reader( data_file )

    # Score each alignment
    for maf in mafs:
        ints = rp.mapping.DNA.translate_list( [ c.text for c in maf.components ] )
        if mapping: ints = mapping.translate( ints )
        score_windows( maf, array.array( 'i', list( ints ) ), model, out_file, window, shift )

def score_windows( maf, string, model, out, window, shift ):
    if maf.text_size < window: return
    half_window = window // 2
    text = maf.components[0].text
    # Output position is middle of window
    abs_pos = maf.components[0].start + ( half_window - text.count( '-' ) ) 
    last_pos = None
    for i, c in enumerate( text ):
        if i + window >= len( text ): break
        if c != '-': abs_pos += 1
        if abs_pos % shift == 0:
	    score = model.score( string, i, window )
	    if score is not None:
		if abs_pos == last_pos: continue
		print >>out, abs_pos, score
		last_pos = abs_pos

def getopt( options, name, default ):
    v = getattr( options, name )
    if v is None: return default
    return v

def main():

    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    #try:
    if 1:
        data_fname, model_fname, out_fname = args
        window = int( getopt( options, 'window', 100 ) )
        shift = int( getopt( options, 'shift', 5 ) )
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
    #except:
    #    cookbook.doc_optparse.exit()

    out = open( out_fname, "w" )
    run( open( data_fname ), open( model_fname ), out, mapping, window, shift )
    out.close()

if __name__ == "__main__": main()
