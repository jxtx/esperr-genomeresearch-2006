#!/usr/bin/env python2.3

"""
Score a set of int seqs, using an accompanying BED file for position information

usage: %prog data.ints data.bed score_matrix out [options]
   -w, --window=N:     Size of window to scroll over sequence (default 100)
   -s, --shift=N:      Amount to shift window (deafult 5)
"""

from itertools import *

import align.maf
import array
import cookbook.doc_optparse
import seq_numarray
import sys
import traceback

import rp.io 
import rp.mapping
import rp.models.standard

def run( data_file, bed_file, model_file, out_file, window, shift ):

    # Read model
    model = rp.models.standard.from_file( model_file )
    radix = model.get_radix()

    # Open maf file
    data = imap( lambda x: map( int, x.split() ), data_file )

    # Score each alignment
    for ints, bed_line in zip( data, bed_file ):
        chr, start, end = bed_line.split()
        # Wiggle header
        if '.' in chr: chr = chr.split('.')[1]
        print >>out_file, "variableStep chrom=" + chr
        score_windows( int( start ), array.array( 'i', list( ints ) ), model, out_file, window, shift )

def score_windows( start, string, model, out, window, shift ):
    half_window = window // 2
    # Output position is middle of window
    abs_pos = start + half_window 
    last_pos = None
    for i in range( len( string ) - window ):
        abs_pos += 1
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
        data_fname, bed_fname, model_fname, out_fname = args
        window = int( getopt( options, 'window', 100 ) )
        shift = int( getopt( options, 'shift', 5 ) )
    #except:
    #    cookbook.doc_optparse.exit()

    out = open( out_fname, "w" )
    run( open( data_fname ), open( bed_fname ), open( model_fname ), out, window, shift )
    out.close()

if __name__ == "__main__": main()
