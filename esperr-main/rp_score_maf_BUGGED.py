#!/usr/bin/env python2.4

"""
Score a set of alignments (MAF format) using a model

usage: %prog data score_matrix out [options]
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -M, --model=MODEL:  Name of model to use
   -w, --window=N:     Size of window to scroll over sequence (default 100)
   -s, --shift=N:      Amount to shift window (deafult 5)
   -b, --low=N:        Truncate to this minimum score
   -e, --high=N:       Truncate to this maximum score
   -r, --reorder=0,2,1:Reorder the species in each block before scoring.
"""

from __future__ import division

try: 
    import psyco
    psyco.full()
except: 
    pass

from Numeric import *

import pkg_resources
pkg_resources.require( "bx-python" )

import bx.align.maf
import cookbook.doc_optparse
import sys
import traceback

import rp.io 
import rp.mapping
import rp.models

def run( data_file, modname, model_file, out_file, mapping, window, shift, low, high, reorder ):

    # Read model
    model = rp.models.get( modname ).from_file( model_file )
    radix = model.get_radix()

    # Open maf file
    mafs = bx.align.maf.Reader( data_file )

    # Score each alignment
    for i, maf in enumerate( mafs ):
        if reorder: components = [ maf.components[ i ] for i in reorder ]
        else: components = maf.components
        ints = rp.mapping.DNA.translate_list( [ c.text for c in components ] )
        if mapping: ints = mapping.translate( ints )
        # print i
        score_windows( maf, ints, model, out_file, window, shift, low, high )

def score_windows( maf, string, model, out, window, shift, low, high ):
    if maf.text_size < window: return
    half_window = window // 2
    rc = maf.components[0] 
    text = rc.text
    # Output position is middle of window
    abs_pos = rc.start + ( half_window - text.count( '-', 0, half_window ) ) 
    last_pos = None
    chrom = rc.src
    if '.' in chrom: chrom = chrom.split('.')[1]
    scores = array( [ float("nan") ] * len( text ), typecode="f" )
    model.score_positions( string, scores )
    # Build cumulative sum of scores AND of number of good words per window (note: nan!=nan)
    goodwords = cumsum( equal(scores,scores) )
    putmask( scores, not_equal(scores,scores) , 0 )
    old_scores = scores
    scores = cumsum( scores.astype( Float64 ) )
    import pickle
    f = open( "foo.bin", "w" )
    pickle.dump( old_scores, f )
    f.close()
    need_header = True
    for i, c in enumerate( text ):
        if i + window >= len( text ): break
        if c != '-': 
            abs_pos += 1
        if abs_pos % shift == 0:
            ngood = goodwords[i+window-1]
            if i > 0: 
                ngood -= goodwords[i-1]
            if ngood < 1:
               if abs_pos != last_pos:
                   need_header = True
            elif abs_pos == last_pos:
                pass
            else:
                sumscore = scores[i+window-1]
                if i > 0: sumscore -= scores[i-1]
                score = sumscore / ngood
                #print "error: ", sum( old_scores[:i+window-1] ), sum( abs( old_scores[:i+window-1] ) )
                #print i, abs_pos, score, sumscore, ngood, model.score( string, i, window )
                #assert round( score, 5 ) == round( model.score( string, i, window ), 5 )
                if score > high: score = high
                elif score < low: score = low
                if need_header:
                    print >>out, "fixedStep chrom=%s start=%d step=%d" % ( chrom, abs_pos, shift )
                    need_header = False
                print >>out, abs_pos, round( score, 6 )
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
        low = float( getopt( options, 'low', -1.0 ) )
        high = float( getopt( options, 'high', 1.0 ) )
        if options.mapping:
            align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
        else:
            mapping = None
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        reorder = getopt( options, 'reorder', None )
        if reorder: reorder = map( int, reorder.split( ',' ) )
    #except:
    #    cookbook.doc_optparse.exit()

    out = open( out_fname, "w" )
    run( open( data_fname ), modname, open( model_fname ), out, mapping, window, shift, low, high, reorder )
    out.close()

if __name__ == "__main__": main()
