#!/usr/bin/env python

"""
Train a model from two mafs and a mapping, then dump the score of each 
alignment word in the training data (in order of occurence)

usage: %prog pos_data neg_data [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -o, --order=N:      Order
   -M, --model=name:   Name of model to train (default 'standard')
"""

from numpy import *

import cookbook.doc_optparse
import sys
import traceback
import itertools

import bx.align.maf

import rp.io
import rp.mapping
import rp.models

from sets import Set
set = Set

def run( pos_file, neg_file, format, mapping, order, modname ):

    pos_blocks = list( bx.align.maf.Reader( pos_file ) )
    neg_blocks = list( bx.align.maf.Reader( neg_file ) )
    
    # Read integer sequences
    pos_strings = list( [ mapping.translate( rp.mapping.DNA.translate_list( [ c.text for c in block.components ] ) ) for block in pos_blocks ] )
    neg_strings = list( [ mapping.translate( rp.mapping.DNA.translate_list( [ c.text for c in block.components ] ) ) for block in neg_blocks ] )
    
    # Determine radix
    radix = max( map( max, pos_strings ) + map( max, neg_strings ) ) + 1
               
    # Build model
    model = rp.models.train( modname, order, radix, pos_strings, neg_strings )

    # Find all words in the training data
    words = words_from_blocks( itertools.chain( pos_blocks, neg_blocks ), order )
    for word, count in words:
        ints = rp.mapping.DNA.translate_list( word )
        ints = mapping.translate( ints )
        assert len( ints ) == order + 1
        scores = array( [ float("nan") ] * len( ints ), typecode="f" )
        model.score_positions( ints, scores )
        print "%s\t%d\t%0.6f" % ( "|".join( word ), count, scores[-1] )
    

def words_from_blocks( blocks, order ):
    all_words = []
    all_words_2 = dict()
    for block in blocks:
        for i in range( 0, block.text_size - order ):
            rows = []
            for c in block.components:
                rows.append( c.text[i:i+order+1].upper() )
            rows = tuple( rows )
            if rows not in all_words_2:
                all_words.append( rows )
                all_words_2[ rows ] = 1
            else:
                all_words_2[ rows ] += 1
    return [ ( word, all_words_2[word] ) for word in all_words ]

def main():

    # Parse command line
    #try:
    if 1:
        options, args = cookbook.doc_optparse.parse( __doc__ )
        pos_fname, neg_fname = args
        order = int( getattr( options, 'order' ) )
        format = getattr( options, 'format' )
        modname = getattr( options, 'model' )
        if modname is None: modname = 'standard'
        align_count, mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ) )
    #except:
    #    cookbook.doc_optparse.exit()

    run( open( pos_fname ), open( neg_fname ), format, mapping, order, modname )

if __name__ == "__main__": main()
