#!/usr/bin/env python2.3

"""
For each block in a maf file (read from stdin) write a sequence of ints 
corresponding to the columns of the block after applying the provided mapping.

The 'correct' number of species is determined by the mapping file, blocks not having
this number of species will be ignored.

usage: %prog mapping_file
"""

from __future__ import division

import psyco_full

from numarray import *

import align.maf
import alphabet
import seq_numarray
import string
import sys

def main():

    alpha_map = alphabet.Mapping( file( sys.argv[1] ) )
    symbol_count = alpha_map.symbol_count
    align_count = alpha_map.align_count

    for maf in align.maf.Reader( sys.stdin ):
        # Translate alignment to ints
        int_seq = seq_numarray.DNA.translate_alignment( [ c.text for c in maf.components ] )
        # Apply mapping 
        int_seq = alpha_map.translate( int_seq )
        # Write ints separated by spaces
        for i in int_seq: print foo[i],
        print

if __name__ == "__main__": main()
