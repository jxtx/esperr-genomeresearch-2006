#!/usr/bin/env python2.3

"""
usage: %prog motif_len mapping_file
"""

from __future__ import division

import psyco_full

import bx.align.maf
from bx import seqmapping
import string
import sys
from Numeric import *

def main():

    word_length = int( sys.argv[1] )

    # Read mapping
    radix = 0
    mapping = dict()
    for line in open( sys.argv[2] ):
        col, val = line.split()
        val = int( val ) 
        radix = max( radix, val + 1 )
        if val in mapping:
            mapping[val].append(col)
        else:
            mapping[val] = [ col ]

    row = 0

    for i in range( radix ):
        for j in range( radix ):
            for k in range( radix ):
                print "%d\t%d %d %d" % ( row, k, j, i )
                # Triple check
                word = i, j, k
                index = 0
                factor = 1
                for x in range( len( word ) ):
                    letter = word[ len(word) - 1 - x ]
                    index += letter * factor
                    factor *= radix
                assert index == row
                row += 1

if __name__ == "__main__": main()
