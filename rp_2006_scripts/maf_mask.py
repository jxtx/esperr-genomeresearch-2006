#!/usr/bin/env python2.4

"""
"""

import psyco_full

import cookbook.doc_optparse

import bx.align.maf
from bx.align import core
from bx import intervals
from bx.bitset_builders import *
import sys

def coord_to_col( start, text, pos ):
    col = 0
    while start < pos:
        assert 0 <= col < len( text ), "0 <= %d < %d" % ( col, len( text ) )
        if text[col] != '-': 
            start += 1
        col += 1 
    while text[col] == '-':
        col += 1
    return col

def __main__():

    # Parse Command Line

    maf_fname, mask_fname = sys.argv[1], sys.argv[2] 

    # Load Intervals

    bitsets = binned_bitsets_from_file( open( mask_fname ) )

    # Start MAF on stdout

    out = bx.align.maf.Writer( sys.stdout )

    # Iterate over input MAF

    for maf in bx.align.maf.Reader( open( maf_fname ) ):
        rc = maf.components[ 0 ]
        chr = rc.src.split( '.' )[1]
        if chr in bitsets:
            positions_to_mask = [ 0 ] * len( rc.text )
            for i in range( rc.forward_strand_start, rc.forward_strand_end ):
                if bitsets[ chr ][ i ]:
                    if rc.strand == "-":
                        i = ( rc.src_size - i - 1 )
                    col = coord_to_col( rc.start, rc.text, i )
                    positions_to_mask[col] = 1
                    
            in_gap = False        
            mask_left = False
            for i in range( len( rc.text ) ):
                if rc.text[i] == '-': 
                    if not in_gap:
                        in_gap = True
                        gap_start = i
                        mask_left = ( i > 0 ) or positions_to_mask[i-1] 
                elif in_gap:
                    if positions_to_mask[i]:
                        for j in range( gap_start, i ):
                            positions_to_mask[j] = 1
                    in_gap = False        
            if in_gap and mask_left:
                for j in range( gap_start, len( rc.text ) ):
                    positions_to_mask[j] = 1
            for comp in maf.components:                
                l = list( comp.text )
                for i in range( len( comp.text ) ):
                    if positions_to_mask[ i ]:
                        l[ i ] = 'X'
                comp.text = ''.join( l )

        out.write( maf )    
         
    out.close()

if __name__ == "__main__": __main__()
