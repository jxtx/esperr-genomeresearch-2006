#!/usr/bin/env python2.3

"""
Read a maf file from stdin and write out a new maf with only blocks having all of
the passed in species, after dropping any other species and removing  columns 
containing only gaps.

usage: %prog species,species2,... < maf
"""

import psyco_full

import bx.align.maf
import bx.align
import copy
import sys

from itertools import *

from cookbook import doc_optparse

def fuse_list( mafs ):
    rval = []
    last = None
    for m in mafs:
        if last is None:
            last = m
        else:
            fused = fuse( last, m )
            if fused:
                last = fused
            else:
                rval.append( last )
                last = m
    if last: rval.append( last )
    return rval

def fuse( m1, m2 ):
    # Check if the blocks are adjacent, return none if not.
    if len( m1.components ) != len( m2.components ): return None
    for c1, c2 in izip( m1.components, m2.components ):
        if c1.src != c2.src: return None
        # SPECIAL CASE: if the right and left components are both fake we'll
        # happily fuse them, no need to check orientation or contiguity
        if not( c1.src.endswith( ".fake" ) and c2.src.endswith( ".fake" ) ):
            if c1.strand != c2.strand: return None
            if c1.end != c2.start: return None
    # Try to fuse:
    n = copy.deepcopy( m1 )
    for c1, c2 in izip( n.components, m2.components ):
        c1.text += c2.text
        c1.size += c2.size
    n.text_size = len( n.components[0].text )
    return n

def _remove_all_gap_columns( texts ):
    """
    Remove any columns containing only gaps from alignment texts
    """
    seqs = [ list( t ) for t in texts ]
    i = 0
    text_size = len( texts[0] )
    while i < text_size:
        all_gap = True
        for seq in seqs:
            if seq[i] not in ( '-', '#', '*', '=', 'X', '@' ): 
                all_gap = False
        if all_gap:
            for seq in seqs: 
                del seq[i]
            text_size -= 1
        else:
            i += 1
    return [ ''.join( s ) for s in seqs ]
    
def remove_all_gap_columns( comps ):
    texts = [ c.text for c in comps ]
    texts = _remove_all_gap_columns( texts )
    for c, t  in zip( comps, texts ):
        c.text = t
    return comps


def main():

    species = []
    optional = []
    
    for s in sys.argv[1].split( ',' ):
        if s.endswith( '*' ):
            species.append( ( s[:-1], True ) )
        else:
            species.append( ( s, False ) )

    maf_reader = bx.align.maf.Reader( sys.stdin )
    maf_writer = bx.align.maf.Writer( sys.stdout )
    
    blocks = []

    for m in maf_reader:        
        new_components = []
        good = True
        for s, opt in species:
            comp = m.get_component_by_src_start( s )
            if not comp and opt:
                fake_len = len( new_components[0].text )
                new_components.append( bx.align.Component( s + ".fake", 0, fake_len, "?", fake_len, "*" * fake_len ) )
            elif not comp:
                good = False
                break
            else:
                new_components.append( comp )
                
        if good:
            m.components = remove_all_gap_columns( new_components )
            blocks.append( m )
            
    for m in fuse_list( blocks ):
        maf_writer.write( m )
        
    maf_reader.close()
    maf_writer.close()

if __name__ == "__main__": 
    main()
