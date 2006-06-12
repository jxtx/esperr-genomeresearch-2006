#!/usr/bin/env python2.4

"""
Read a (hopefully small) MAF file and a set of intervals, and write a new
maf with 1 block for each interval (fusing any adjacent blocks in the input, 
and joing non adjacent blocks with a pad of all gaps)
"""

import sys
import bx.align.maf
import bx.align as align
from bx import intervals
from itertools import *

def main():

    intersecters = dict()    

    nspecies = None
    species = None

    for maf in bx.align.maf.Reader( open( sys.argv[1] ) ):
        rc = maf.components[ 0 ]
        if nspecies:
            assert len( maf.components ) == nspecies, "All blocks must have same number of species"
        else:
            nspecies = len( maf.components )
            species = [ c.src.split( "." )[0] for c in maf.components ]
        if not rc.src in intersecters: 
            intersecters[ rc.src ] = intervals.Intersecter()
        intersecters[ rc.src ].add_interval( intervals.Interval( int( rc.start ), int( rc.end ), maf ) )

    out = bx.align.maf.Writer( sys.stdout )
    
    for line in open( sys.argv[2] ):
        fields = line.split()
        src, start, end = fields[0], int( fields[1] ), int( fields[2] )
        intersection = None
        if src in intersecters:
            intersection = intersecters[ src ].find( start, end )
        if not intersection:
            a = align.Alignment()
            for i, name in enumerate( species ):
                text = ""
                size = len( text )
                c = align.Component( name + ".fake", 0, size, "?", size, text )
                a.add_component( c )
            intersection = [ a ]
        else:
            intersection.sort()
            intersection = [ i.value for i in intersection ]
            # Try to fuse
            fuse_list( intersection )
            # Still need to join some
            if len( intersection ) > 1:
                texts = [ [] for i in range( nspecies ) ]
                for block in intersection:
                    for t, c in izip( texts, block.components ):
                        t.append( c.text )
                a = align.Alignment()
                for i, name in enumerate( species ):
                    text = "----------".join( texts[i] )
                    size = len( text ) - text.count( "-" )
                    c = align.Component( name + ".fake", 0, size, "?", size, text )
                    a.add_component( c )
                intersection = [ a ]
        # Write
        out.write( intersection[0] )

    out.close()
        
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
        if c1.strand != c2.strand: return None
        if c1.end != c2.start: return None
    # Try to fuse:
    n = copy.deepcopy( m1 )
    for c1, c2 in izip( n.components, m2.components ):
        c1.text += c2.text
        c1.size += c2.size
    n.text_size = len( n.components[0].text )
    return n
    
if __name__ == "__main__":
    main()
