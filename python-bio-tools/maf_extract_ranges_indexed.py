#!/usr/bin/env python2.3

"""
Reads a list of intervals and a maf. Produces a new maf containing the
blocks or parts of blocks in the original that overlapped the intervals

NOTE: May write the same block more than once, fixme?

usage: %prog maf_file index_file [options] < interval_file
   -m, --mincols=10: Minimum length (columns) required for alignment to be output
   -c, --chop:       Should blocks be chopped to only portion overlapping (no by default)
   -s, --src=s:      Use this src for all intervals
   -p, --prefix=p:   Prepend this to each src before lookup
"""

import psyco_full

import cookbook.doc_optparse

import align.maf
import intervals
import misc
import ranges
import sys


def __main__():

    # Parse Command Line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        maf_files = args
        if options.mincols: mincols = int( options.mincols )
        else: mincols = 10
        if options.src: fixed_src = options.src
        else: fixed_src = None
        if options.prefix: prefix = options.prefix
        else: prefix = None
        chop = bool( options.chop )
    except:
        cookbook.doc_optparse.exit()

    # Open indexed access to mafs
    indexes = [ align.maf.Indexed( maf_file, maf_file + ".index" ) for maf_file in maf_files ]

    # Start MAF on stdout

    out = align.maf.Writer( sys.stdout )

    # Iterate over input ranges 

    for line in sys.stdin:
        fields = line.split()
        if fixed_src:
            src, start, end = fixed_src, int( fields[0] ), int( fields[1] )
        else:
            src, start, end = fields[0], int( fields[1] ), int( fields[2] )
            if prefix: src = prefix + src
        # Find overlap with reference component
        blocks = []
        for index in indexes: blocks += index.get( src, start, end )
        # Write each intersecting block
        if chop:
            for block in blocks: 
                ref = block.get_component_by_src( src )
                slice_start = max( start, ref.start )
                slice_end = min( end, ref.end )
                sliced = block.slice_by_component( ref, slice_start, slice_end ) 
                good = True
                for c in sliced.components: 
                    if c.size < 1: 
                        good = False
                if good and sliced.text_size > mincols: 
                    out.write( sliced )
        else:
            for block in blocks:
                out.write( block )
         
    # Close output MAF

    out.close()

if __name__ == "__main__": __main__()
