#!/usr/bin/env python2.3

"""
Reads a list of intervals and a maf. Produces a new maf containing the
portions of the original that overlapped the intervals

usage: %prog interval_file refindex [options] < maf_file
   -m, --mincols=10: Minimum length (columns) required for alignment to be output
"""

import psyco_full

import cookbook.doc_optparse

from bx import align.maf
from bx import intervals
import sys


def __main__():

    # Parse Command Line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        range_filename = args[ 0 ]
        refindex = int( args[ 1 ] )
        if options.mincols: mincols = int( options.mincols )
        else: mincols = 10
    except:
        cookbook.doc_optparse.exit()

    # Load Intervals

    intersecter = intervals.Intersecter()
    for line in file( range_filename ):
        fields = line.split()
        intersecter.add_interval( intervals.Interval( int( fields[0] ), int( fields[1] ) ) )

    # Start MAF on stdout

    out = align.maf.Writer( sys.stdout )

    # Iterate over input MAF

    for maf in align.maf.Reader( sys.stdin ):
        ref_component = maf.components[ refindex ]
        # Find overlap with reference component
        intersections = intersecter.find( ref_component.start, ref_component.end )
        # Keep output maf ordered
        intersections.sort()
        # Write each intersecting block
        for interval in intersections: 
            start = max( interval.start, ref_component.start )
            end = min( interval.end, ref_component.end )
            sliced = maf.slice_by_component( refindex, start, end ) 
            good = True
            for c in sliced.components: 
                if c.size < 1: 
                    good = False
            if good and sliced.text_size > mincols: out.write( sliced )
         
    # Close output MAF

    out.close()

if __name__ == "__main__": __main__()
