#!/usr/bin/env python2.3

# Reads a list of ranges and a maf. Produces a new maf containing the
# portions of the original that overlapped the ranges

import ranges, sys
from align import maf
from optparse import OptionParser

try:
    import psyco
    psyco.profile()
except:
    pass

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "-i", "--refindex", action="store", type="int", default=0, help="Index of reference sequence in alignments" )
    parser.add_option( "-m", "--mincols", action="store", type="int", default=0, help="Minimum length (cols) to output" )

    ( options, args ) = parser.parse_args()

    try:
        range_filename = args[0]
        refindex = options.refindex
        mincols = options.mincols
    except:
        parser.error( "'ranges' option is required" )

    range_reader = ranges.Reader( file( range_filename ) )
    maf_reader = maf.Reader( sys.stdin )

    out = maf.Writer( sys.stdout )

    current_range = range_reader.next()
    current_maf = maf_reader.next()

    overlapping_ranges = []

    while 1:

        if not current_maf:
            break

        current_maf_ref = current_maf.components[ refindex ]
        
        if not current_range:
            write_intersection( current_maf, refindex, overlapping_ranges, out, mincols )
            overlapping_ranges = []
            break

        # Range before maf, advance range
        elif current_range.end < current_maf_ref.start:
            current_range = range_reader.next()

        # Range after maf, advance maf
        elif current_range.start > current_maf_ref.get_end():
            write_intersection( current_maf, refindex, overlapping_ranges, out, mincols )
            overlapping_ranges = []
            current_maf = maf_reader.next()

        # Overlap, save range and advance
        else:
            overlapping_ranges.append( current_range )
            current_range = range_reader.next()

    # Finish reading streams

    while current_maf: current_maf = maf_reader.next()
    while current_range: current_range = range_reader.next()

    range_reader.close()
    maf_reader.close()

    out.close()

def write_intersection( maf, refindex, ranges, out, mincols ):

    for range in ranges:
        start = max( range.start, maf.components[ refindex ].start )
        end = min( range.end, maf.components[ refindex ].start + maf.components[ refindex ].size )
        sliced = maf.slice_by_component( refindex, start, end ) 
        if sliced.text_size > mincols: out.write( sliced )

if __name__ == "__main__": __main__()
