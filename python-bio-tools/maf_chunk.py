#!/usr/bin/env python2.3

"""
Read a MAF from stdin and break into a set of mafs containing 
no more than a certain number of columns
"""

usage = "usage: %prog chunk_size out_dir"

import sys, align.maf
from optparse import OptionParser

try:
    import psyco
    psyco.profile()
except:
    pass

INF="inf"

def __main__():

    # Parse command line arguments

    parser = OptionParser( usage=usage )
    ( options, args ) = parser.parse_args()

    chunk_size = int( args[0] )
    out_dir = args[1]

    maf_reader = align.maf.Reader( sys.stdin )

    maf_writer = None

    count = 0
    current_chunk = -1

    chunk_min = INF
    chunk_max = 0

    interval_file = file( "%s/intervals.txt" % out_dir, "w" )	

    for m in maf_reader:
        chunk_min = min( chunk_min, m.components[0].start )
        chunk_max = max( chunk_max, m.components[0].get_end() )
        if not maf_writer or count + m.text_size > chunk_size:
            current_chunk += 1
            if maf_writer: 
                maf_writer.close()
                interval_file.write( "%s %s\n" % ( chunk_min, chunk_max ) )
                chunk_min = INF
                chunk_max = 0
            maf_writer = align.maf.Writer( file( "%s/%09d.maf" % ( out_dir, current_chunk ), "w" ) )
            count = 0
        maf_writer.write( m )
        count += m.text_size

    if maf_writer:
        maf_writer.close()
        interval_file.write( "%s %s\n" % ( chunk_min, chunk_max ) )
        interval_file.close()

if __name__ == "__main__": __main__()
