#!/usr/bin/env python2.3

"""
Read a MAF from stdin and break into a set of mafs containing 
no more than a certain number of columns
"""

usage = "usage: %prog chunk_size out_dir"

import sys
from bx import align.maf
from optparse import OptionParser
import psyco_full
import random

INF="inf"

def __main__():

    # Parse command line arguments

    parser = OptionParser( usage=usage )
    parser.add_option( "--prob", action="store", default=None, type="float", 
                       help="Probability of writing a given chunk" )
    
    ( options, args ) = parser.parse_args()

    chunk_size = int( args[0] )
    out_dir = args[1]
    prob = options.prob

    maf_reader = align.maf.Reader( sys.stdin )

    maf_writer = None

    count = 0
    current_chunk = -1

    chunk_min = INF
    chunk_max = 0

    write_current_chunk = True

    interval_file = file( "%s/intervals.txt" % out_dir, "w" )	

    for m in maf_reader:
        chunk_min = min( chunk_min, m.components[0].start )
        chunk_max = max( chunk_max, m.components[0].end )
        if not maf_writer or count + m.text_size > chunk_size:
            current_chunk += 1
            # Finish the last chunk            
            if maf_writer: 
                maf_writer.close()
                interval_file.write( "%s %s\n" % ( chunk_min, chunk_max ) )
                chunk_min = INF
                chunk_max = 0
            # Decide if the new chunk will be written     
            if prob: write_current_chunk = bool( random.random() <= prob )
            else: write_current_chunk = True
            if write_current_chunk:
                maf_writer = align.maf.Writer( file( "%s/%09d.maf" % ( out_dir, current_chunk ), "w" ) )
            else:
                maf_writer = None
            count = 0
        if maf_writer: maf_writer.write( m )
        #count += m.text_size
        count += m.components[0].size
    
    if maf_writer:
        maf_writer.close()
        interval_file.write( "%s %s\n" % ( chunk_min, chunk_max ) )
        interval_file.close()

if __name__ == "__main__": __main__()
