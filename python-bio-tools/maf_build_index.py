#!/usr/bin/env python2.3

"""
Build an index file for a set of MAF alignment blocks
"""

import psyco_full

import interval_index_file
import sys

from align import maf

def main():

    maf_file = sys.argv[1]

    if len( sys.argv ) > 2: index_file = sys.argv[2]
    else: index_file = maf_file + ".index" 

    maf_reader = maf.Reader( open( maf_file ) )

    indexes = interval_index_file.Indexes()

    # Need to be a bit tricky in our iteration here to get the 'tells' right
    while 1:
        pos = maf_reader.file.tell()
        block = maf_reader.next()
        if block is None: break
        for c in block.components:
            indexes.add( c.src, c.start, c.end, pos )

    out = open( index_file, 'w' )
    indexes.write( out )
    out.close()

if __name__ == "__main__": main()
