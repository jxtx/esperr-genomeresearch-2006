#!/usr/bin/env python2.4

import sys

desired_build = sys.argv[1]
window_size = int( sys.argv[2] )

for line in sys.stdin:
    name, build, chrom, size, fname = line.split( "," )
    size = int( size )
    if build == desired_build:
        start = 0
        while 1:
            if start + window_size > size:
                print chrom, start, start + size
                break
            else:
                print chrom, start, start + window_size
                start += window_size
                
