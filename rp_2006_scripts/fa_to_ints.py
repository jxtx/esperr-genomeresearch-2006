#!/usr/bin/env python2.4

import sys
from bx import seqmapping

def read_fa( f ):
    last = None
    last_header = None
    for line in f:
        if line.startswith( ">" ):
            if last_header:
                yield last_header, ''.join( last )
            last_header = line[1:]
            last = []
        else:
            last.append( line.strip() )
    if last_header:
        yield last_header, ''.join( last )

def main():
    for h, s in read_fa( sys.stdin ):
        print ' '.join( map( str, seqmapping.DNA.translate( s ) ) )

if __name__ == "__main__":
    main()
