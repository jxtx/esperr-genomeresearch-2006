#!/usr/bin/env python2.3

import sys

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
    parser.add_option( "--component_count", action="store", default=None, type="int", help="" )
    parser.add_option( "--min_cols", action="store", default=None, type="int", help="" )
    parser.add_option( "-e", "--expr", action="store", default=None )

    ( options, args ) = parser.parse_args()

    component_count = options.component_count
    min_cols = options.min_cols
    expr = options.expr

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:

        if component_count and len( m.components ) != component_count: continue
        if min_cols and m.text_size < min_cols: continue
        if expr and not bool( eval( expr, { "m": m, "maf": m } ) ): continue

        maf_writer.write( m )

if __name__ == "__main__": __main__()
