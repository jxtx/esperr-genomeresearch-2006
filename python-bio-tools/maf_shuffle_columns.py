#!/usr/bin/env python2.3

import psyco_full

import sys

import ranges, sys
import align
from align import maf

def __main__():

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:
    
        align.shuffle_columns( m )

        maf_writer.write( m )

if __name__ == "__main__": __main__()
