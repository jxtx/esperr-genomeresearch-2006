#!/usr/bin/env python

import psyco_full

import sys

import ranges, sys
import bx.align.maf
from optparse import OptionParser

def __main__():

    maf_reader = bx.align.maf.Reader( sys.stdin )
    maf_writer = bx.align.maf.Writer( sys.stdout )

    group = sys.argv[1]

    for line, m in zip( open( sys.argv[2] ), maf_reader ):
        if line.strip() == group:
            maf_writer.write( m )

    maf_writer.close()        

if __name__ == "__main__": __main__()
