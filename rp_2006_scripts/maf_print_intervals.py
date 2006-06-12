#!/usr/bin/env python2.3

"""
Read a maf from stdin and print the chromosome number or each alignment

usage: %prog refindex [options]
"""

from __future__ import division

import sys
import cookbook.doc_optparse
from bx.align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments
    options, args = cookbook.doc_optparse.parse( __doc__ )

    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader: 
        c = m.components[ 0 ]
        print "%s\t%d\t%d" % ( c.src.split('.')[1], c.forward_strand_start, c.forward_strand_end )

if __name__ == "__main__": __main__()
