#!/usr/bin/env python2.3

"""
Read a MAF from standard input and print average GC content of each alignment

usage: %prog [options]
"""

from __future__ import division

import cookbook.doc_optparse
import sys

import psyco_full

from align import maf


def __main__():

    # Parse command line arguments
    # options, args = cookbook.doc_optparse.parse( __doc__ )
    
    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader:
        match = 0
        for i in range( 0, m.text_size ):
            if m.components[0].text[i].lower() == m.components[1].text[i].lower():
                match += 1
        print match / m.text_size


if __name__ == "__main__": __main__()
