#!/usr/bin/env python2.3

"""
Read a maf and print the text as a fasta file.
"""

from __future__ import division

import textwrap
import sys
import cookbook.doc_optparse
from align import maf
from optparse import OptionParser

def __main__():

    maf_reader = maf.Reader( sys.stdin )

    for i, m in enumerate( maf_reader ):
        for c in m.components:
            print ">%s:%d-%d" % ( c.src, c.start, c.end )
            print_n( c.text, 50 )

def print_n( s, n, f = sys.stdout ):
    p = 0
    while p < len( s ):
        print >> f, s[p:min(p+n,len(s))]
        p += n

if __name__ == "__main__": __main__()
