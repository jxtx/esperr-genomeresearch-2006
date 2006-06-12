#!/usr/bin/env python2.3

"""
Read a set of ranges and a nib file, print portions of nib overlapping
those ranges to stdout

usage: %prog nib_dir < range_file
"""

from __future__ import division

import cookbook.doc_optparse
import bx.seq.nib
import string
import sys

def __main__():

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        nib_dir = args[0] 
    except:
        cookbook.doc_optparse.exit()

    nibs = {}

    for line in sys.stdin: 
        fields = line.split()
        chrom, start, end = fields[0], int( fields[1] ), int( fields[2] ) 
        if chrom in nibs:
            nib = nibs[chrom]
        else:
            nibs[chrom] = nib = bx.seq.nib.NibFile( file( "%s/%s.nib" % ( nib_dir, chrom ) ) )
        seq = nib.get( start, end - start ).upper()
        n = 0
        for i in range( 0, len( seq ) - 1 ):
            if seq[i] == 'C' and seq[i+1] == 'G':
                n += 1
        print n / len( seq )        
            

if __name__ == "__main__": __main__()
