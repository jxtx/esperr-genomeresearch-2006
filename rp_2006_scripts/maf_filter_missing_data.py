#!/usr/bin/env python2.4

from __future__ import division

import psyco_full

import sys

import ranges, sys
from bx.align import maf
from optparse import OptionParser

MIN_PER_GROUP = 1
MIN_REAL = .50

def frac_real_bases( text ):
    u = text.upper()
    return ( u.count( 'A' ) + u.count( 'C' ) + u.count( 'T' ) + u.count( 'G' ) ) / len( text )

def main():

    # Groups is a list like "hg17,panTro1;mm3,rn3". A valid block will have 
    # sequence in at least 1 component from each group
    groups = map( lambda s: s.split( "," ), sys.argv[1].split( ";" ) )

    if len( sys.argv ) > 2:
        min_length = int( sys.argv[2] )
    else:    
        min_length = 100

    ## print >> sys.stderr, groups

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:
        good = True
        for group in groups:
            number_in_group = 0
            for src in group:
                c = m.get_component_by_src_start( src )
                ## print >>sys.stderr, src, c.text, frac_real_bases( c.text )
                if c and frac_real_bases( c.text ) >= MIN_REAL:
                    number_in_group += 1
            if number_in_group < MIN_PER_GROUP:
                good = False
                break
        if good and len( m.components[0].text ) >= min_length:
            maf_writer.write( m )

if __name__ == "__main__": 
        main()
