#!/usr/bin/env python2.3

import sys

import ranges, sys, random
from align import maf
from math import *
from optparse import OptionParser

def __main__():

    if len( sys.argv ) > 1: fraction = float( sys.argv[1] )

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    mafs = list( maf_reader )

    # for m in maf_reader: mafs.append( m )

    random.shuffle( mafs )

    count = len( mafs )
    if fraction: count = int( floor( count * fraction ) )

    for i in range( 0, count ): maf_writer.write( mafs[ i ] )
    
if __name__ == "__main__": __main__()
