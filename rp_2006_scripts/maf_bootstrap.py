#!/usr/bin/env python2.3

import sys

import ranges, sys, random
from bx.align import maf
from math import *
from optparse import OptionParser

def main():

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    mafs = list( maf_reader )

    for i in range( 0, len( mafs ) ):
        maf_writer.write( random.choice( mafs ) )
    
if __name__ == "__main__": 
    main()
