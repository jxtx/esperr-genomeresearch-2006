#!/usr/bin/env python2.3

import sys

import ranges, sys
from align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "-m", "--maf",   action="store", help="" )

    ( options, args ) = parser.parse_args()

    maf_reader = maf.Reader( file( options.maf ) )

    for m in maf_reader: print m.score

if __name__ == "__main__": __main__()
