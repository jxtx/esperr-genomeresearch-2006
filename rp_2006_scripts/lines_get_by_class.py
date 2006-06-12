#!/usr/bin/env python

import psyco_full

import sys

import ranges, sys
import bx.align.maf
from optparse import OptionParser
from itertools import *

def __main__():

    group = sys.argv[1]

    for line, m in izip( open( sys.argv[2] ), sys.stdin ):
        if line.strip() == group:
            print m, 

if __name__ == "__main__": __main__()
