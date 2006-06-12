#!/usr/bin/env python2.3

import random, sys

alphabet_size = int( sys.argv[1] )
symbols = range( alphabet_size )

for line in sys.stdin:
    col, val = line.split()
    val = int( val )
    # If old value was 'valid', sample a new one
    if val > 0: val = random.choice( symbols )
    print col, val
