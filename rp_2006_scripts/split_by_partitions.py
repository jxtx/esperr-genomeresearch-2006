#!/usr/bin/env python2.4

from __future__ import division

import sys

import pkg_resources; pkg_resources.require("bx-python")
from bx.bitset import *
from bx.bitset_builders import *

# Load up the partitions
partitions = list()
for line in open( sys.argv[1] ):
    name, score, filename = line.split()
    partitions.append( ( name, score, binned_bitsets_from_file( open( filename ) ) ) )


for line in sys.stdin:
    fields = line.rstrip().split( )
    chr, start, end = fields[0], int( fields[1] ), int( fields[2] )
    if len( fields ) > 3: 
        label = fields[3]
    else: 
        label = ""
    if len( fields ) > 4:
        strand = fields[4]
    else:
        strand = "+"
    # Find which partition it overlaps
    for name, score, bb in partitions:
        # Is there at least 1bp overlap?
        if chr in bb:
            overlap = bb[chr].count_range( start, end-start )
            if overlap > 0:
                break
    else:
        # No overlap with any partition? For now throw this since the 
        # partitions tile the encode regions completely, indicate an interval
        # that does not even overlap an encode region
        # print >> sys.stderr, "warning: Interval (%s, %d, %d) does not overlap any partition" % ( chr, start, end )
        name = "OTHER"
        score = -1
        overlap = (end-start)
    # Annotate with the name of the partition
    frac_overlap = overlap / (end-start)
    # BED6 plus?
    print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%0.4f" % ( chr, start, end, label, score, strand, name, frac_overlap )
