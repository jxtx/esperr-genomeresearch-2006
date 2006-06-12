#!/usr/bin/env python2.4

import bx.align.maf
import sys

from itertools import *

from cookbook.progress_bar import ProgressBar

def get_texts( block, species ):
    texts = []
    for s in species:
        comp = block.get_component_by_src_start( s )
        if comp is None:
            assert len( texts ) > 0, "First species must always be present (%s)" % s
            texts.append( "*" * len( texts[0] ) )
        else:
            texts.append( comp.text.upper() )
    return texts
            
def update_counts( texts, counts ):
    for col in zip( * texts ):
        col = "".join( col )
        if col in counts:
            counts[ col ] += 1
        else:
            counts[ col ] = 1
            
        #try: counts[ col ] += 1
        #except: counts[ col ] = 1

counts = {}

species = sys.argv[1].split( "," )
maf_fnames = sys.argv[2:]

for fname in maf_fnames:
    print >> sys.stderr, "Processing", fname
    f = open( fname )
    # Determine file size
    f.seek( 0, 2 )
    file_size = f.tell()
    f.seek( 0, 0 )
    bar = ProgressBar( 0, file_size, 80 )
    for i, block in enumerate( bx.align.maf.Reader( f ) ):
        texts = get_texts( block, species )
        # Increment count for each column
        update_counts( texts, counts )
        if i % 100 == 0:
            bar.update_and_print( f.tell(), sys.stderr )
    print >>sys.stderr, "Done."

print "NSEQS =", len( species )
print "LENGTH = 2886607813"
print "TUPLE_SIZE = 1"
print "NTUPLES = ", len( counts )
print "NAMES = ", ",".join( species )
print "ALPHABET = ACGT"
print "NCATS = -1"
print

for i, ( col, count ) in enumerate( counts.iteritems() ):
    print i, col, count