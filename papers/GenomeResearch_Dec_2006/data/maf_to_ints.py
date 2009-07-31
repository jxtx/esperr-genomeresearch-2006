#!/usr/bin/env python

import sys
import bx.align.maf
import rp.mapping 

size, mapping = rp.mapping.alignment_mapping_from_file( file( sys.argv[1] ) )

for maf in bx.align.maf.Reader( sys.stdin ):
    ints = rp.mapping.DNA.translate_list( [ c.text for c in maf.components ] )
    ints = mapping.translate( ints )
    print ' '.join( map( str, ints ) )
