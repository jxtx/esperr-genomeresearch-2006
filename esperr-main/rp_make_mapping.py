#!/usr/bin/env python2.3

"""
Make an identity mapping for some number of species 
(each column mapped to unique symbol)
"""

import sys
import rp.mapping

align_count = int( sys.argv[1] )

mapping = rp.mapping.identity_mapping( rp.mapping.DNA.get_out_size() ** align_count )

for i, symbol in enumerate( mapping.get_table() ):
    print  str.join( '', rp.mapping.DNA.reverse_map( i, align_count ) ), symbol
