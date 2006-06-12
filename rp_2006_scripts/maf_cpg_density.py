#!/usr/bin/env python2.4

from __future__ import division

import sys
import bx.align.maf

for m in bx.align.maf.Reader( sys.stdin ):
    if m.text_size == 0:
        print "NA"
        continue
    c = m.components[0]
    text = m.components[0].text.upper()
    cpg_count = 0
    for i in range( len( text ) - 1 ):
        if text[i] == "C" and text[i+1] == "G":
            cpg_count += 1
    print cpg_count / c.size
