#!/usr/bin/env python

import seq.nib, sys

nib = seq.nib.NibFile( file( sys.argv[1] ) )
print nib.length
