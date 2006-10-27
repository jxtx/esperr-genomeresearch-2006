#!/usr/bin/env python2.4

import sys
for line in sys.stdin:
    if line.startswith( "variable" ):
        print line,
        continue
    fields = line.split()
    if fields[1][0] == "-":
        print fields[0] + "\t" + "0.0"
    else:
        print line,
