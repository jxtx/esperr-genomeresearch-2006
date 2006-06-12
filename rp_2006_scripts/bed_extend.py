#!/usr/bin/env python2.4

import sys

padl = int( sys.argv[1] )
padr = int( sys.argv[2] )

for line in sys.stdin:
    fields = line.rstrip( "\r\n" ).split("\t")
    fields[1] = str( max( 0, int( fields[1] ) - padl ) )
    fields[2] = str( int( fields[2] ) + padr )
    print "\t".join( fields )
    
