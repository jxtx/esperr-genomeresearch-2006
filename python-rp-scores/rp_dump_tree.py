#!/usr/bin/env python

import sys
from elementtree.ElementTree import ElementTree, Element, parse

def draw( element, level ):
    if not element: return
    symbol = element.get('symbol')
    vals = map( float, element.find( "vals" ).text.split( ',' ) )
    # Draw ascii treeness
    if level > 1:
        sys.stdout.write( " |  " * ( level - 1 ) )
    if level > 0:
        sys.stdout.write( " +--" )
    # Symbol and vals for level
    print "[" + str( symbol ) + "] (",
    for val in vals: print "%0.8f" % val,
    print "|", sum(vals),
    print ")"
    for child in element.find( "children" ):
        draw( child, level+1 )
    
if __name__ == "__main__":
    tree = parse( sys.argv[1] )
    root = tree.getroot()
    draw( root.find( "node" ), 0 )
