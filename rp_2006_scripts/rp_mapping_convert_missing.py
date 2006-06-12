#!/usr/bin/env python2.4

import sys

GAP_CHAR = "-"
MISSING_CHAR = "*"

def main():

    max_missing = int( sys.argv[1] )
    replace_missing = "--replace" in sys.argv[2:]

    keys = []
    values = {}

    for line in sys.stdin:
        fields = line.split()
        keys.append( fields[0] )
        values[ fields[0] ] = fields[1]

    for key in keys:    
        if key.count( MISSING_CHAR ) > max_missing:
            values[key] = -1
        elif replace_missing and key.count( MISSING_CHAR ) > 0:
            values[ key ] = values[ key.replace( MISSING_CHAR, GAP_CHAR ) ]

    for key in keys:
        print key, values[key]

if __name__ == "__main__":
    main()
