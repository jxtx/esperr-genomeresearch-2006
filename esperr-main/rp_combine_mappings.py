#!/usr/bin/env python

import sys

def main():
    # Load mappings to be merged
    sets_1, forward_1 = load( open( sys.argv[1] ) )
    sets_2, forward_2 = load( open( sys.argv[2] ) )
    # Split each set in the first according to the second
    all_new_sets = []
    for set in sets_1.values():
        new_sets = dict()
        for item in set:
            try: new_sets[ forward_2[ item ] ].append( item )
            except: new_sets[ forward_2[ item ] ] = [ item ]
        all_new_sets += new_sets.values()
    print >>sys.stderr, len( all_new_sets )
    new_forward = dict()
    for i, set in enumerate( all_new_sets ):
        for col in set:
            new_forward[col] = i
    for line in open( sys.argv[1] ):
        col, val = line.split()
        val = int(val)
        if val < 0:
            print col, val
        else:
            print col, new_forward[col] 

def load( f ):
    sets = dict()
    forward = dict()
    for line in f:
        col, val = line.split()
        val = int(val)
        if val >= 0:
            try: sets[val].append( col )
            except: sets[val] = [ col ]
            forward[col] = val
    return sets, forward

if __name__ == "__main__": main()

