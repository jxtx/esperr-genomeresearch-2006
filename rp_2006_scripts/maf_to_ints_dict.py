#!/usr/bin/env python2.3

"""
For each block in a maf file (read from stdin) write a sequence of ints 
corresponding to the columns of the block after applying the provided mapping.

The 'correct' number of species is determined by the mapping file, blocks not having
this number of species will be ignored.

usage: %prog mapping_file
"""

from __future__ import division

import psyco_full

import bx.align.maf
import string
import sys

def main():

    mapping = dict()
    for line in open( sys.argv[1] ):
        fields = line.split()
        mapping[fields[0]] = fields[1]

    for maf in bx.align.maf.Reader( sys.stdin ):
        for col in maf.column_iter():
            key = "".join( col ).upper()
            try:
                print mapping[ key ],
            except KeyError:
                print -1,
        print

if __name__ == "__main__": main()
