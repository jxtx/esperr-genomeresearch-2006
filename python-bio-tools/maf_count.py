#!/usr/bin/env python2.3

"""
Read a MAF from standard input and count alignments, bases, or columns. 
"""

import sys

from align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "-c", "--cols",  action="store_const", dest="action", const="cols" )
    parser.add_option( "-b", "--bases", action="store_const", dest="action", const="bases" )
    parser.add_option( "-e", "--each", action="store_true", dest="each" )

    ( options, args ) = parser.parse_args()

    action = "aligns"
    if options.action: action = options.action
    
    print_each = bool( options.each )

    maf_reader = maf.Reader( sys.stdin )

    count = 0

    for m in maf_reader:
        
        if action == "aligns": count += 1
        elif action == "cols": count += m.text_size
        elif action == "bases": count += m.components[0].size

        if print_each: 
            print count
            count = 0

    if not print_each: print count

if __name__ == "__main__": __main__()
