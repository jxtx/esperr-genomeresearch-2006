#!/usr/bin/env python2.3

from numarray import *
from cookbook.attribute import *

class SeqTranslater( object ):
    attribute( table=None, base=0 )
    """Encapsulates a mapping from chars to integers"""
    def __init__( self, mapping ):
        """Initialize translation table from a list in which the 
           characters of each list item map to the same number"""
        table = zeros( 256 )
        value = 1
        for chars in mapping:
            for ch in chars:
                put( table, array( ch, 'b' ), value )
            value += 1
        self.table = table 
        self.base = len( mapping ) + 1
    def translate( self, seq ):
        """Convert a character sequence to a single integer array"""
        return take( self.table, array( seq, 'b' ) )
    def translate_alignment_column( self, col ):
        value = 0
        factor = 1
        for ch in col:
            index = array( ch, 'b' )[0]
            value += ( self.table[ index ] * factor )
            factor *= self.base
        return value

DNA = SeqTranslater( ( 'Aa', 'Cc', 'Gg', 'Tt' ) )
