#!/usr/bin/env python2.3

from numarray import *
from cookbook.attribute import *

import string

class SeqTranslater( object ):
    attribute( table=None, base=0 )
    """Encapsulates a mapping from chars to integers"""
    def __init__( self, mapping ):
        """Initialize translation table from a list in which the 
           characters of each list item map to the same number"""
        table = zeros( 256 )
        reverse_mapping = [ '-' ]
        value = 1
        for chars in mapping:
            for ch in chars:
                put( table, array( ch, 'b' ), value )
            reverse_mapping.append( chars[ 0 ] )
            value += 1
        self.table = table 
        self.reverse_mapping = reverse_mapping
        self.base = len( mapping ) + 1
    def translate( self, seq ):
        """Convert a character sequence to a single integer array"""
        return take( self.table, array( seq, 'b' ) )
    def translate_alignment( self, seqs ):   
        """Convert the rows of a multiple alignment to a single integer array"""
        if len( seqs ) < 1: return None
        rval = zeros( len( seqs[ 0 ] ) )
        factor = 1
        for seq in seqs:
            rval += ( take( self.table, array( seq, 'b' ) ) * factor )
            factor *= self.base
        return rval

    def translate_alignment_column( self, col ):
        value = 0
        factor = 1
        for ch in col:
            index = array( ch, 'b' )[0]
            value += ( self.table[ index ] * factor )
            factor *= self.base
        return value
    def reverse_alignment_column( self, align_count, value ):
        col = []
        for i in range( 0, align_count ):
            index = ( value / self.base ** i ) % self.base 
            col.append( self.reverse_mapping[ index ] )
        return string.join( col, '' )

DNA = SeqTranslater( ( 'Aa', 'Cc', 'Gg', 'Tt' ) )
