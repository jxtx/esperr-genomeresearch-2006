from numarray import *

import seq_numarray

DNA_BASE=5

class Mapping( object ):
    def __init__( self, f=None ):
        if f: self.read_from_file( f )

    def read_from_file( self, f ):
        align_count = None
        max_symbol = 0
        reverse_map = {}
        for line in f:
            ( key, value ) = line.split()
            if not align_count: 
                align_count = len( key )
                self.table = zeros( DNA_BASE ** align_count )
            else:
                assert align_count == len( key )
            index = seq_numarray.DNA.translate_alignment_column( key )
            self.table[ index ] = int( value )
            reverse_map.setdefault( int( value ), [] ).append( key )
            max_symbol = max( self.table[ index ], max_symbol )
        self.align_count = align_count
        self.symbol_count = max_symbol + 1
        self.reverse_table = [ reverse_map[ index ] for index in range( 0, self.symbol_count ) ]

    def translate_alignment( self, seqs ):
        return self.translate( seq_numarray.DNA.translate_alignment( seqs ) )

    def translate( self, seq ):
        return take( self.table, seq )

    def reverse( self, seq ):
        return [ self.reverse_table[ index ] for index in seq ]
