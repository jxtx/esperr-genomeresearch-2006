from mapping_helper import *

# Char->Int mapping for DNA characters
                
DNA_BASE=5
                
DNA = CharToIntArrayMapping()
DNA.set_mapping( "a", 0 )
DNA.set_mapping( "A", 0 )
DNA.set_mapping( "c", 1 )
DNA.set_mapping( "C", 1 )
DNA.set_mapping( "g", 2 )
DNA.set_mapping( "G", 2 )
DNA.set_mapping( "t", 3 )
DNA.set_mapping( "T", 3 )
DNA.set_mapping( "N", 4 )
DNA.set_mapping( "n", 4 )
DNA.set_mapping( "-", 4 )

# Creating mappings

def alignment_mapping_from_file( f ):
    """Create a mapping from a file of alignment columns"""
        
    columns, symbols = [], []
    for line in f:
        column, symbol = line.split()
        columns.append( column )
        symbols.append( int( symbol ) )
                
    align_count = len( columns[0] )
        
    mapping = IntToIntMapping( DNA_BASE ** align_count )
        
    for column, symbol in zip( columns, symbols ):
        index = DNA.translate_list( list( column ) )[0]
        mapping.set_mapping( index, symbol )

    return align_count, mapping

def identity_mapping( size ):
    mapping = IntToIntMapping( size )
    for i in range( size-1 ):
        mapping.set_mapping( i, i )
    return mapping    
    
