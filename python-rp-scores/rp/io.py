import align.maf
import array
import alphabet
import seq_numarray
import sys

def get_reader( f, format=None, mapping=None ):
    # Open correct reader for format
    if format == 'maf':
        source = MAFSource( f )
    else:
        source = IntSource( f )
    # Wrap in mapping if needed
    if mapping:
        source = MappingSource( source, mapping )
    # Force result to return int arrays
    while 1:
        x = source.next()
        if not x: break
        yield array.array( 'i', list( x ) )

class Reader( object ):
    def __init__( self, source ):
        self.source = source
    def __iter__( self ):
        return self
    def next( self ):
        x = self.source.next()
        if x is None: raise StopIteration
        return array.array( 'i', list( self.source.next() ) )

class MappingSource( object ):
    def __init__( self, other, mapping ):
        self.other = other
        self.mapping = mapping
    def next( self ):
        x = self.other.next()
        if x is None: return None
        return self.mapping.translate( self.other.next() )

class IntSource( object ):
    def __init__( self, f ):
        self.f = f
    def next( self ):
        line = self.f.readline()
        if not line: return None
        return map( int, line.split() )

class MAFSource( object ):
    def __init__( self, f ):
        self.reader = align.maf.Reader( f )
    def next( self ):
        maf = self.reader.next()
        if maf is None: return None
        ints = seq_numarray.DNA.translate_alignment( [ c.text for c in maf.components ] )
        return ints
