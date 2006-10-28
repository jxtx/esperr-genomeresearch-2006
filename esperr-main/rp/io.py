import bx.align.maf
import array
import sys

from numpy import array

import rp.mapping

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
        if x is None: break
        yield x

class MappingSource( object ):
    def __init__( self, other, mapping ):
        self.other = other
        self.mapping = mapping
    def __iter__( self ):
        return self
    def next( self ):
        x = self.other.next()
        if x is None: return None
        return self.mapping.translate( x )

class IntSource( object ):
    def __init__( self, f ):
        self.f = f
    def next( self ):
        line = self.f.readline()
        if not line: return None
        return array( map( int, line.split() ), 'i' )

class MAFSource( object ):
    def __init__( self, f ):
        self.reader = bx.align.maf.Reader( f )
    def next( self ):
        maf = self.reader.next()
        if maf is None: return None
        ints = rp.mapping.DNA.translate_list( [ c.text for c in maf.components ] )
        return ints
