from __future__ import division

import sys, struct, string, math

NIB_MAGIC_NUMBER = 0x6BE93D3A
NIB_HEADER_SIZE = 8
NIB_I2C_TABLE = "TCAGNXXXtcagnxxx"

class NibFile( object ):
    def __init__( self, file ):
        ( magic, length ) = struct.unpack( "LL", file.read( NIB_HEADER_SIZE ) )
        if magic != NIB_MAGIC_NUMBER: raise "Not a NIB file"
        self.file = file
        self.magic = magic
        self.length = length
    def get( self, start, len ):
        assert start >= 0
        assert start + len - 1 < self.length
        # Read block of bytes containing sequence
        block_start = int( math.floor( start / 2 ) )
        block_len = int( math.ceil( len / 2 ) )
        self.file.seek( NIB_HEADER_SIZE + block_start )
        data = struct.unpack( "%dB" % block_len, self.file.read( block_len ) )
        # Translate to character representation
        result = []
        for value in data:
            result.append( NIB_I2C_TABLE[ ( value >> 4 ) & 0xF ] )
            result.append( NIB_I2C_TABLE[ ( value >> 0 ) & 0xF ] )
        # Trim if start / end are odd 
        if start & 1: del result[ 0 ]
        if len & 1: del result[ -1 ]
        # Return as string
        return string.join( result, '' )
