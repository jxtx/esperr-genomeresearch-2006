from align import *

# Tools for dealing with multiple alignments in MAF format

class Reader:

    def __init__( self, file ):
        self.file = file
        # Read and verify maf header, store any attributes
        fields = self.file.readline().split()
        if fields[0] != '##maf': raise "File does not have MAF header"
        self.attributes = parse_attributes( fields[1:] )

    def next( self ):
        alignment = Alignment()
        # Attributes line
        line = readline( self.file, skip_blank=True )
        if not line: return None
        fields = line.split() 
        if fields[0] != 'a': raise "Expected 'a ...' line"
        alignment.attributes = parse_attributes( fields[1:] )
        alignment.score = alignment.attributes['score']
        del alignment.attributes['score']
        # Sequence lines
        while 1:
            line = readline( self.file )
            # EOF or Blank line terminates alignment components
            if not line or line.isspace(): break
            if line.isspace(): break 
            # Verify
            fields = line.split()
            if fields[0] != 's': raise "Expected 's ...' line"
            # Parse 
            component = Component()
            component.src = fields[1]
            component.start = int( fields[2] )
            component.size = int( fields[3] )
            component.strand = fields[4]
            component.src_size = int( fields[5] )
            if len(fields) > 6: component.text = fields[6].strip()
            # Add to set
            alignment.add_component( component )
        return alignment

    def __iter__( self ):
        while 1:
            v = self.next()
            if not v: break
            yield v

    def close( self ):
        self.file.close()

class Writer:

    def __init__( self, file, attributes={} ):
        self.file = file
        # Write header, Webb's maf code wants version first, we accomodate
        if not attributes.has_key('version'): attributes['version'] = 1
        self.file.write( "##maf version=%s" % attributes['version'] )
        for key in attributes: 
            if key == 'version': continue
            self.file.writelines( " %s=%s" % ( key, attributes[key] ) )
        self.file.write( "\n" )

    def write( self, alignment ):
        self.file.write( "a score=" + str( alignment.score ) )
        for key in alignment.attributes:
            self.file.write( " %s=%s" % ( key, alignment.attributes[key] ) )
        self.file.write( "\n" )
        # Components
        rows = []
        for c in alignment.components:
            rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text ) )
        self.file.write( format_tabular( rows, "llrrrrr" ) )
        self.file.write( "\n" )

    def close( self ):
        self.file.close()

# Helper methods

def readline( file, skip_blank=False ):
    """Read a line from provided file, skipping any blank or comment lines"""
    while 1:
        line = file.readline()
        if not line: return None 
        if line[0] != '#' and not ( skip_blank and line.isspace() ):
            return line

def parse_attributes( fields ):
    """Parse list of key=value strings into a dict"""
    attributes = {}
    for field in fields:
        pair = field.split( '=' )
        attributes[ pair[0] ] = pair[1]
    return attributes

def format_tabular( rows, align=None ):
    if len( rows ) == 0: return ""
    lengths = [ len( col ) for col in rows[ 0 ] ]
    for row in rows[1:]:
        for i in range( 0, len( row ) ):
            lengths[ i ] = max( lengths[ i ], len( row[ i ] ) )
    rval = ""
    for row in rows:
        for i in range( 0, len( row ) ):
            if align and align[ i ] == "l":
                rval += row[ i ].ljust( lengths[ i ] )
            else:
                rval += row[ i ].rjust( lengths[ i ] )
            rval += " "
        rval += "\n"
    return rval
        
