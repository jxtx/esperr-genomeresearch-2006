#import weave

# Tools for dealing with multiple alignments in MAF format

class Alignment:

    def __init__( self, score=0, attributes={} ):
        self.score = 0
        self.text_size = 0
        self.attributes = {}
        self.components = []

    def add_component( self, component ):
        self.components.append( component )
        if self.text_size == 0: self.text_size = len( component.text )
        elif self.text_size != len( component.text ): raise "Components must have same text length"

    def __str__( self ):
        s = "a score=" + str( self.score )
        for key in self.attributes: 
            s += " %s=%s" % ( key, self.attributes[key] )
        s += "\n"
        # Components
        rows = []
        for c in self.components:
            rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text ) )
        s += format_tabular( rows, "llrrrrr" )
        s += "\n"
        return s

    def get_component_by_src( self, src ):
        for c in self.components:
            if c.src == src: return c
        return None

    def slice( self, start, end ):
        new = Alignment( score=self.score, attributes=self.attributes )
        for component in self.components:
            new.components.append( component.slice( start, end ) )
        new.text_size = end - start
        return new

    def slice_by_component( self, component_index, start, end ):
        ref = self.components[ component_index ]
        start_col = ref.coord_to_col( start )  
        end_col = ref.coord_to_col( end )  
        return self.slice( start_col, end_col )
    
class Component:

    def __init__( self, src='', start=0, size=0, strand=None, src_size=0, text='' ):
        self.src = src
        self.start = start
        self.size = size
        self.strand = strand
        self.src_size = src_size
        self.text = text

    def __str__( self ):
        return "s %s %d %d %s %d %s" % ( self.src, self.start, 
                                           self.size, self.strand, 
                                           self.src_size, self.text )

    def get_end( self ):
        return self.start + self.size
    end = property( fget=get_end )

    def slice( self, start, end ):
        new = Component( src=self.src, start=self.start, strand=self.strand, src_size=self.src_size )
        new.text = self.text[start:end]

        #for i in range( 0, start ):
        #    if self.text[i] != '-': new.start += 1
        #for c in new.text:
        #    if c != '-': new.size += 1
        new.start += start - self.text.count( '-', 0, start )
        new.size = len( new.text ) - new.text.count( '-' )

        return new

    def slice_by_coord( self, start, end ):
        start_col = self.coord_to_col( start )  
        end_col = self.coord_to_col( end )  
        return self.slice( start_col, end_col )
    
    def coord_to_col( self, pos ):
        if pos < self.start or pos > self.get_end():
            raise "Range error: %d not in %d-%d" % ( pos, self.start, self.get_end() )
        return self.py_coord_to_col( pos )

    def weave_coord_to_col( self, pos ):
        text = self.text
        text_size = len( self.text )
        start = self.start
        pos = pos
        return weave.inline( """
                                int col;
                                int i;
                                const char * ctext = text.c_str();
                                for ( col = 0, i = start - 1; col < text_size; ++col )
                                    if ( text[ col ] != '-' && ++i == pos ) break;
                                return_val = col;
                             """, 
                             ['text', 'text_size', 'start', 'pos' ] )

    def py_coord_to_col( self, pos ):
        if pos < self.start or pos > self.get_end():
            raise "Range error: %d not in %d-%d" % ( pos, self.start, self.get_end() )
        i = self.start
        col = 0
        text = self.text
        while i < pos:
            if text[col] != '-': i += 1
            col += 1 
        return col

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
            if len(fields) > 6: component.text = fields[6]
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
        self.file.write( str( alignment ) )

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
        
