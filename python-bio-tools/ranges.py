"""Some classes for working with ranges. This needs to be reconciled with 'intervals'"""

class Range:
    def __init__( self, start, end ):
        self.start = start
        self.end = end
    def __repr__( self ):
        return "Range( %d, %d )" % ( self.start, self.end )
    def __cmp__( self, other ):
        return cmp( self.start, other.start ) or cmp( self.end, other.end )

class Reader:
    def __init__( self, file, indexes=(0,1), sep=None ):
        self.file = file
        self.indexes = indexes
        self.sep = sep
    def next( self ):
        line = self.file.readline()
        if not line: return None
        fields = line.split( self.sep )
        return Range( int( fields[ self.indexes[0] ] ), int( fields[ self.indexes[1] ] ) )
    def close( self ):
        self.file.close()
    def __iter__( self ):
        while 1:
            v = self.next()
            if not v: break 
            yield v

def merge_ranges( range_set ):

    ranges = iter( range_set )
    rval = []

    first = ranges.next()

    start = first.start
    end = first.end

    for range in ranges:
        # Ensure sorted
        if range.start < start: raise "Input must be sorted"
        # If this range overlaps a previous one, just extend
        if range.start <= end:
            end = max( end, range.end )
        # Otherwise, right out the range
        else:
            rval.append( Range( start, end ) )
            start = range.start
            end = range.end

    if start > 0: rval.append( Range( start, end ) )

    return rval

def invert_ranges( ranges, start, end ):
    rval = []
    last = start
    for range in ranges:
        if range.start > last:
            rval.append( Range( last, range.start ) )
        last = range.end
    if last < end:
        rval.append( Range( last, end ) )
    return rval



