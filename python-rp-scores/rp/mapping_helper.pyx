cdef extern from "stdlib.h":
    void* malloc( int )
    void free( void* )

cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object, void **, int *) except -1
    int PyObject_AsWriteBuffer(object, void **, int *) except -1
    int PyString_AsStringAndSize(object, char **, int *) except -1

from Numeric import zeros

cdef class CharToIntArrayMapping:
    """Mapping for converting strings to int arrays"""
    
    cdef int table[256]
    cdef int out_size
    
    def __new__( self ):
        """Init empty mapping (all characters map to -1)"""
        cdef int i
        for i from 0 <= i < 256: self.table[i] = -1
        self.out_size = 0

    def set_mapping( self, char char, int symbol ):
        """Modify mapping so 'chars' map to 'symbol'"""
        self.table[ char ] = symbol
        if self.out_size <= symbol:
            self.out_size = symbol + 1

    def translate( self, string ):
        """Translate 'string' and return as int array"""
        cdef int s_len, t_len
        cdef unsigned char * s_buf
        cdef int * t_buf
        # Get direct access to string
        PyString_AsStringAndSize( string, <char **> &s_buf, &s_len )
        # Initialize empty array
        rval = zeros( s_len, 'i' )
        PyObject_AsWriteBuffer( rval, <void **> &t_buf, &t_len ) 
        # Translate
        for i from 0 <= i < s_len:
            t_buf[i] = self.table[ s_buf[ i ] ]
        # Done
        return rval
        
    def translate_list( self, strings ):
        """Translate a list of strings into an int array"""
        cdef int text_len, factor, i
        cdef int s_len, t_len
        cdef unsigned char * s_buf
        cdef int * t_buf

        # No input, no output
        if len( strings ) < 1: return None
        
        # Length of result
        text_len = len( strings[0] )
        
        # Init result array
        rval = zeros( text_len, 'i' )
        PyObject_AsWriteBuffer( rval, <void **> &t_buf, &t_len )  
              
        # Loop over seqs and accumulate result values
        factor = 1
        for string in strings:
            PyString_AsStringAndSize( string, <char **> &s_buf, &s_len )
            for i from 0 <= i < text_len:
                t_buf[i] = t_buf[i] + ( self.table[ s_buf[i] ] * factor )
            factor = factor * self.out_size
        return rval
         
cdef class IntToIntMapping:
    
    cdef int* table
    cdef int in_size
    cdef int out_size
    
    def __new__( self, int in_size ):
        self.in_size = in_size
        self.table = <int *> malloc( in_size * sizeof( int ) )
        if self.table == NULL: raise "Malloc Failed"
        for i from 0 <= i < in_size: self.table[i] = -1
        self.out_size = 0
        
    def __dealloc__( self ):
        free( self.table )

    def set_mapping( self, int index, int symbol ):
        assert ( 0 <= index < self.in_size )
        self.table[index] = symbol
        if self.out_size <= symbol:
            self.out_size = symbol + 1

    def translate( self, src ):
        """Translate 'string' and return as int array"""
        cdef int s_len, t_len
        cdef int *s_buf, *t_buf
        # Get direct access to string
        PyObject_AsReadBuffer( src, <void **> &s_buf, &s_len )
        s_len = s_len / sizeof( int )
        # Initialize empty array
        rval = zeros( s_len, 'i' )
        PyObject_AsWriteBuffer( rval, <void **> &t_buf, &t_len )
        # Translate
        for i from 0 <= i < s_len:
            t_buf[i] = self.table[ s_buf[ i ] ]
        # Done
        return rval

    def collapse( self, int a, int b ):
        cdef int i
        cdef IntToIntMapping copy
        copy = IntToIntMapping( self.in_size )
        copy.out_size = self.out_size - 1
        for i from 0 <= i < self.in_size:
            if self.table[i] == b: copy.table[i] = a
            elif self.table[i] == copy.out_size: copy.table[i] = b
            else: copy.table[i] = self.table[i]
        return copy

    def get_out_size( self ): return self.out_size
    
    def get_table( self ):
        rval = zeros( self.in_size, 'i' )
        for i in range( self.in_size ):
            rval[i] = self.table[i]
        return rval

