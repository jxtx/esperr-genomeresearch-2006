cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, void **buffer, int *buffer_len) except -1
    int PyObject_AsWriteBuffer(object obj, void **buffer, int *buffer_len) except -1

cdef extern from "standard_core.h":
    ctypedef double real
    ctypedef int bool
    int** new_counts( int order, int radix )
    void free_counts( int** v, int order )
    void fill_in_counts( int order, int radix, int** counts, int* string, int string_len )
    real** counts_to_probs( int order, int radix, int** counts, bool average, bool backoff )
    void free_probs( real** probs, int order )
    real* new_real_array( int size )
    real* probs_to_score_matrix( int order, int radix, real** pos_probs, real** neg_probs, bool averaging )
    void free_scores( real* scores )
    bool score_string( int order, int radix, real* score_matrix, int* string, int start, int length, real* rval )
    bool score_string_positions( int order, int radix, real* score_matrix, int* string, float* target, int start, int length )
    

import sys
from array import array

cdef class StandardModel:
    cdef int order
    cdef int radix
    cdef real* scores
    cdef int scores_len

    cdef init( self, int order, int radix, real* scores ):
        self.order = order
        self.radix = radix
        self.scores = scores
        self.scores_len = <int> ( ( <float> radix ) ** ( order + 1 ) )

    def get_order( self ):
        return self.order

    def get_radix( self ):
        return self.radix

    def get_scores( self ):
        a = array( 'd' )
        for i from 0 <= i < self.scores_len:
            a.append( self.scores[i] )
        return a

    def score( self, string, int start=0, int length=-1 ):
        cdef int* buf
        cdef int buf_len
        cdef real rval
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        buf_len = buf_len / sizeof( int )
        if length < 0:
            length = buf_len
        else:
            assert start + length <= buf_len
        if score_string( self.order, self.radix, self.scores, buf, start, length, &rval ):
            return rval
        else:
            return None

    def score_positions( self, string, target ):
        assert string.typecode() == "i", "String must be int array"
        assert target.typecode() == "f", "Target must be float array"
        cdef int* buf
        cdef float* t_buf
        cdef int buf_len, t_buf_len
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        buf_len = buf_len / sizeof( int )
        PyObject_AsWriteBuffer( target, <void**> &t_buf, &t_buf_len )
        t_buf_len = t_buf_len / sizeof( float )
        assert buf_len == t_buf_len, "String and target should have same size"        
        score_string_positions( self.order, self.radix, self.scores, buf, t_buf, 0, buf_len )

    def to_file( self, file ):
        file.write( "order: " + str( self.order ) + "\n" )
        file.write( "radix: " + str( self.radix ) + "\n" )
        for i from 0 <= i < self.scores_len:
            file.write( str( self.scores[i] ) + "\n" )

    def __dealloc__( self ):
        free_scores( self.scores )

def from_file( f ):
    cdef int order
    cdef int size
    cdef int radix
    cdef real* scores
    cdef StandardModel rval
    # Parse order
    try:
        fields = f.next().split()
        assert fields[0] == 'order:'
        order = int( fields[1] )
    except:
        raise "Expected 'order:'"
    # Parse radix
    try:
        fields = f.next().split()
        assert fields[0] == 'radix:'
        radix = int( fields[1] )
    except:
        raise "Expected 'radix:'"
    # Parse values
    size = <int> ( ( <float> radix ) ** ( order + 1 ) )
    scores = new_real_array( size )
    for i from 0 <= i < size:
        scores[i] = float( f.next() )
    # Build result
    rval = StandardModel()
    rval.init( order, radix, scores )
    return rval

def train( int order, int radix, pos_strings, neg_strings, **kwargs ):
    
    cdef int** pos_counts
    cdef int** neg_counts

    cdef real** pos_probs
    cdef real** neg_probs

    cdef real* scores

    cdef StandardModel rval

    averaging = int( kwargs.get( 'averaging', 0 ) )
    backoff = int( kwargs.get( 'backoff', 1 ) )

    pos_counts = new_counts( order, radix )
    neg_counts = new_counts( order, radix )
    
    if pos_counts == NULL or neg_counts == NULL: 
        raise "Malloc failed (counts)"

    fill_in( order, radix, pos_counts, pos_strings )
    fill_in( order, radix, neg_counts, neg_strings )
    
    pos_probs = counts_to_probs( order, radix, pos_counts, averaging, backoff )
    neg_probs = counts_to_probs( order, radix, neg_counts, averaging, backoff )

    if pos_counts == NULL or neg_counts == NULL: 
        raise "Malloc failed (probs)"

    free_counts( pos_counts, order )
    free_counts( neg_counts, order )

    scores = probs_to_score_matrix( order, radix, 
                                    pos_probs, neg_probs, averaging )

    if scores == NULL: 
        raise "Malloc failed (scores)"

    free_probs( pos_probs, order )
    free_probs( neg_probs, order )

    rval = StandardModel()
    rval.init( order, radix, scores )

    return rval

cdef fill_in( int order, int radix, int** counts, strings ):
    cdef int* buf
    cdef int buf_len
    for string in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        fill_in_counts( order, radix, counts, buf, buf_len / sizeof( int ) )
