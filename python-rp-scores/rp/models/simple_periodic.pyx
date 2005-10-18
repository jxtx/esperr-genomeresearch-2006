cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, void **buffer, int *buffer_len) except -1
    int PyObject_AsWriteBuffer(object obj, void **buffer, int *buffer_len) except -1

cdef extern from "stdlib.h":    
    void* malloc( int )
    void free( void* )

cdef extern from "simple_periodic_core.h":
    ctypedef double real
    ctypedef int bool
    int* new_counts( int order, int radix )
    void free_counts( int* v, int order )
    void fill_in_counts( int order, int radix, int* counts, int* string, int string_len, int phase, int period )
    real* counts_to_probs( int order, int radix, int* counts )
    void free_probs( real* probs, int order )
    real* new_real_array( int size )
    real* probs_to_score_matrix( int order, int radix, real* pos_probs, real* neg_probs )
    void free_scores( real* scores )
    bool score_string( int order, int radix, real** score_matrix, int* string, int start, int length, real* rval, int period )
    # bool score_string_positions( int order, int radix, real* score_matrix, int* string, float* target, int start, int length )
    

import sys
from array import array

cdef class Model:
    cdef int order
    cdef int radix
    cdef real** scores
    cdef int scores_len
    cdef int period

    cdef init( self, int order, int radix, real** scores, int period ):
        self.order = order
        self.radix = radix
        self.scores = scores
        self.scores_len = radix ** ( order + 1 )
        self.period = period

    def get_order( self ):
        return self.order

    def get_radix( self ):
        return self.radix

#     def get_scores( self ):
#         a = array( 'd' )
#         for i from 0 <= i < self.scores_len:
#             a.append( self.scores[i] )
#         return a

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
        if score_string( self.order, self.radix, self.scores, buf, start, length, &rval, self.period ):
            return rval
        else:
            return None

#     def score_positions( self, string, target ):
#         assert string.typecode() == "i", "String must be int array"
#         assert target.typecode() == "f", "Target must be float array"
#         cdef int* buf
#         cdef float* t_buf
#         cdef int buf_len, t_buf_len
#         PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
#         buf_len = buf_len / sizeof( int )
#         PyObject_AsWriteBuffer( target, <void**> &t_buf, &t_buf_len )
#         t_buf_len = t_buf_len / sizeof( float )
#         assert buf_len == t_buf_len, "String and target should have same size"        
#         score_string_positions( self.order, self.radix, self.scores, buf, t_buf, 0, buf_len )

    def to_file( self, file ):
        file.write( "order: " + str( self.order ) + "\n" )
        file.write( "radix: " + str( self.radix ) + "\n" )
        file.write( "period: " + str( self.period ) + "\n" )
        for phase in range( self.period ):
            file.write( "model_for_phase: " + str( phase ) + "\n" )
            for i from 0 <= i < self.scores_len:
                file.write( str( self.scores[phase][i] ) + "\n" )

    def __dealloc__( self ):
        for i from 0 <= i < self.period:
            free_scores( self.scores[i] )
        free( self.scores )

# def from_file( f ):
#     cdef int order
#     cdef int size
#     cdef int radix
#     cdef real* scores
#     cdef Model rval
#     # Parse order
#     try:
#         fields = f.next().split()
#         assert fields[0] == 'order:'
#         order = int( fields[1] )
#     except:
#         raise "Expected 'order:'"
#     # Parse radix
#     try:
#         fields = f.next().split()
#         assert fields[0] == 'radix:'
#         radix = int( fields[1] )
#     except:
#         raise "Expected 'radix:'"
#     # Parse values
#     size = radix ** ( order + 1 )
#     scores = new_real_array( size )
#     for i from 0 <= i < size:
#         scores[i] = float( f.next() )
#     # Build result
#     rval = StandardModel()
#     rval.init( order, radix, scores )
#     return rval

def train( int order, int radix, pos_strings, neg_strings, **kwargs ):
    
    cdef int* pos_counts
    cdef int* neg_counts

    cdef real* pos_probs
    cdef real* neg_probs

    cdef real* scores

    cdef real** scores_by_period

    cdef Model rval
    
    cdef int phase
    cdef int period
    
    pos_period = int( kwargs.get( 'pos_period', 1 ) )

    scores_by_period = <real **> malloc( pos_period * sizeof( real* ) )

    ### Negative model is NOT periodic here! we can reuse it. 
    ### No wait we can't, it breaks the output which means something is wrong. Ick. 
    # neg_counts = new_counts( order, radix )
    # if neg_counts == NULL: raise "Malloc failed (neg counts)"
    # fill_in( order, radix, neg_counts, neg_strings, 0, 1 )
    # neg_probs = counts_to_probs( order, radix, neg_counts )
    # if neg_probs == NULL: raise "Malloc failed (probs)"
    # free_counts( neg_counts, order )

    for phase in range( pos_period ):
        
        neg_counts = new_counts( order, radix )
        if neg_counts == NULL: raise "Malloc failed (neg counts)"
        fill_in( order, radix, neg_counts, neg_strings, 0, 1 )
        neg_probs = counts_to_probs( order, radix, neg_counts )
        if neg_probs == NULL: raise "Malloc failed (probs)"
        free_counts( neg_counts, order )

        pos_counts = new_counts( order, radix )
        if pos_counts == NULL: raise "Malloc failed (counts)"
        fill_in( order, radix, pos_counts, pos_strings, phase, pos_period )
        pos_probs = counts_to_probs( order, radix, pos_counts )
        if pos_probs == NULL: raise "Malloc failed (probs)"
        free_counts( pos_counts, order )
        scores_by_period[phase] = probs_to_score_matrix( order, radix, 
                                                         pos_probs, neg_probs )
        if scores_by_period[phase] == NULL: raise "Malloc failed (scores)"
        free_probs( pos_probs, order )
    
    # free_probs( neg_probs, order )

    rval = Model()
    rval.init( order, radix, scores_by_period, pos_period )

    return rval

cdef fill_in( int order, int radix, int* counts, strings, int phase, int period ):
    cdef int* buf
    cdef int buf_len
    for string in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        fill_in_counts( order, radix, counts, buf, buf_len / sizeof( int ), phase, period )
