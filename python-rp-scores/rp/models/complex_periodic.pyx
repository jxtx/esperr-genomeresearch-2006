cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, void **buffer, int *buffer_len) except -1
    int PyObject_AsWriteBuffer(object obj, void **buffer, int *buffer_len) except -1

cdef extern from "stdlib.h":    
    void* malloc( int )
    void free( void* )

cdef extern from "simple_periodic_core.h":
    ctypedef double real
    ctypedef int bool
    int matrix_size( int order, int radix )
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
    bool score_string_2( int order1, int order2, int radix1, int radix2, 
                     real** score_matrix1, real *** score_matrix2,
                     int* string1, int* string2, 
                     int start, int length, real* rval, int period )
    void fill_in_counts_2( int order1, int order2, int radix1, int radix2, 
                       int** counts, int* string1, int* string2,
                       int string_len, int phase, int period )
    int** new_counts_2( int order1, int order2, int radix1, int radix2 )

import sys
from array import array

cdef class Model:
    cdef int order1
    cdef int order2
    cdef int radix1
    cdef int radix2
    cdef real** scores1
    cdef real*** scores2
    cdef int nmodels
    cdef int scores_len1
    cdef int scores_len2
    cdef int period

    cdef init( self, int order1, int order2, int radix1, int radix2, real** scores1, real*** scores2, int period ):
        self.order1 = order1
        self.order2 = order2
        self.radix1 = radix1
        self.radix2 = radix2
        self.scores1 = scores1
        self.scores2 = scores2
        self.nmodels = matrix_size( order2, radix1 )
        self.scores_len1 = radix1 ** ( order1 + 1 )
        self.scores_len2 = radix2 ** ( order2 + 1 )
        self.period = period

    # def get_order( self ):
    #     return self.order
    
    # def get_radix( self ):
    #     return self.radix
    
    # def get_scores( self ):
    #     a = array( 'd' )
    #     for i from 0 <= i < self.scores_len:
    #         a.append( self.scores[i] )
    #     return a

    def score( self, string, int start=0, int length=-1 ):
        cdef int* buf1
        cdef int buf_len1
        cdef int* buf2
        cdef int buf_len2
        cdef real rval
        PyObject_AsReadBuffer( string[0], <void**> &buf1, &buf_len1 )
        PyObject_AsReadBuffer( string[1], <void**> &buf2, &buf_len2 )
        buf_len1 = buf_len1 / sizeof( int )
        buf_len2 = buf_len2 / sizeof( int )
        assert buf_len1 == buf_len2
        if length < 0:
            length = buf_len1
        else:
            assert start + length <= buf_len1
        if score_string_2( self.order1, self.order2, self.radix1, self.radix2, 
                           self.scores1, self.scores2, buf1, buf2, start, length, &rval, self.period ):
            return rval
        else:
            return None

    # def score_positions( self, string, target ):
    #     assert string.typecode() == "i", "String must be int array"
    #     assert target.typecode() == "f", "Target must be float array"
    #     cdef int* buf
    #     cdef float* t_buf
    #     cdef int buf_len, t_buf_len
    #     PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
    #     buf_len = buf_len / sizeof( int )
    #     PyObject_AsWriteBuffer( target, <void**> &t_buf, &t_buf_len )
    #     t_buf_len = t_buf_len / sizeof( float )
    #     assert buf_len == t_buf_len, "String and target should have same size"        
    #     score_string_positions( self.order, self.radix, self.scores, buf, t_buf, 0, buf_len )

    # def to_file( self, file ):
    #     file.write( "order: " + str( self.order ) + "\n" )
    #     file.write( "radix: " + str( self.radix ) + "\n" )
    #     file.write( "period: " + str( self.period ) + "\n" )
    #     for phase in range( self.period ):
    #         file.write( "model_for_phase: " + str( phase ) + "\n" )
    #         for i from 0 <= i < self.scores_len:
    #             file.write( str( self.scores[phase][i] ) + "\n" )

    def __dealloc__( self ):
        for i from 0 <= i < self.period:
            free_scores( self.scores1[i] )
            for j from 0 <= j < self.nmodels:
                free_scores( self.scores2[i][j] )
            free( self.scores2[i] )
        free( self.scores1 )
        free( self.scores2 )

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

def train( int order1, int order2, int radix1, int radix2, pos_strings, neg_strings, **kwargs ):
    
    cdef int* pos_counts1
    cdef int* neg_counts1

    cdef int** pos_counts2
    cdef int** neg_counts2

    cdef real* pos_probs1
    cdef real* neg_probs1
    
    cdef real* pos_probs2
    cdef real* neg_probs2

    cdef real* scores
    cdef real** scores2

    cdef real** scores_by_period1
    cdef real*** scores_by_period2

    cdef Model rval
    
    cdef int phase
    cdef int period
    
    cdef int nmodels
    
    # Number of different models for the second sequence, based on patterns 
    # in the first sequence
    nmodels = matrix_size( order2, radix1 )
    
    pos_period = int( kwargs.get( 'pos_period', 1 ) )

    scores_by_period1 = <real **> malloc( pos_period * sizeof( real* ) )
    scores_by_period2 = <real ***> malloc( pos_period * sizeof( real** ) )
    
    for phase from 0 <= phase < pos_period:
        
        # Positive probs for model 1
        pos_counts1 = new_counts( order1, radix1 )
        if pos_counts1 == NULL: raise "Malloc failed (counts)"
        fill_in( order1, radix1, pos_counts1, pos_strings, phase, pos_period )
        pos_probs1 = counts_to_probs( order1, radix1, pos_counts1 )
        if pos_probs1 == NULL: raise "Malloc failed (probs)"
        free_counts( pos_counts1, order1 )
        
        # Negative probs for model 1
        neg_counts1 = new_counts( order1, radix1 )
        if neg_counts1 == NULL: raise "Malloc failed (neg counts)"
        fill_in( order1, radix1, neg_counts1, neg_strings, 0, 1 )
        neg_probs1 = counts_to_probs( order1, radix1, neg_counts1 )
        if neg_probs1 == NULL: raise "Malloc failed (probs)"
        free_counts( neg_counts1, order1 )
        
        # Scores for model 1
        scores_by_period1[phase] = probs_to_score_matrix( order1, radix1, 
                                                          pos_probs1, neg_probs1 )
        if scores_by_period1[phase] == NULL: raise "Malloc failed (scores)"
        free_probs( pos_probs1, order1 )
        free_probs( neg_probs1, order1 )

        # Positive and negative counts for model 2
        pos_counts2 = new_counts_2( order1, order2, radix1, radix2 )
        if pos_counts2 == NULL: raise "Malloc failed (counts)"
        fill_in_2( order1, order2, radix1, radix2, pos_counts2, pos_strings, phase, pos_period )
        neg_counts2 = new_counts_2( order1, order2, radix1, radix2 )
        if neg_counts2 == NULL: raise "Malloc failed (counts)"
        fill_in_2( order1, order2, radix1, radix2, neg_counts2, neg_strings, 0, 1 )
        
        # Probs and scores for each sub model 2
        scores_by_period2[phase] = <real **> malloc( nmodels * sizeof( real* ) )
        for i from 0 <= i < nmodels:
            
            pos_probs2 = counts_to_probs( order2, radix2, pos_counts2[i] )
            if pos_probs2 == NULL: raise "Malloc failed (probs)"
            neg_probs2 = counts_to_probs( order2, radix2, neg_counts2[i] )
            if neg_probs2 == NULL: raise "Malloc failed (probs)"
                
            scores_by_period2[phase][i] = probs_to_score_matrix( order2, radix2, 
                                                                 pos_probs2, neg_probs2 )
            if scores_by_period2[phase][i] == NULL: raise "Malloc failed (scores)"
            free_probs( pos_probs2, order2 )
            free_probs( neg_probs2, order2 )

    rval = Model()
    rval.init( order1, order2, radix1, radix2, 
               scores_by_period1, scores_by_period2, pos_period )

    return rval

cdef fill_in( int order, int radix, int* counts, strings, int phase, int period ):
    cdef int* buf
    cdef int buf_len
    for string, _ in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        fill_in_counts( order, radix, counts, buf, buf_len / sizeof( int ), phase, period )

cdef fill_in_2( int order1, int order2, int radix1, int radix2, 
                int** counts, strings, int phase, int period ):
    cdef int* buf1
    cdef int* buf2
    cdef int buf_len1
    cdef int buf_len2
    for string1, string2 in strings:
        PyObject_AsReadBuffer( string1, <void**> &buf1, &buf_len1 )
        PyObject_AsReadBuffer( string2, <void**> &buf2, &buf_len2 )
        fill_in_counts_2( order1, order2, radix1, radix2, counts, buf1, buf2, buf_len1 / sizeof( int ), phase, period )



