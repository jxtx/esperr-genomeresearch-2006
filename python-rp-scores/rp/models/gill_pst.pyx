import sys

cdef extern from "stdio.h":
    ctypedef struct FILE
    int fprintf( FILE*, char*, ... )
    int fscanf( FILE*, char*, ... )

cdef extern from "stdlib.h":    
    void* malloc( int )
    void free( void* )
    
cdef extern from "string.h":
    char* strdup( char* src )

cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object, void **, int *) except -1
    int PyObject_AsWriteBuffer(object, void **, int *) except -1
    int PyString_AsStringAndSize(object, char **, int *) except -1
    FILE* PyFile_AsFile( object ) except? NULL

cdef extern from "pst.h":
    cdef struct pst_struct:
        char smooth
        double *p_c
        pst_struct **sonsp
    ctypedef pst_struct *pst_type
    double log10like_on_pst( char* alphabet, pst_type tree, char* string )
    void add_probs_pst( char *AB, int absize, int max_string, int m, char **ds_ptr, pst_type T )
    void smoothen_pst( int absize, double min_prob, pst_type T )

cdef extern from "learn_algs.h":
    pst_type grow_tree1 (char *AB, int absize, int m, char **ts_ptr, int *li, int *pot_nodesp, pst_type T0)
    
# NOTE: Parameters for learning are GLOBAL VARIABLES which we need to set
#       before calling the above. Consequently training this model is 
#       NOT THREAD SAFE!!
cdef extern double p_min
cdef extern double alpha
cdef extern double gamma_min
cdef extern double p_ratio
cdef extern int L_max

cdef extern from "file_handle.h":
    void rec_save_pst ( FILE *fp, int absize, pst_type T )
    void rec_load_pst ( FILE *fp, int absize, pst_type *Tp )

# Character to start at when encoding int arrays as strings
cdef int offset
offset = 33

cdef class PSTModel:
    cdef int radix
    cdef char* alphabet
    cdef int l_max
    cdef pst_type pos_tree, neg_tree

    cdef init( self, int radix, char* alphabet, int l_max, pst_type pos_tree, pst_type neg_tree ):
        self.radix = radix
        self.alphabet = strdup( alphabet )
        self.l_max = l_max
        self.pos_tree = pos_tree
        self.neg_tree = neg_tree
        
    def from_file( self, infile ):
        cdef int nradix
        cdef char* nalphabet
        cdef int nl_max
        
        my_load_pst( infile, &(self.alphabet), &(self.radix), &(self.l_max), &(self.pos_tree) )
        my_load_pst( infile, &(nalphabet), &(nradix), &(nl_max), &(self.neg_tree) )
        
        assert self.radix == nradix
        assert str( self.alphabet ) == str( nalphabet ) == str.join( '', map( chr, range( offset, offset + nradix ) ) )
        assert self.l_max == nl_max

    def get_radix( self ):
        return self.radix

    def score( self, string, int start=0, int length=-1 ):
        cdef double pos_log
        cdef double neg_log
        cdef char* s
        
        if length > 0:
            assert start + length < len( string )
        else:
            length = len( string )
            
        s = ints_to_chars( string )
        s[start+length+1] = 0
        
        pos_log = log10like_on_pst( self.alphabet, self.pos_tree, s + start )
        neg_log = log10like_on_pst( self.alphabet, self.neg_tree, s + start )
        
        free( s )
        
        return pos_log - neg_log

    def to_file( self, file ):
        my_save_pst( file, self.alphabet, self.radix, self.l_max, self.pos_tree )
        my_save_pst( file, self.alphabet, self.radix, self.l_max, self.neg_tree )
        
    def __dealloc__( self ):
        free( self.alphabet )
        # FIXME: must free trees!

def from_file( f ):
    rval = PSTModel()
    rval.from_file( f )
    return rval

# save_pst always clobbers the file, we want to be a bit nicer
cdef my_save_pst( outfile, char* alphabet, int radix, int l_max, pst_type tree ):
    cdef FILE* fp
    fp = PyFile_AsFile( outfile )
    fprintf( fp, "%d %s %d ", radix, alphabet, l_max )
    rec_save_pst( fp, radix, tree )
    
cdef my_load_pst( infile, char **alphabetp, int *radixp, int *l_maxp, pst_type *treep ):
    cdef FILE* fp
    fp = PyFile_AsFile( infile )
    fscanf( fp, "%d ", radixp )
    alphabetp[0] = <char *> malloc ( (radixp[0]) + 1 )
    if alphabetp[0] == NULL: raise "Malloc failed (alphabet)"
    fscanf( fp, "%s %d ", alphabetp[0], l_maxp );
    rec_load_pst( fp, radixp[0], treep );

cdef char* ints_to_chars( string ):
    cdef int s_len
    cdef int *s_buf
    cdef char *rval
    cdef int i
    PyObject_AsReadBuffer( string, <void **> &s_buf, &s_len )
    s_len = s_len / sizeof( int )
    rval = <char*> malloc( ( s_len + 1 ) * sizeof( char ) )
    rval[s_len] = 0
    
    for i from 0 <= i < s_len:
        rval[i] = offset + s_buf[i]
    return rval

def train( int order, int radix, pos_strings, neg_strings ):
    
    cdef char * alphabet
    cdef pst_type pos_tree, neg_tree
    cdef PSTModel rval
    
    # Create alphabet string
    s = str.join( '', map( chr, range( offset, offset + radix ) ) )
    alphabet = strdup( s )
    
    # Train two models
    pos_tree = train_pst( order, radix, alphabet, pos_strings )
    neg_tree = train_pst( order, radix, alphabet, neg_strings )
    
    # Set training parameters to something reasonable
    # FIXME: make configurable
    
    global p_min, alpha, gamma_min, p_ratio, L_max
    p_min = 0.001
    alpha = 0
    gamma_min = 0.001
    p_ratio = 1.05
    L_max = 5
    
    # Return PSTModel object
    rval = PSTModel()
    rval.init( radix, alphabet, L_max, pos_tree, neg_tree )
    return rval

cdef pst_type train_pst( int order, int radix, char* alphabet, strings ):
    
    cdef int strings_len 
    cdef char ** strings_ptr
    cdef int * lengths_ptr
    cdef int i
    
    cdef pst_type tree
    cdef int pot_nodes
    
    global gamma_min
    
    # Convert training data to character strings
    ## The pst code wastes one string pointer at the front for some reason
    strings_len = len( strings )
    strings_ptr = <char **> malloc( ( strings_len + 1 ) * sizeof( char* ) )
    lengths_ptr = <int *> malloc( ( strings_len + 1 ) * sizeof( int ) )
    strings_ptr[0] = NULL
    lengths_ptr[0] = 0
    for i from 0 <= i < strings_len:
        lengths_ptr[i+1] = len( strings[i] )
        strings_ptr[i+1] = ints_to_chars( strings[i] )
        
    # Build new tree 
    tree = NULL
    tree = grow_tree1( alphabet, radix, strings_len, strings_ptr, lengths_ptr, &pot_nodes, tree )
    add_probs_pst( alphabet, radix, L_max, strings_len, strings_ptr, tree )
    ## FIXME: this seems to screw things up (scores == nan)
    ## smoothen_pst( radix, gamma_min, tree )
    
    # Free training set
    for i from 0 <= i < strings_len + 1:
        free( strings_ptr[i] )
    free( strings_ptr )
    free( lengths_ptr )
    
    # Return tree
    return tree