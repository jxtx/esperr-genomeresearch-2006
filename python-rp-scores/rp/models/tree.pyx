cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, void **buffer, int *buffer_len) except -1
    object PyFloat_FromDouble( double v )

cdef extern from "stdlib.h":    
    void* malloc( int )
    void free( void* )
    
cdef extern from "stdio.h":
    int printf( char*, ... )

cdef extern from "math.h":
    double log( double )
    
cdef struct Node:
    float* vals
    Node** children
    
cdef Node* new_node( int radix ):
    cdef Node* node
    cdef int i
    node = <Node*> malloc( sizeof( Node ) )
    if node == NULL: return NULL
    node.children = <Node**> malloc( radix * sizeof( Node* ) )
    node.vals = <float*> malloc( radix * sizeof( float ) )
    for i from 0 <= i < radix: 
        node.children[i] = NULL
        node.vals[i] = 0
    return node
    
cdef Node* clone_node( int radix, Node* node ):
    cdef Node* rval
    cdef int i
    rval = new_node( radix )
    if rval == NULL: return NULL
    for i from 0 <= i < radix:
        rval.vals[i] = node.vals[i]

cdef free_node( Node* node, int radix ):
    cdef int i
    if node == NULL: return
    for i from 0 <= i < radix: 
        if node.children[i] != NULL:
            free( node.children[i] )
    free( node.vals )
    free( node )
    
cdef int count( int order, int radix, int* s, int slen, Node* root ):
    cdef Node *cur
    cdef int i, j, symbol, prefix_symbol
    # Loop over every prefix of the string
    for i from 0 <= i < slen:
        cur = root
        symbol = s[ i ]
        root.vals[ symbol ] = root.vals[ symbol ] + 1
        # Walk back along the prefix up to 'order' symbols
        for j from 1 <= j <= order:
            if i - j < 0: break
            prefix_symbol = s[i-j]
            if cur.children[prefix_symbol] == NULL:
                cur.children[prefix_symbol] = new_node( radix )
                if cur.children[prefix_symbol] == NULL: return 0
            cur = cur.children[prefix_symbol]
            cur.vals[ symbol ] = cur.vals[ symbol ] + 1
    return 1

cdef int fill_in_counts( int order, int radix, Node* counts, strings ):
    cdef int* buf
    cdef int buf_len
    for string in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        if count( order, radix, buf, buf_len / sizeof( int ), counts ) == 0:
            return 0
    return 1

cdef to_probs( int radix, Node* node ):
    cdef int i, total, some_zero
    total = 0
    some_zero = 0
    for i from 0 <= i < radix:
        total = total + node.vals[i]
        some_zero = some_zero or ( node.vals[i] == 0 )
    if some_zero: 
        total = total + radix
    for i from 0 <= i < radix:
        if some_zero:
            node.vals[i] = node.vals[i] + 1
        node.vals[i] = node.vals[i] / total
    for i from 0 <= i < radix:
        if node.children[i] != NULL:
            to_probs( radix, node.children[i] )

cdef Node* to_scores( int radix, Node* probs1, Node* probs2 ):
    cdef Node* rval
    cdef int i, j
    rval = new_node( radix )
    if rval == NULL:
        return NULL
    for i from 0 <= i < radix:
        rval.vals[i] = ( log( probs1.vals[i] ) - log( probs2.vals[i] ) )
    for i from 0 <= i < radix:
        if probs1.children[i] == NULL and probs2.children[i] == NULL:
            continue
        elif probs1.children[i] == NULL:
            probs1.children[i] = new_node( radix )
            if probs1.children[i] == NULL:
                free_node( rval, radix )
                return NULL
            for j from 0 <= j < radix:
                probs1.children[i].vals[j] = probs1.vals[j]
        elif probs2.children[i] == NULL:
            probs2.children[i] = new_node( radix )
            if probs2.children[i] == NULL:
                free_node( rval, radix )
                return NULL
            for j from 0 <= j < radix:
                probs2.children[i].vals[j] = probs2.vals[j]
        rval.children[i] = to_scores( radix, probs1.children[i], probs2.children[i] )
    return rval

cdef score_string( int order, int radix, Node* tree, int* text, int start, int length ):
    cdef int i, j, good, words, symbol, prefix_symbol
    cdef float score
    cdef Node* cur
    score = 0
    words = 0
    for i from start <= i < start + length:
        # First, is it a valid word -- may be too stringent but consistent with fixed order
        if i - order < 0: continue
        good = 1
        for j from 0 <= j <= order:
            if text[i-j] < 0: 
                good = 0
        if good == 0: continue
        # Now walk back and score
        cur = tree        
        symbol = text[i]
        for j from 1 <= j <= order:
            if i - j < 0: break
            prefix_symbol = text[i-j]
            # If we can't go back any further, use this context
            if cur.children[prefix_symbol] == NULL: break 
            # Otherwise we step back another symbol
            cur = cur.children[prefix_symbol]
        score = score + cur.vals[ symbol ]
        words = words + 1
    if words > 0:
        return score / <float> words
    else:
        return None

cdef print_node( int level, int radix, Node* node ):
    print "Node: [", 
    for i in range( 0, radix ): print str(node.vals[i]) + ",",
    print "],"
    for i in range( 0, radix ):
        if node.children[i] != NULL:
            for j from 0 <= j < level: printf( "      " );
            print i, "->",
            print_node( level+1, radix, node.children[i] )

cdef class Model:
    cdef int order
    cdef int radix
    cdef Node* tree

    cdef init( self, int order, int radix, Node* tree ):
        self.order = order
        self.radix = radix
        self.tree = tree

    def get_order( self ):
        return self.order

    def get_radix( self ):
        return self.radix
    
    def score( self, string, int start=0, int length=-1 ):
        cdef int* buf
        cdef int buf_len
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        buf_len = buf_len / sizeof( int )
        if length < 0:
            length = buf_len
        else:
            assert start + length <= buf_len
        s = score_string( self.order, self.radix, self.tree, buf, start, length )
        if s is None: raise "No valid data in region to be scored"
        return s

    def to_file( self, file ):
        raise "Not Yet Implemented"
        #file.write( "order: " + str( self.order ) + "\n" )
        #file.write( "radix: " + str( self.radix ) + "\n" )

    def __dealloc__( self ):
        free_node( self.tree, self.radix )
    
def train( int order, int radix, pos_strings, neg_strings ):
    cdef Node *pos_node, *neg_node, *scores
    cdef Model rval

    pos_node = new_node( radix )
    if pos_node == NULL: raise "Malloc failed creating pos root"
    if fill_in_counts( order, radix, pos_node, pos_strings ) == 0:
        raise "Malloc failed, filling in pos counts"
    to_probs( radix, pos_node )

    neg_node = new_node( radix )
    if neg_node == NULL: raise "Malloc failed creating neg root"
    if fill_in_counts( order, radix, neg_node, neg_strings ) == 0:
        raise "Malloc failed, filling in neg counts"
    to_probs( radix, neg_node )

    scores = to_scores( radix, pos_node, neg_node )
    if scores == NULL: raise "Malloc failed build score tree"

    free_node( pos_node, radix )
    free_node( neg_node, radix )

    rval = Model()
    rval.init( order, radix, scores )

    return rval


def test():
    cdef int text1[6]
    cdef int text2[6]
    cdef Node* node1
    cdef Node* node2
    cdef Node* scores
    
    text1[0] = 0
    text1[1] = 0
    text1[2] = 0
    text1[3] = 1
    text1[4] = 2
    text1[5] = 2

    text2[0] = 1
    text2[1] = 0
    text2[2] = 2
    text2[3] = 1
    text2[4] = 1
    text2[5] = 2
    
    node1 = new_node( 3 )
    node2 = new_node( 3 )

    count( 3, 3, text1, 6, node1 )    
    to_probs( 3, node1 )
    print_node( 1, 3, node1 )

    count( 3, 3, text2, 6, node2 )    
    to_probs( 3, node2 )
    print_node( 1, 3, node2 )

    scores = to_scores( 3, node1, node2 )
    print_node( 1, 3, scores )
