from elementtree.ElementTree import ElementTree, Element, parse
import sys

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
            free_node( node.children[i], radix )
    free( node.children )
    free( node.vals )
    free( node )
    
cdef int count( int order, int radix, int* s, int slen, Node* root ):
    cdef Node *cur
    cdef int o, i, j, symbol
    # First, fill in the order 0 counts
    for i from order <= i < slen:
        symbol = s[i]
        if symbol < 0: continue
        root.vals[ symbol ] = root.vals[ symbol ] + 1
    # Now fill in progressively longer orders
    for o from 1 <= o <= order:
        for i from order <= i < slen:
            cur = root
            symbol = s[i]
            if symbol < 0: continue
            # Walk back until next-to-last step
            for j from 1 <= j < o:
                prefix_symbol = s[i-j]
                if prefix_symbol < 0:
                    cur = NULL
                else:
                    cur = cur.children[ prefix_symbol ]
                if cur == NULL: break
            if cur == NULL: continue
            # Extend to last node of context and increment
            prefix_symbol = s[i-o]
            if prefix_symbol < 0: continue
            # Pruning, if we only saw this context N times on the previous pass, don't extend
            if cur.vals[prefix_symbol] < radix:
                continue
            if cur.children[prefix_symbol] == NULL:
                cur.children[prefix_symbol] = new_node( radix )
                if cur.children[prefix_symbol] == NULL: return 0
            cur.children[prefix_symbol].vals[symbol] = cur.children[prefix_symbol].vals[symbol] + 1    
    return 1

cdef int fill_in_counts( int order, int radix, Node* counts, strings ):
    cdef int* buf
    cdef int buf_len
    for string in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        if count( order, radix, buf, buf_len / sizeof( int ), counts ) == 0:
            return 0
    return 1

cdef not_to_probs( int radix, Node* node ):
    """Laplace version"""
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

cdef to_probs( int radix, Node* node ):
    """Discount Version"""
    cdef int i, total, some_zero
    cdef float discount, fudge_zero, fudge_nonzero
    total = 0
    num_zero = 0
    # Determine total and number of zero nodes
    for i from 0 <= i < radix:
        if node.vals[i] == 0:
            num_zero = num_zero + 1
        else:
            total = total + node.vals[i]
    # Spread discount among nodes
    discount = 1.0 / radix 
    if num_zero > 0:
        fudge_zero = discount / num_zero
        fudge_nonzero = (1-discount)
    else:
        fudge_zero = 0
        fudge_nonzero = 1
    # Now do probs
    for i from 0 <= i < radix:
        if node.vals[i] == 0:
            node.vals[i] = fudge_zero
        else:
            node.vals[i] = ( node.vals[i] / total ) * fudge_nonzero
    # Recursively visit children        
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
        root = Element( "root", order=str(self.order), radix=str(self.radix) )
        root.append( self.node_to_element( self.tree ) )
        ElementTree( root ).write( file )

    cdef node_to_element( self, Node* node ):
        e = Element( "node" )
        if node == NULL: return e

        v = Element( "vals" )
        s = []
        for i in range( self.radix ): s.append( str( node.vals[i] ) )
        v.text = ','.join( s )
        e.append( v )

        c = Element( "children" )
        for i in range( self.radix ):
            c.append( self.node_to_element( node.children[i] ) )
        e.append( c )

        return e

    def __dealloc__( self ):
        # sys.stderr.write( "freeing tree_prned_1.Model\n" ); sys.stderr.flush()
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

def from_file( f ):
    cdef int order
    cdef int radix
    cdef Node* node
    cdef Model rval

    tree = parse( f )
    root = tree.getroot()
    order = int( root.get( 'order' ) )
    radix = int( root.get( 'radix' ) )

    node = parse_node( root.find( "node" ), radix )

    rval = Model()
    rval.init( order, radix, node )
    return rval

cdef Node* parse_node( element, int radix ):
    cdef Node* rval
    # Empty element corresponds to NULL pointer
    if not element: return NULL
    # Create new node
    rval = new_node( radix )
    # Parse vals
    vals = map( float, element.find( "vals" ).text.split( ',' ) )
    assert len( vals ) == radix
    for i in range( radix ):
        rval.vals[i] = vals[i]
    # Parse children
    children = element.find( "children" )
    assert len( children ) == radix
    for i in range( radix ):
        rval.children[i] = parse_node( children[i], radix )
    return rval
