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
    float fabsf( float )
    
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
   
cdef int count_for_order( int order, int max_order, int radix, int* s, int slen, Node* root ):
    """Fill in counts _just_ for the specified order."""
    cdef Node* cur
    cdef int i, j, symbol, prefix_symbol
    # If order is 0, just fill in the counts for the root node
    if order == 0:
        for i from max_order <= i < slen:
            symbol = s[i]
            if symbol < 0: continue
            root.vals[ symbol ] = root.vals[ symbol ] + 1
        return 1
    # With longer orders, we walk the tree for the context at each position
    for i from max_order <= i < slen:
        cur = root
        symbol = s[i]
        if symbol < 0: continue
        # Walk back until next-to-last step
        for j from 1 <= j < order:
            prefix_symbol = s[i-j]
            if prefix_symbol < 0: 
                cur = NULL
            else: 
                cur = cur.children[ prefix_symbol ]
            if cur == NULL: break
        # If any context was missing on the way back we do not extend
        if cur == NULL: continue
        # Extend to last node of context and increment
        prefix_symbol = s[i-order]
        if prefix_symbol < 0: continue
        # Create the context we have just seen if needed
        if cur.children[prefix_symbol] == NULL:
            cur.children[prefix_symbol] = new_node( radix )
            if cur.children[prefix_symbol] == NULL: return 0
        # Increment the count for the context / symbol we just observed
        cur.children[prefix_symbol].vals[symbol] = cur.children[prefix_symbol].vals[symbol] + 1     
    return 1

#cdef to_probs( int radix, Node* node ):
#    """Laplace version"""
#    cdef int i, total, some_zero
#    total = 0
#    some_zero = 0
#    for i from 0 <= i < radix:
#        total = total + node.vals[i]
#        some_zero = some_zero or ( node.vals[i] == 0 )
#    if some_zero: 
#        total = total + radix
#    for i from 0 <= i < radix:
#        if some_zero:
#            node.vals[i] = node.vals[i] + 1
#        node.vals[i] = node.vals[i] / total
#    for i from 0 <= i < radix:
#        if node.children[i] != NULL:
#            to_probs( radix, node.children[i] )

cdef to_probs( int radix, Node* node, Node* parent, float discount ):
    """Discount Version: Proportional take, Previous order give"""
    cdef int i, total, num_zero
    total = 0
    num_zero = 0
    # Determine total and number of zero nodes
    for i from 0 <= i < radix:
        if node.vals[i] == 0:
            num_zero = num_zero + 1
        else:
            total = total + node.vals[i]
    # Make probabilities
    if num_zero == 0:
        # No zero nodes, just straight probabilities
        for i from 0 <= i < radix: 
            node.vals[i] = node.vals[i] / total        
    else:
        # Spread discount among nodes
        # print "spreading", discount
        for i from 0 <= i < radix:
            if parent == NULL:
                node.vals[i] = ( (1-discount) * (node.vals[i]/total) ) + ( discount / radix )
            else:
                node.vals[i] = ( (1-discount) * (node.vals[i]/total) ) \
                               + ( (discount) * parent.vals[i] )       
    # Recursively visit children        
    for i from 0 <= i < radix:
        if node.children[i] != NULL:
            to_probs( radix, node.children[i], node, discount )

cdef prune( int current_depth, int desired_depth, int radix, Node* node, int N ):
    cdef int i, j
    # Can't prune the root node
    if desired_depth == 0: 
        return
    elif desired_depth == current_depth + 1:
        for i from 0 <= i < radix:
            # Simple count based pruning -- seen less than N times, prune
            if node.vals[i] < N and node.children[i] != NULL:
                free_node( node.children[i], radix )
                node.children[i] = NULL
    else:
        for i from 0 <= i < radix:
            if node.children[i] != NULL:
                prune( current_depth + 1, desired_depth, radix, node.children[i], N )
                
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
# ------ In this case, if we observe the node in _either_ set we extend
#        elif probs1.children[i] == NULL:
#            probs1.children[i] = new_node( radix )
#            if probs1.children[i] == NULL:
#                free_node( rval, radix )
#                return NULL
#            for j from 0 <= j < radix:
#                probs1.children[i].vals[j] = probs1.vals[j]
#        elif probs2.children[i] == NULL:
#            probs2.children[i] = new_node( radix )
#            if probs2.children[i] == NULL:
#                free_node( rval, radix )
#                return NULL
#            for j from 0 <= j < radix:
#                probs2.children[i].vals[j] = probs2.vals[j]
# ------- Alternate strategy, require it to be observed in both sets to extend -- if this works, integrate it into probs to make faster
        elif probs1.children[i] == NULL:
            free_node( probs2.children[i], radix )
            probs2.children[i] = NULL
            continue
        elif probs2.children[i] == NULL:
            free_node( probs1.children[i], radix )
            probs1.children[i] = NULL
            continue
# ------- End alternate

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
        root.append( node_to_element( self.order, self.radix, -1, self.tree ) )
        ElementTree( root ).write( file )

    def __dealloc__( self ):
        # sys.stderr.write( "freeing tree_prned_1.Model\n" ); sys.stderr.flush()
        free_node( self.tree, self.radix )

def train( int order, int radix, pos_strings, neg_strings, **kwargs ):
    cdef Node *pos_node, *neg_node, *scores
    cdef Model rval
    cdef int* buf
    cdef int buf_len
    cdef int i
    cdef int D, N

    
    # Convert keyword parameters
    try: d = float( kwargs['D'] )
    except: d = .10
    assert 0.0 < d < 1.0, "Discount must be between 0 and 1"
    
    try: N = int( kwargs['N'] )
    except: N = 5
    assert N >= 0, "N must be non-negative"

    # Create root nodes for count/prob trees
    pos_node = new_node( radix )
    if pos_node == NULL: raise "Malloc failed creating pos root"
    neg_node = new_node( radix )
    if neg_node == NULL: raise "Malloc failed creating neg root"

    # Fill in counts one order at a time
    for i from 0 <= i <= order:
        # Need to loop over the training set completely for each order so we can
        # see the counts at that order and decide what nodes get extended
        for string in pos_strings:
            PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
            if count_for_order( i, order, radix, buf, buf_len / sizeof( int ), pos_node ) == 0:
                raise "Failed while adding to pos counts"
        for string in neg_strings:
            PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
            if count_for_order( i, order, radix, buf, buf_len / sizeof( int ), neg_node ) == 0:
                raise "Failed while adding to neg counts"
        #to_file( order, radix, pos_node, str( i ) + ".before.pos_node.debug" )
        prune( 0, i, radix, pos_node, N )
        #to_file( order, radix, pos_node, str( i ) + ".after.pos_node.debug" )
        prune( 0, i, radix, neg_node, N )
    if 'dump' in kwargs:
        to_file( order, radix, pos_node, "pos_node.debug" )
        to_file( order, radix, neg_node, "neg_node.debug" )
    to_probs( radix, pos_node, NULL, d )
    to_probs( radix, neg_node, NULL, d )

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

cdef to_file( int order, int radix, Node* node, filename ):
    root = Element( "root", order=str(order), radix=str(radix) )
    root.append( node_to_element( order, radix, -1, node) )
    out = open( filename, 'w' )
    ElementTree( root ).write( out )
    out.close()
    
cdef node_to_element( int order, int radix, int symbol, Node* node ):
    e = Element( "node", symbol=str(symbol) )
    if node == NULL: return e

    v = Element( "vals" )
    s = []
    for i in range( radix ): s.append( str( node.vals[i] ) )
    v.text = ','.join( s )
    e.append( v )

    c = Element( "children" )
    for i in range( radix ):
        c.append( node_to_element( order, radix, i, node.children[i] ) )
    e.append( c )

    return e

cdef Node* parse_node( element, int radix ):
    cdef Node* rval
    cdef int index
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
    for child in element.find( "children" ):
        index = int( child.get( 'symbol' ) )
        assert index < radix
        rval.children[ index ] = parse_node( child, radix )
    return rval
