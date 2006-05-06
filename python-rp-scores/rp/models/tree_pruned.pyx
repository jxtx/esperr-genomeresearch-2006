"""
Variable order markov model stored as a tree with pruning of contexts seen
less than a fixed number of times (N) and discount smoothing that distributes
a fixed amount of probability mass (D) across nodes with any zero counts
"""

from elementtree.ElementTree import ElementTree, Element, parse
import sys

cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, void **buffer, int *buffer_len) except -1
    int PyObject_AsWriteBuffer(object obj, void **buffer, int *buffer_len) except -1
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
    float total
    float* vals
    Node** children

cdef Node* new_node( int radix ):
    """
    Creates a new empty node capable of holding 'radix' values and children.
    """
    cdef Node* node
    cdef int i
    node = <Node*> malloc( sizeof( Node ) )
    if node == NULL: return NULL
    node.children = <Node**> malloc( radix * sizeof( Node* ) )
    node.vals = <float*> malloc( radix * sizeof( float ) )
    for i from 0 <= i < radix: 
        node.children[i] = NULL
        node.vals[i] = 0
    node.total = 0
    return node
    
cdef free_node( Node* node, int radix ):
    """
    Free a node and recursively free all of its children.
    """
    cdef int i
    if node == NULL: return
    for i from 0 <= i < radix: 
        if node.children[i] != NULL:
            free_node( node.children[i], radix )
    free( node.children )
    free( node.vals )
    free( node )
   
cdef int count_for_order( int order, int max_order, int radix, int* s, int slen, Node* root ):
    """
    Fill in counts just for the specified order. Extend the context tree as
    needed to accomodate the counts. Starts at position 'max_order' in the
    string of ints 's' and counts the next 'slen' positions.
    """
    cdef Node* cur
    cdef int i, j, symbol, prefix_symbol
    # If order is 0, just fill in the counts for the root node
    if order == 0:
        for i from order <= i < slen:
            symbol = s[i]
            if symbol < 0: continue
            root.vals[ symbol ] = root.vals[ symbol ] + 1
        return 1
    # With longer orders, we walk the tree for the context at each position
    for i from order <= i < slen:
        # Start the cursor at the root of the tree
        cur = root
        # Get the symbol we are 'transitioning to'
        symbol = s[i]
        # If it is an invalid character skip this word
        if symbol < 0: continue
        # Walk down to the node just above the one that will hold this orders counts
        for j from 1 <= j < order:
            # Grab the next symbol in the current context
            prefix_symbol = s[i-j]
            # Skip this entire word if an invalid character was encountered
            if prefix_symbol < 0: cur = NULL
            # Move our cursor to the child corresponding to the context symbol
            else: cur = cur.children[ prefix_symbol ]
            # If the child was null, stop here, we will not extend
            if cur == NULL: break
        # If any context was missing on the way back we do not extend
        if cur == NULL: continue
        # Extend to last node of context and increment
        prefix_symbol = s[i-order]
        # Skip invalid characters as usual
        if prefix_symbol < 0: continue
        # Create the context we have just seen if needed
        if cur.children[prefix_symbol] == NULL:
            cur.children[prefix_symbol] = new_node( radix )
            if cur.children[prefix_symbol] == NULL: return 0
        # Increment the count for the context / symbol we just observed
        cur.children[prefix_symbol].vals[symbol] = cur.children[prefix_symbol].vals[symbol] + 1
        cur.children[prefix_symbol].total = cur.children[prefix_symbol].total + 1
    return 1

cdef to_probs( int radix, Node* node, Node* parent, float discount ):
    """
    Recursively visit the tree and convert counts to probabilities. If at least
    some symbols have zero counts the 'discount' will be used to allocate some
    probability according to the parents distribution, or uniformly if 'parent'
    is NULL. 
    """
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
#     if num_zero == 0:
#         # No zero nodes, just straight probabilities
#         for i from 0 <= i < radix: 
#             node.vals[i] = node.vals[i] / total        
#     else:
#         # Spread discount among nodes
#         for i from 0 <= i < radix:
#             if parent == NULL:
#                 node.vals[i] = ( (1-discount) * (node.vals[i]/total) ) + ( discount / radix )
#             else:
#                 node.vals[i] = ( (1-discount) * (node.vals[i]/total) ) \
#                                + ( (discount) * parent.vals[i] )       
     
    # Always smooth
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

cdef to_logs( int radix, Node* node  ):
    cdef int i
    for i from 0 <= i < radix:
        node.vals[ i ] = log( node.vals[ i ] )
        if node.children[i] != NULL:
            to_logs( radix, node.children[i] )

cdef prune( int current_depth, int desired_depth, int radix, Node* node, int N ):
    """
    Prune from the tree any context at level 'desired_depth' that has been
    seen less than 'N' times. This should be run AFTER counts have been
    filled in for level 'desired_depth'. 
    """
    cdef int i, j
    # Can't prune the root node
    if desired_depth == 0: 
        return
    # We do stop visiting at the parents of the desired depth and
    # examine each child (since we have to touch the parent to remove
    # the reference this is easier than going all the way down).
    elif desired_depth == current_depth + 1:
        for i from 0 <= i < radix:
            # Check if the child exists and if so if it should be pruned
            if node.children[i] != NULL and node.children[i].total < N:
                # Remove the child
                free_node( node.children[i], radix )
                node.children[i] = NULL
    # Not at parent of 'desired_depth', so recursively visit each child
    else:
        for i from 0 <= i < radix:
            if node.children[i] != NULL:
                prune( current_depth + 1, desired_depth, radix, node.children[i], N )
                
cdef Node* to_scores( int radix, Node* probs1, Node* probs2 ):
    """
    Take two trees of nodes that have been converted to probs and create a tree
    of log-odds scores from them. The resulting tree toplogy is the union of the
    two prob trees (every context that is in either of the prob trees has a node
    in the score tree).
    """
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
        # In this case, if we observe the node in either set we extend the
        # set in which it was not observed. FIXME: no reason to do all this
        # copying here, we should be able to just reuse the parent node
        # by passing down some flag
# ------- First strategy, require it to be observed in either to extend 
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
# ------- Alternate strategy, require it to be observed in both sets to extend 
# ------- if this works, integrate it into probs to make faster  
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

cdef score_string( int order, int radix, Node* tree, int* text, int start, int length, int norm ):
    """
    Score a string of integers using a log-odds scoring tree 'tree' as produced
    by to_scores. For each position in string[start:start+length] walk the tree
    to the deepest possible node that can be reached from the context at that
    position and add the score for the symbol at that position to the running
    total. Return the total score / the number of positions scored, or 'None'
    if there were no valid words in the string.
    """
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
        if norm:
            return score / <float> words
        else:
            return score
    else:
        return None

cdef score_string_positions( int order, int radix, Node* tree, int* text, float* target, int start, int length ):
    """
    Fill into 'target' scores for each position in 'text' If there is no valid score
    for a position the value in 'target' wll not be changed
    """
    cdef int i, j, good, words, symbol, prefix_symbol
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
        target[i] = cur.vals[ symbol ]

cdef print_node( int level, int radix, Node* node ):
    """
    For debugging, dumps a tree display of a node and all of its children
    including their values.
    """
    print "Node: [", 
    for i in range( 0, radix ): print str(node.vals[i]) + ",",
    print "],"
    for i in range( 0, radix ):
        if node.children[i] != NULL:
            for j from 0 <= j < level: printf( "      " );
            print i, "->",
            print_node( level+1, radix, node.children[i] )

cdef int count_transition_probs( int radix, Node* node):
    cdef int rval
    rval = radix
    for i from  0 <= i < radix:
        if node.children[i] != NULL:
            rval = rval + count_transition_probs( radix, node.children[i] ) - 1
    return rval

cdef class Model:
    """
    Python class that provides the 'Model' interface expected by various RP
    programs. Wraps a 'Node' tree as produced by to_scores, and allows
    python int buffers to be scored with the model. Can be saved to
    a file (in XML using ElementTree).
    """
    cdef int order
    cdef int radix
    cdef Node* tree

    cdef init( self, int order, int radix, Node* tree ):
        """
        Create new model having max order 'order', alphabet size 'radix' and
        log-odds score tree 'tree'.
        """
        self.order = order
        self.radix = radix
        self.tree = tree

    def get_order( self ):
        return self.order

    def get_radix( self ):
        return self.radix
    
    def score( self, string, int start=0, int length=-1 ):
        """
        Score string[start:start+length] under the model. If start/length are not
        specified they default to 0 and the length of the buffer respectively.
        """
        cdef int* buf
        cdef int buf_len
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        buf_len = buf_len / sizeof( int )
        if length < 0:
            length = buf_len
        else:
            assert start + length <= buf_len
        s = score_string( self.order, self.radix, self.tree, buf, start, length, 1 )
        # if s is None: raise "No valid data in region to be scored"
        return s

    def score_positions( self, string, target ):
        """
        Score string[start:start+length] under the model. If start/length are not
        specified they default to 0 and the length of the buffer respectively.
        """
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
        score_string_positions( self.order, self.radix, self.tree, buf, t_buf, 0, buf_len )
        
    def to_file( self, file ):
        """
        Write the model to the open file-like object 'file'.
        """
        root = Element( "root", order=str(self.order), radix=str(self.radix) )
        root.append( node_to_element( self.order, self.radix, -1, self.tree ) )
        ElementTree( root ).write( file )
        
    def count_transition_probs( self ):
        return count_transition_probs( self.radix, self.tree )

    def __dealloc__( self ):
        """
        When this object is destroyed it needs to free the tree which was created
        outside the python memory allocator.
        """
        free_node( self.tree, self.radix )
        
cdef class ProbModel:
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
        s = score_string( self.order, self.radix, self.tree, buf, start, length, 1 )
        # if s is None: raise "No valid data in region to be scored"
        return s
    def to_file( self, file ):
        root = Element( "root", order=str(self.order), radix=str(self.radix) )
        root.append( node_to_element( self.order, self.radix, -1, self.tree ) )
        ElementTree( root ).write( file )
    ## def count_transition_probs( self ):
    ##     return count_transition_probs( self.radix, self.tree )
    def __dealloc__( self ):
        free_node( self.tree, self.radix )
        
def prob_train( int order, int radix, pos_strings, **kwargs ):
    """
    NOTE: This one is just for a one class model which can calculate a (log) prob
    """
    cdef Node *pos_node
    cdef ProbModel rval
    cdef int* buf
    cdef int buf_len
    cdef int i
    cdef int D, N
    # Convert discount parameter
    try: 
        d = float( kwargs['D'] )
    except: 
        d = .10
    assert 0.0 < d < 1.0, "Discount must be between 0 and 1"
    # Convert pruning parameter
    try: 
        N = int( kwargs['N'] )
    except: 
        N = 5
    assert N >= 0, "N must be non-negative"
    # Create root nodes for count/prob trees
    pos_node = new_node( radix )
    if pos_node == NULL: 
        raise "Malloc failed creating pos root"
    # Fill in counts one order at a time
    for i from 0 <= i <= order:
        # Need to loop over the training set completely for each order so we can
        # see the counts at that order and decide what nodes get extended
        for string in pos_strings:
            PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
            if count_for_order( i, order, radix, buf, buf_len / sizeof( int ), pos_node ) == 0:
                raise "Failed while adding to pos counts"
        prune( 0, i, radix, pos_node, N )        
    # Allow model to be dumped to a file for debugging    
    if 'dump' in kwargs:
        to_file( order, radix, pos_node, "pos_node.debug" )
    # Convert counts to probs    
    to_probs( radix, pos_node, NULL, d )
    to_logs( radix, pos_node )
    # Create and return ProbsModel object
    rval = ProbModel()
    rval.init( order, radix, pos_node )
    return rval

def train( int order, int radix, pos_strings, neg_strings, **kwargs ):
    """
    Train a new model with max order 'order' and alphabet size 'radix' from the
    positive/negative training data in 'pos_strings'/'neg_strings'. Valid keyword
    arguments are 'D' for affecting the discount amount used in 'to_probs' and 'N'
    for affecting the number of times a node must be seen to avoid pruning in
    'prune'.
    """
    cdef Node *pos_node, *neg_node, *scores
    cdef Model rval
    cdef int* buf
    cdef int buf_len
    cdef int i
    cdef int D, N

    # Convert discount parameter
    try: d = float( kwargs['D'] )
    except: d = .10
    assert 0.0 < d < 1.0, "Discount must be between 0 and 1"

    # Convert pruning parameter
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
        prune( 0, i, radix, pos_node, N )
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
    """
    Read a 'Model' from the open file-like object 'f'
    """
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

def prob_from_file( f ):
    """
    Read a 'Model' from the open file-like object 'f'
    """
    cdef int order
    cdef int radix
    cdef Node* node
    cdef ProbModel rval

    tree = parse( f )
    root = tree.getroot()
    order = int( root.get( 'order' ) )
    radix = int( root.get( 'radix' ) )

    node = parse_node( root.find( "node" ), radix )

    rval = ProbModel()
    rval.init( order, radix, node )
    return rval

cdef to_file( int order, int radix, Node* node, filename ):
    """
    Write the tree under 'node' to the file named 'filename'. This is used
    only for debugging purposes. You can dump a count/prob/score tree at
    any time with this and view it with 'rp_dump_tree.py'.
    """
    root = Element( "root", order=str(order), radix=str(radix) )
    root.append( node_to_element( order, radix, -1, node) )
    out = open( filename, 'w' )
    ElementTree( root ).write( out )
    out.close()
    
cdef node_to_element( int order, int radix, int symbol, Node* node ):
    """
    Build an XML/ElementTree element for a node (recursively handling
    children). Used for writing to files
    """
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
    """
    Parse an XML/ElementTree element into a Node object, for reading
    models from files.
    """
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
