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
    for i from 0 <= i < radix:
        rval.vals[i] = node.vals[i]

cdef free_node( Node* node, int radix ):
    cdef int i
    for i from 0 <= i < radix: 
        if node.children[i] != NULL:
            free( node.children[i] )
    free( node.vals )
    free( node )
    
cdef count( int order, int radix, int* s, int slen, Node* root ):
    cdef Node *cur
    cdef int i, j, symbol, prefix_symbol
    # Loop over every prefix of the string
    for i from 0 <= i < slen:
        cur = root
        symbol = s[ i ]
        root.vals[ symbol ] = root.vals[ symbol ] + 1
        # Walk back along the prefix up to 'order' symbols
        for j from 1 <= j < order:
            if i - j < 0: break
            prefix_symbol = s[i-j]
            if cur.children[prefix_symbol] == NULL:
                cur.children[prefix_symbol] = new_node( radix )
            cur = cur.children[prefix_symbol]
            cur.vals[ symbol ] = cur.vals[ symbol ] + 1

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
    for i from 0 <= i < radix:
        rval.vals[i] = ( log( probs1.vals[i] ) - log( probs2.vals[i] ) )
    for i from 0 <= i < radix:
        if probs1.children[i] == NULL and probs2.children[i] == NULL:
            continue
        elif probs1.children[i] == NULL:
            probs1.children[i] = new_node( radix )
            for j from 0 <= j < radix:
                probs1.children[i].vals[j] = probs1.vals[j]
        elif probs2.children[i] == NULL:
            probs2.children[i] = new_node( radix )
            for j from 0 <= j < radix:
                probs2.children[i].vals[j] = probs2.vals[j]
        rval.children[i] = to_scores( radix, probs1.children[i], probs2.children[i] )
    return rval

cdef print_node( int level, int radix, Node* node ):
    
    print "Node: [", 
    for i in range( 0, radix ): print str(node.vals[i]) + ",",
    print "],"
    for i in range( 0, radix ):
        if node.children[i] != NULL:
            for j from 0 <= j < level: printf( "      " );
            print i, "->",
            print_node( level+1, radix, node.children[i] )
    
def train( int order, int radix, pos_strings, neg_strings ):
    cdef Node *pos_node *neg_node *scores

    pos_node = new_node()
    fill_in( order, radix, pos_node, pos_strings )
    to_probs( radix, pos_node )

    neg_node = new_node()
    fill_in( order, radix, neg_node, neg_strings )
    to_probs( radix, neg_node )

    scores = to_scores( radix, pos_node, neg_node )

    free_node( pos_node )
    free_node( neg_node )

    


cdef fill_in( int order, int radix, Node* counts, strings ):
    cdef int* buf
    cdef int buf_len
    for string in strings:
        PyObject_AsReadBuffer( string, <void**> &buf, &buf_len )
        count( order, radix, buf, buf_len / sizeof( int ), counts )


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
