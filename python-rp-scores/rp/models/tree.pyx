cdef extern from "stdlib.h":    
    void* malloc( int )
    void free( void* )
    
cdef extern from "stdio.h":
    int printf( char*, ... )
    
cdef struct Node:
    int* counts
    Node** children
    
    
cdef Node* new_node( int radix ):
    cdef Node* node
    cdef int i
    node = <Node*> malloc( sizeof( Node ) )
    node.children = <Node**> malloc( radix * sizeof( Node* ) )
    node.counts = <int*> malloc( radix * sizeof( int ) )
    for i from 0 <= i < radix: 
        node.children[i] = NULL
        node.counts[i] = 0
    return node
    
cdef free_node( Node* node, int radix ):
    cdef int i
    for i from 0 <= i < radix: 
        if node.children[i] != NULL:
            free( node.children[i] )
    free( node )
    
cdef Node* count( int order, int radix, int* s, int slen ):
    cdef Node *root, *cur
    cdef int i, j, ch
    # Initialize root
    root = new_node( radix )
    # Loop over every prefix of the string
    for i from 0 <= i < slen:
        root.count = root.count + 1
        # Walk back along the prefix up to 'order' symbols
        cur = root
        for j from 0 <= j < order:
            if i - j < 0: break
            ch = s[i-j]
            if cur.children[ch] == NULL:
                cur.children[ch] = new_node( radix )
            cur = cur.children[ch]
            cur.count = cur.count + 1
    return root     
    
cdef print_node( int level, int radix, Node* node ):
    cdef int i, j
    printf( "Node, count=%d\n", node.count )
    for i from 0 <= i < radix:
        if node.children[i] != NULL:
            for j from 0 <= j < level: printf( " " );
            printf( "%d -> ", i )
            print_node( level+1, radix, node.children[i] )
       
    
def test():
    cdef int text[6]
    cdef Node* node
    
    text[0] = 0
    text[1] = 0
    text[2] = 0
    text[3] = 1
    text[4] = 2
    text[5] = 2
    
    node = count( 3, 3, text, 6 )
    
    print_node( 0, 3, node )
