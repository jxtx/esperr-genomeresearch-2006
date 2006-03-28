#!/usr/bin/env python2.4

from __future__ import division

import sys

def read_dend( f ):
    by_stage = dict()
    last = None
    for i, line in enumerate( f ): 
        fields = line.split()
        left, right = map( int, fields[1:3] )
        stage = i + 1
        sim = float( fields[3] )
        # Integers less than zero indicate rows of the original distance 
        # matrix, greate than zero indicate mergers from previous stages
        if right >= 0:
            right = by_stage[ right ]
        else:
            right = right * -1
        if left >= 0:
            left = by_stage[ left ]
        else:
            left = left * -1
        by_stage[ stage ] = last = ( left, right, sim )
    # The last stage should be the root
    return last
    
class DendCutter( object ):
    def __init__( self ):
        self.clusters = {}
        self.last_cluster_id = 0
        self.max_element_id = 0
    def cut_dend( self, node, level ):
        if type( node ) == int:
            self.max_element_id = max( self.max_element_id, node )
            self.clusters[ node ] = self.last_cluster_id
            self.last_cluster_id += 1
        else:        
            left, right, sim = node
            if sim <= level:
                self.cut_dend( left, level )
                self.cut_dend( right, level )
            else:
                new_cluster_id = self.last_cluster_id
                self.last_cluster_id += 1
                self.fill_cluster( node, new_cluster_id )
    def fill_cluster( self, node, new_cluster_id ):
        left, right, sim = node 
        # Left subtree   
        if type( left ) == int:
            self.max_element_id = max( self.max_element_id, left )
            self.clusters[ left ] = new_cluster_id
        else:
            self.fill_cluster( left, new_cluster_id )
        # Right subtree
        if type( right ) == int:
            self.max_element_id = max( self.max_element_id, right )
            self.clusters[ right ] = new_cluster_id
        else:
            self.fill_cluster( right, new_cluster_id )
    
     
def main():
    root = read_dend( sys.stdin )
    level = float( sys.argv[1] )
    dc = DendCutter()
    dc.cut_dend( root, level )
    for i in range( dc.max_element_id ):
        print dc.clusters[i+1]
        
if __name__ == "__main__": main()