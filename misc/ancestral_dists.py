#!/usr/bin/env python2.4

"""
Construct ancestral probability distribution for a set of alignment columns

usage: %prog tree_mod species_names col_count_file

`tree_mod`: A tree model as produced by phyloFit from PHAST.
`species_names`: a comma seperated list of the species names in the order
                 they appear in the alignment
`col_file`: if present, a list of alignment columns and counts for which
            distributions will be computed. Otherwise distributions will
            be computed for all possible columns. Column will be taken from
            the first field in this file, and all other fields will be kept.
"""

from __future__ import division

import cookbook
import sys
from Numeric import *

from newick import tree
import tree_model

def main():

    tm = tree_model.from_file( open( sys.argv[1] ) )
    t = tree.parse_tree( tm.tree )

    names = sys.argv[2].split( ',' )

    nucs = tm.alphabet
    background = tm.background

    print "#", ','.join( names )

    if len( sys.argv ) > 2:
        
        for line in open( sys.argv[3] ):
            
            fields = line.rstrip( "\r\n" ).split( "\t" )
            col = fields[0]
            cols.append( ( fields[0], int( fields[1]) ) )   
            a = dict( zip( names, col ) )
            lik = felsen( t, a, tm )
            prob = lik_to_prob( t, lik, background, tm ) 
            print "\t".join( fields + ( ' '.join( map( str, prob ) ) ) )
    else:
        
        for rows in cookbook.cross_lists( *( [ nucs ] * len( names ) ) ):
            
            a = dict( zip( names, rows ) )
            lik = felsen( t, a, tm )
            prob = lik_to_prob( t, lik, background, tm ) 
            print ''.join( rows ), ' '.join( map( str, prob ) )

def felsen( edge, column, tm ):
    """Reconstruct likelihood using felsen's algorithm"""
    nucs = tm.alphabet
    if edge.__class__ == tree.Leaf:
        symbol = column[ edge.id ]
        if symbol == '?':
            # Uniform distribution
            ##return [ 1 / len( nucs ) ] # len( nucs ) 
            # Equilibrium (\pi) distribution
            return tm.background
        return [ int( symbol == nuc ) for nuc in nucs ] 
    else:
        # Traverse children and determine likelihoods
        l_children = [ felsen( child[0], column, tm ) for child in edge.edges ]
        l_self = []
        # Determine liklihood for each possible 'ancestral' sequence
        for i_y in range( len( nucs ) ):
            y = nucs[ i_y ]
            p_y = 1
            # For each child, sum over all paths
            for i in range( len( edge.edges ) ):
                child, bs, length = edge.edges[i]
                cl = l_children[i]
                psub = tm.matrix_for_time( length )
                p_y_c = 0
                for i_a in range( len( nucs ) ):
                    p_y_c += psub[i_y,i_a] * cl[ i_a ]
                # Product over all children
                p_y *= p_y_c
            l_self.append( p_y )
        return l_self

def lik_to_prob( edge, lik, p_parent, tm ):
    nucs = tm.alphabet
    p = []
    for i in range( len( nucs ) ):
        p.append( lik[i] * p_parent[i] )
    p_x = sum( p )
    # print >>sys.stderr, lik
    rval = [ top / p_x for top in p ]
    # print >>sys.stderr, rval
    return rval
        
if __name__ == "__main__":
    main()