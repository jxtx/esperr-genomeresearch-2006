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

from bx import cookbook
import sys
from numpy import *
from scipy.linalg import expm

import bx.phylo.newick
import bx.phylo.phast

def main():

    tm = bx.phylo.phast.TreeModel.from_file( open( sys.argv[1] ) )
    t = bx.phylo.newick.newick_parser.parse_string( tm.tree )

    names = sys.argv[2].split( ',' )

    nucs = tm.alphabet
    background = tm.background

    print "#", ','.join( names )

    if len( sys.argv ) > 3:
       
        cols = []
	   
        for line in open( sys.argv[3] ):
            
            fields = line.rstrip( "\r\n" ).split()
            col = fields[0]
            # HACK!
            col = col.replace( "N", "*" )
            cols.append( ( fields[0], int( fields[1]) ) )   
            a = dict( zip( names, col ) )
            lik = felsen( t, a, tm )
            # print >>sys.stderr, col
            prob = lik_to_prob( lik, background, tm ) 
            print "\t".join( fields + [ ' '.join( map( str, prob ) ) ] )
    else:
        
        for rows in cookbook.cross_lists( *( [ nucs ] * len( names ) ) ):
            
            a = dict( zip( names, rows ) )
            lik = felsen( t, a, tm )
            prob = lik_to_prob( lik, background, tm ) 
            print "\t".join( [ ''.join(rows), "?" , ' '.join( map( str, prob ) ) ] )

matrix_by_time_cache = {}

def felsen( node, column, tm ):
    """Reconstruct likelihood using felsen's algorithm"""
    nucs = tm.alphabet
    # Is it a leaf node?
    if node.edges is None:
        symbol = column[ node.label ]
        if symbol == '*':
            ## Uniform distribution
            # return [ 1 / len( nucs ) ] # len( nucs ) 
            ## Equilibrium (\pi) distribution
            #return tm.background
            # Eliminate entirely
            return None
        return [ int( symbol == nuc ) for nuc in nucs ] 
    else:
        # Traverse children and determine likelihoods
        ## l_children = [ felsen( edge.tip, column, tm ) for edge in node.edges ]
        l_children = []
        for edge in node.edges:
            t = felsen( edge.tip, column, tm )
            if t is not None:
                l_children.append( ( t, edge ) )
        if not l_children: return None
        l_self = []
        # Determine liklihood for each possible 'ancestral' sequence
        for i_y in range( len( nucs ) ):
            y = nucs[ i_y ]
            p_y = 1
            # For each child, sum over all paths
            for cl, edge in l_children:
                if edge.length not in matrix_by_time_cache:
                    matrix_by_time_cache[ edge.length ] = matrix_for_time( tm.matrix, edge.length )
                psub = matrix_by_time_cache[ edge.length ]
                p_y_c = 0
                for i_a in range( len( nucs ) ):
                    p_y_c += psub[i_y,i_a] * cl[ i_a ]
                # Product over all children
                p_y *= p_y_c
            l_self.append( p_y )
        return l_self

def lik_to_prob( lik, p_parent, tm ):
    nucs = tm.alphabet
    p = []
    for i in range( len( nucs ) ):
        p.append( lik[i] * p_parent[i] )
    p_x = sum( p )
    # print >>sys.stderr, lik
    # print >>sys.stderr, p
    rval = [ top / p_x for top in p ]
    # print >>sys.stderr, rval
    return rval

def matrix_for_time( matrix, t ):
    return expm( matrix * t )
        
if __name__ == "__main__":
    main()
