#!/usr/bin/env python2.4

from __future__ import division

from pylab import *
from matplotlib.lines import Line2D

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
     
class PlotBuilder( object ):
    def __init__( self ):
        self.last_pos = 0
        self.labels = []
        self.lines = [] 
    def plot_subtree( self, node ):
        left, right, sim = node 
        # Left subtree   
        if type( left ) == int:
            left_level = 0
            left_pos = self.last_pos
            self.labels.append( left )
            self.last_pos += 1
        else:
            left_level, left_pos = self.plot_subtree( left )
        # Right subtree
        if type( right ) == int:
            right_level = 0
            right_pos = self.last_pos
            self.labels.append( right )
            self.last_pos += 1
        else:
            right_level, right_pos = self.plot_subtree( right )
        # Add lines
        # Hack
        sim = .001 - sim
        self.lines.append( ( left_pos, left_pos, sim, left_level ) )
        self.lines.append( ( right_pos, right_pos, sim, right_level ) )
        self.lines.append( ( left_pos, right_pos, sim, sim ) )
        return sim, ( ( left_pos + right_pos ) / 2 )
    def plot_to_axes( self, ax ):
        for left, right, top, bottom in self.lines:
            ax.add_line( Line2D( ( left, right), ( top, bottom ) ) )

def main():
    pb = PlotBuilder()
    root = read_dend( sys.stdin )
    pb.plot_subtree( root )
    ax = gca()
    pb.plot_to_axes( ax )
    xticks( range( 0, len( pb.labels ) ), [ str( l ) for l in pb.labels ] )
    ylim( 0, .001 )
    show()
    
if __name__ == "__main__": main()
        