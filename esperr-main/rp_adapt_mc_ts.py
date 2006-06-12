#!/usr/bin/env python

"""
Collapse from a starting alphabet to a series of nested alphabets using
cross validation as 'merit'.

usage: %prog pos_data neg_data out_dir [options]
   -a, --atoms=FILE:   A mapping specifying the largest set of symbols (these never get broken)
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence
"""

from __future__ import division

import bx.align.maf

import cookbook.doc_optparse
import os.path
import random
import sys
import time
import traceback

from Numeric import *
from cookbook.progress_bar import *
from itertools import *
from operator import *

import rp.cv
import rp.io
import rp.models
import rp.mapping

import rp.models.complex_periodic

stop_size = 5
fold = 5
passes = 5

samp_size_collapse = 30
samp_size_expand = 10

# This should be in a module somewhere
def all_pairs( n ):
    rval = []
    for i in range( 0, n ):
        for j in range( i + 1, n ):
            rval.append( ( i, j ) )
    return rval

# A mapping for the alphabet in which twinscan style alignment rows are 
# spelled (specifically we add a symbol for 'not aligned')
TS_DNA = rp.mapping.CharToIntArrayMapping()
TS_DNA.set_mapping( "a", 0 )
TS_DNA.set_mapping( "A", 0 )
TS_DNA.set_mapping( "c", 1 )
TS_DNA.set_mapping( "C", 1 )
TS_DNA.set_mapping( "g", 2 )
TS_DNA.set_mapping( "G", 2 )
TS_DNA.set_mapping( "t", 3 )
TS_DNA.set_mapping( "T", 3 )
TS_DNA.set_mapping( "-", 4 )
TS_DNA.set_mapping( "*", 5 )

TS_DNA_BASE = TS_DNA.get_out_size()

# PERIOD=1

# class ProductModel( object ):
#     def __init__( self, radix, pos_strings, neg_strings ):
#         self.genome_model = rp.models.simple_periodic.train( 5, 4, imap( itemgetter(0), pos_strings ), imap( itemgetter(0), neg_strings ), pos_period=PERIOD )
#         self.align_model = rp.models.simple_periodic.train( 1, radix, imap( itemgetter(1), pos_strings ), imap( itemgetter(1), neg_strings ), pos_period=PERIOD )
#     def score( self, string ):
#         #return self.genome_model.score( string[0] ) + self.align_model.score( string[1] )
#         return self.align_model.score( string[1] )

class Main( object ):

    def read_maf( self, fname ):
        all_texts = []
        for block in bx.align.maf.Reader( open( fname ) ):
            assert len( block.components ) - 1 == self.align_count, \
                "Alignments should all have %d rows" % self.align_count + 1
            # Build human sequence and alignment column sequence, dropping
            # any columns that are gaps in human
            texts = [ [] for i in range( self.align_count + 1 ) ]
            for i in range( len( block.components[0].text ) ):
                if block.components[0].text[i] != '-':
                    for j in range( self.align_count + 1 ):
                        texts[j].append( block.components[j].text[i] )
            all_texts.append( texts )
            #print >>sys.stderr, block, texts
        genome_seqs = []
        align_seqs = []
        for texts in all_texts:
            genome_seqs.append( TS_DNA.translate( ''.join( texts[0] ) ) )
            align_seqs.append( self.atom_mapping.translate( TS_DNA.translate_list( [ ''.join( text ) for text in texts[1:] ] ) ) )
        return genome_seqs, align_seqs
            
    def run( self ):

        mapping = self.starting_mapping

        # Open merit output
        merit_out = open( os.path.join( self.out_dir, 'merits.txt' ), 'w' )

        # Read training data and build seperate sets for sequence (translated 
        # to ints in range(0,4)) and alignment (translated with self.atom_mapping)
        print >>sys.stderr, "Loading training data"
        self.pos_genome_seqs, self.pos_strings = self.read_maf( self.pos_fname )
        self.neg_genome_seqs, self.neg_strings = self.read_maf( self.neg_fname )
        
        # Count how many times each atom appears in the training data
        atom_counts = zeros( self.atom_mapping.get_out_size() )
        for string in chain( self.pos_strings, self.neg_strings ):
            for val in string:
                atom_counts[ val ] += 1

        # Valid candiates for expansion must occur more than 10 times in the training data
        can_expand = compress( atom_counts > 10, arange( len( atom_counts ) ) )

        # Handling bad columns in the training data is not obvious, so don't do it for now
        # for string in chain( pos_strings, neg_strings ):
        #    assert -1 not in string, "Cannot have invalid columns (map to -1) in training data"

        best_merit_overall = 0
        best_mapping_overall = None
        index_best_merit_overall = 0
        out_counter = 0

        step_counter = 0

        # 
        last_force_counter = 0

        print >>sys.stderr, "Searching"

        # Collapse
        while 1:

            symbol_count = mapping.get_out_size()


            best_merit = 0
            best_mapping = None

            clock = time.clock()
            cv_runs = 0

            # First try a bunch of collapses
            if symbol_count > stop_size:
                pairs = all_pairs( symbol_count )
                if len( pairs ) > samp_size_collapse: pairs = random.sample( pairs, samp_size_collapse ) 
                for i, j in pairs:
                    new_mapping = mapping.collapse( i, j )
                    merit = self.calc_merit( new_mapping )
                    cv_runs += 1
                    if merit > best_merit:
                        best_merit = merit
                        best_mapping = new_mapping

            # Also try a bunch of expansions
            elements = random.sample( can_expand, samp_size_expand )
            for i in elements:
                new_mapping = mapping.expand( i )
                if new_mapping.get_out_size() == symbol_count: continue
                merit = self.calc_merit( new_mapping )
                cv_runs += 1
                if merit > best_merit:
                    best_merit = merit
                    best_mapping = new_mapping

            clock = time.clock() - clock

            mapping = best_mapping

            # Append merit to merit output        
            print >>merit_out, step_counter, symbol_count, best_merit
            merit_out.flush()

            if best_merit >= best_merit_overall:
                best_merit_overall = best_merit
                best_mapping_overall = best_mapping
                # So we know what step the best mapping was encountered at
                best_merit_overall_index = step_counter
                restart_counter = step_counter
                # Reset the counter we use to force expansions
                last_force_counter = step_counter
                # Write best mapping to a file
                mapping_out = open( os.path.join( self.out_dir, "%03d.mapping" % out_counter ), 'w' )
                for i, symbol in enumerate( self.atom_mapping.get_table() ): 
                    # Apply the 'second' mapping to the atom symbol
                    if symbol >= 0: symbol = mapping[ symbol ]
                    print >>mapping_out, str.join( '', TS_DNA.reverse_map( i, self.align_count ) ), symbol
                mapping_out.close()
                out_counter += 1

            print >>sys.stderr, "%06d, New best merit: %2.2f%%, size: %d, overall best: %2.2f%% at %06d, cvs per sec: %f" \
                % ( step_counter, best_merit * 100, mapping.get_out_size(), best_merit_overall * 100, best_merit_overall_index, cv_runs/clock  )

            # If we have gone 50 steps without improving over the best, restart from best
            if step_counter > restart_counter + 50:
                print >>sys.stderr, "Restarting from best mapping"
                print >>merit_out, step_counter, "RESTART"
                mapping = best_mapping_overall
                restart_counter = step_counter
                # Immediately force expansions after restart
                last_force_counter = 0

            if step_counter > last_force_counter + 20:
                last_force_counter = step_counter
                print >>sys.stderr, "Forcing expansions"
                print >>merit_out, step_counter, "FORCED EXPANSIONS"
                for i in range( 5 ):
                    symbol_count = mapping.get_out_size()
                    best_merit = 0
                    best_mapping = None
                    for i in random.sample( can_expand, samp_size_expand ):
                        new_mapping = mapping.expand( i )
                        if new_mapping.get_out_size() == symbol_count: continue
                        merit = self.calc_merit( new_mapping )
                        if merit > best_merit:
                            best_merit = merit
                            best_mapping = new_mapping
                    mapping = best_mapping
            step_counter += 1

    def calc_merit( self, mapping ):
        # Apply mapping to strings
        pos_strings = [ mapping.translate( s ) for s in self.pos_strings ]
        neg_strings = [ mapping.translate( s ) for s in self.neg_strings ]
        # Cross validate using those strings
        radix = mapping.get_out_size()
        ## model_factory = lambda d0, d1: ProductModel( radix, d0, d1 )
        model_factory = lambda d0, d1: rp.models.complex_periodic.train( 5, 1, 4, radix, d0, d1 )
        cv_engine = rp.cv.CV( model_factory, 
                              zip( self.pos_genome_seqs, pos_strings ), 
                              zip( self.neg_genome_seqs, neg_strings ), 
                              fold=fold, passes=passes )
        cv_engine.run()
        # Merit is TP + TN
        ## print "Pos:", cv_engine.cls1
        ## print "Neg:", cv_engine.cls2
        return ( cv_engine.cls1.pos / ( len( pos_strings ) * passes ) + cv_engine.cls2.neg / ( len( neg_strings ) * passes ) ) / 2

    def main( self ):

        # Parse command line
        options, args = cookbook.doc_optparse.parse( __doc__ )

        if 1:
            self.pos_fname, self.neg_fname, self.out_dir = args
            # Load the mapping that specifies the atoms (we never break these apart)
            self.align_count, self.atom_mapping = rp.mapping.alignment_mapping_from_file( file( options.atoms ), char_mapping=TS_DNA )
            # Load a second mapping if provided (must be a partition of the atom mapping)
            if options.mapping: 
                self.starting_mapping = rp.mapping.second_mapping_from_file( file( options.mapping ), self.atom_mapping )
            else: 
                self.starting_mapping = rp.mapping.identity_mapping( self.atom_mapping.get_out_size() )
        
        #except:
        #    cookbook.doc_optparse.exit()


        self.run()


if __name__ == "__main__":
    Main().main()
