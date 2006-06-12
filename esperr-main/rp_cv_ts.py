#!/usr/bin/env python

"""
Use cross validation to evaluate a model for some training data.

usage: %prog pos_data neg_data [options]
   -f, --format=FILE:  Format of input data. 'ints' by default, or 'maf'
   -m, --mapping=FILE: A mapping (alphabet reduction) to apply to each sequence (optional)
   -r, --radix=N:      Radix (optional)
   -o, --orders=N,...: Orders to cross validate over
   -F, --fold=N:       Fold (default 5)
   -M, --model=name:   Name of model to train (default 'standard')
   -l, --loo:          Use leave-one-out cross validation (fold is ignored in this case) 
"""

import bx.align.maf

import array
import cookbook.doc_optparse
import sys
import traceback
import time

import rp.cv
import rp.io
import rp.mapping
import rp.models.complex_periodic

from itertools import *
from operator import *

default_fold = 5

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

# PERIOD=2
#
# class ProductModel( object ):
#     def __init__( self, radix, pos_strings, neg_strings ):
#         self.genome_model = rp.models.simple_periodic.train( 5, 4, imap( itemgetter(0), pos_strings ), imap( itemgetter(0), neg_strings ), pos_period=PERIOD )
#         self.align_model = rp.models.simple_periodic.train( 1, radix, imap( itemgetter(1), pos_strings ), imap( itemgetter(1), neg_strings ), pos_period=PERIOD )
#     def score( self, string ):
#         return self.genome_model.score( string[0] ) + self.align_model.score( string[1] )
#         #return self.genome_model.score( string[0] )
#         #return self.align_model.score( string[1] )

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
            align_seqs.append( self.mapping.translate( TS_DNA.translate_list( [ ''.join( text ) for text in texts[1:] ] ) ) )
        return genome_seqs, align_seqs

    def run( self ):

        # Split up 

        # Read integer sequences
        self.pos_genome_seqs, self.pos_strings = self.read_maf( self.pos_fname )
        self.neg_genome_seqs, self.neg_strings = self.read_maf( self.neg_fname )

        # Determine radix
        radix = self.mapping.get_out_size()

        print "  TP  ~TP  ~FN   FN   FP  ~FP  ~TN   TN       %    time"

        # Cross validate 
        model_factory = lambda d0, d1: rp.models.complex_periodic.train( 5, 1, 4, radix, d0, d1 )
        if self.loo: passes = 1
        else: passes = 5 
        cv_engine = rp.cv.CV( model_factory, 
                      zip( self.pos_genome_seqs, self.pos_strings ), 
                      zip( self.neg_genome_seqs, self.neg_strings ), 
                      fold=self.fold, passes=passes )
        start_time = time.time()
        cv_engine.run()
        seconds = time.time() - start_time

        print cv_engine.cls1, cv_engine.cls2,
        print "  %2.2f    %2.2f" % ( cv_engine.get_success_rate()*100, seconds )
    
    def main( self ):

        # Parse command line

        options, args = cookbook.doc_optparse.parse( __doc__ )

        #try:
        if 1:
            self.pos_fname, self.neg_fname = args
            if options.fold:
                self.fold = int( options.fold )
            else:
                self.fold = default_fold
            self.align_count, self.mapping = rp.mapping.alignment_mapping_from_file( file( options.mapping ), char_mapping=TS_DNA )
            self.loo = bool( options.loo )    
        #except:
        #    cookbook.doc_optparse.exit()

        self.run()


if __name__ == "__main__": Main().main()
