#!/usr/bin/env python2.3

from __future__ import division

from itertools import *
from tempfile import mktemp

import commands
import math
import os
import random
import string
import sys

class CVClassification( object ):
    def __init__( self ):
        self.neg = 0
        self.unc_neg = 0
        self.unc_pos = 0
        self.pos = 0
    def get_total( self ):
        return self.neg + self.unc_neg + self.unc_pos + self.pos
    def __str__( self ):
        return "%4d %4d %4d %4d" % ( self.pos, self.unc_pos, self.unc_neg, self.neg )

class CV( object ):

    def __init__( self, model_class, data1, data2, fold=5, passes=5 ):
        self.model_class = model_class
        self.data1 = data1
        self.data2 = data2
        self.fold = fold
        self.passes = passes

    def get_success_rate( self ):
        return ( float( self.cls1.unc_pos + self.cls1.pos + self.cls2.unc_neg + self.cls2.neg ) /
                 float( self.cls1.get_total() + self.cls2.get_total() ) )

    def run( self ):
        # Initialize classifications
        self.cls1 = CVClassification()
        self.cls2 = CVClassification()
        # Run everything 'passes' times
        for p in range( self.passes ):
            # Create random partitions
            partition1 = [ i % self.fold for i in range( len( self.data1 ) ) ]
            random.shuffle( partition1 )
            partition2 = [ i % self.fold for i in range( len( self.data2 ) ) ]
            random.shuffle( partition2 )
            # Run each fold
            for f in range( self.fold ):
                train1, test1 = self.split_by_partition( self.data1, partition1, f )
                train2, test2 = self.split_by_partition( self.data2, partition2, f )
                self.run_fold( train1, train2, test1, test2 )

    def run_fold( self, train_set_1, train_set_2, test_set_1, test_set_2 ):
        """Run one fold of the cross validation"""
        # Build predictor from training sets
        predictor = self.model_class( train_set_1, train_set_2 )
        # Score all sequences
        train_set_scores_1 = [ predictor.score( x ) for x in train_set_1 ]
        train_set_scores_2 = [ predictor.score( x ) for x in train_set_2 ]
        test_set_scores_1 = [ predictor.score( x ) for x in test_set_1 ]
        test_set_scores_2 = [ predictor.score( x ) for x in test_set_2 ]
        # Determine threshold
        low, mid, high = self.determine_threshold_simple( train_set_scores_1, train_set_scores_2 )
        # Classify
        self.classify( test_set_scores_1, low, mid, high, self.cls1 )
        self.classify( test_set_scores_2, low, mid, high, self.cls2 )

    def split_by_partition( self, set, partition, f ):
        train, test = [], []
        for i in range( len( set ) ):
            if partition[i] == f: test.append( set[i] )
            else: train.append( set[i] )
        return train, test

    def classify( self, scores, low, mid, high, cls ):
        for score in scores:
            if score < low: cls.neg += 1
            elif score < mid: cls.unc_neg += 1
            elif score < high: cls.unc_pos += 1
            else: cls.pos += 1

    def determine_threshold_simple( self, set1, set2 ):
        smallest_pos = min( set1 )
        largest_neg = max( set2 )
        # If completely separated
        if smallest_pos > largest_neg:
            high = smallest_pos
            low = largest_neg + 0.00000000001
            mid = 0
        # Else overlap
        else:
            high = low = mid = 0
        # Return the thresholds
        return low, mid, high

    def determine_threshold( self, set1, set2 ):
        sorted1 = set1[:]; sorted1.sort()
        sorted2 = set2[:]; sorted2.sort()
        # If completely separated
        if sorted1[0] > sorted2[-1]:
            high = sorted1[0]
            low = sorted2[-1] + 0.00000000001
            mid = ( high + low ) / 2.0
        # Else overlap
        else:
            count1, count2 = len( set1 ), len( set2 )
            index1, index2 = 0, 0
            best_qual, best_score = 0.0, 0.0
            while 1:
                current_qual = ( ( float( count1 - index1 ) / float( count1 ) )
                               + ( float( index2 ) / float( count2 ) ) )
                if index2 < count2 and ( ( index1 == count1 ) or ( sorted2[ index2 ] < sorted1[ index1 ] ) ):
                    current_score = sorted2[ index2 ]
                    index2 += 1
                elif index1 < count1:
                    current_score = sorted1[ index1 ]
                    index1 += 1
                else:
                    break
                if current_qual > best_qual:
                    best_score = current_score
                    best_qual = current_qual
            high = low = mid = best_score
        # Return the thresholds
        return low, mid, high
