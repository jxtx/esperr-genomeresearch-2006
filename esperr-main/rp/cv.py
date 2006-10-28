#!/usr/bin/env python2.3

from __future__ import division

from numpy import argmax

from itertools import *
from tempfile import mktemp

import commands
import math
import os
import random
import string
import sys


class MultiCV( object ):
    def __init__( self, model_class, data, fold=5, passes=5, loo=False ):
        self.model_class = model_class
        self.data = data
        self.ndata = len( data )
        self.fold = fold
        self.passes = passes
        self.loo = loo
    def get_success_rate( self ):
        total = 0
        for c in self.cls:
            total += ( c.pos / c.get_total() )
        return total / len( self.cls )
    def get_summary( self ):
        return ", ".join( [ "%d / %d" % ( self.cls[i].pos, self.cls[i].get_total() ) for i in range( self.ndata ) ] )
    def run( self ):
        if self.loo: 
            self.run_loo()
        else: 
            self.run_folds()
    def run_folds( self ):
        # Initialize classifications
        self.cls = [ CVClassification() for i in range( self.ndata ) ]
        # Run everything 'passes' times
        for p in range( self.passes ):
            # Create random partitions 
            partitions = []
            for i, d in enumerate( self.data ):
                partition = [ i % self.fold for i in range( len( d ) ) ]
                random.shuffle( partition )
                partitions.append( partition )
            # Run each fold
            for f in range( self.fold ):
                train_sets = []
                test_sets = []
                for partition, data in izip( partitions, self.data ):
                    train, test = split_by_partition( data, partition, f )
                    train_sets.append( train )
                    test_sets.append( test )
                self.run_fold( train_sets, test_sets )
    def run_loo( self ):
        # Initialize classifications
        self.cls = [ CVClassification() for i in range( self.ndata ) ]
        # Run everything 'passes' times
        for p in range( self.passes ):
            # Run for each item in positive set
            for i, d in enumerate( self.data ):
                for j in range( len( d ) ):                    
                    train_sets = self.data[:]
                    train_sets[i] = self.data[i][:]
                    test_sets = [ [] for _ in range( len( self.data ) ) ]
                    test_sets[i].append( train_sets[i][j] )
                    del train_sets[i][j]
                    self.run_fold( train_sets, test_sets )
    def run_fold( self, train_sets, test_sets ):
        """Run one fold of the cross validation"""
        # Build predictor from training sets
        predictors = []
        for train in train_sets:
            predictors.append( self.model_class( train ) )
        # Score each sequence 
        for i, train in enumerate( test_sets ):
            for s in train:
                scores = [ p.score( s ) for p in predictors ]
                c = argmax( scores )
                if c == i:
                    self.cls[i].pos += 1
                else:
                    self.cls[i].neg += 1

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

    def __init__( self, model_class, data1, data2, fold=5, passes=5, loo=False ):
        self.model_class = model_class
        self.data1 = data1
        self.data2 = data2
        self.fold = fold
        self.passes = passes
        self.loo = loo

    def get_success_rate( self ):
        return ( float( self.cls1.unc_pos + self.cls1.pos + self.cls2.unc_neg + self.cls2.neg ) /
                 float( self.cls1.get_total() + self.cls2.get_total() ) )

    def run( self ):
        if self.loo: self.run_loo()
        else: self.run_folds()

    def run_folds( self ):
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

    def run_loo( self ):
        # Initialize classifications
        self.cls1 = CVClassification()
        self.cls2 = CVClassification()
        # Run everything 'passes' times
        for p in range( self.passes ):
            # Run for each item in positive set
            for i in range( len( self.data1 ) ):
                test = [ self.data1[i] ]
                train = list( self.data1 )
                del train[i]
                self.run_fold( train, self.data2, test, [] )
            # And each item in negative set
            for i in range( len( self.data2 ) ):
                test = [ self.data2[i] ]
                train = list( self.data2 )
                del train[i]
                self.run_fold( self.data1, train, [], test )

    def run_fold( self, train_set_1, train_set_2, test_set_1, test_set_2 ):
        """Run one fold of the cross validation"""
        # Build predictor from training sets
        predictor = self.model_class( train_set_1, train_set_2 )
        # Score all sequences
        train_set_scores_1 = [ predictor.score( x ) for x in train_set_1 ]
        train_set_scores_2 = [ predictor.score( x ) for x in train_set_2 ]
        test_set_scores_1 = [ predictor.score( x ) for x in test_set_1 ]
        test_set_scores_2 = [ predictor.score( x ) for x in test_set_2 ]
        # Free model
        del predictor
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

def split_by_partition( set, partition, f ):
    train, test = [], []
    for i in range( len( set ) ):
        if partition[i] == f: test.append( set[i] )
        else: train.append( set[i] )
    return train, test