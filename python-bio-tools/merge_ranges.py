#!/usr/bin/env python

import sys, ranges


range_reader = ranges.Reader( sys.stdin )

first = range_reader.next()

start = first.start
end = first.end

for range in range_reader:

        # Ensure sorted

        if range.start < start: raise "Input must be sorted"

        # If this range overlaps a previous one, just extend

        if range.start <= end:

            end = max( end, range.end )

        # Otherwise, right out the range

        else: 
        
            print start, end

            start = range.start
            end = range.end

if start > 0: print start, end
