import unittest

from intervals import *

test_intervals = [ 
    Interval( 0,  10,    "1" ),
    Interval( 9,   40,   "2" ),
    Interval( 5,   11,   "3" ),
    Interval( 30,  42,   "4" ),
    Interval( 31,  43,   "5" ),
    Interval( 11,  1000, "6" ),
    Interval( 100, 101,  "7" ),
    Interval( 20,  4000, "8" )
]

class TestCase( unittest.TestCase ):

    def testBasic( self ):
        nx = Intersecter()
        for i in test_intervals: nx.add_interval( i )

        r = nx.find( 5, 35 )
        r.sort()

        assert [ '1', '3', '2', '6', '8', '4', '5' ] == [ interval.value for interval in r ]

if __name__ == "__main__":
    unittest.main()
