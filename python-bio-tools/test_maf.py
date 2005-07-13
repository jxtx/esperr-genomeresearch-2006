import unittest
from bx import align
from bx import align.maf as maf

from StringIO import StringIO

class MAFTestCase( unittest.TestCase ):

    def testReader( self ):

        reader = maf.Reader( file( "test.maf" ) )
        assert reader.attributes["version"] == "1" 
        assert reader.attributes["scoring"] == "humor.v4" 

        a = reader.next()
        assert a.score == "0.128"
        assert len( a.components ) == 3
        check_component( a.components[0], "human_hoxa", 100, 8,  "+", 100257, "ACA-TTACT" )
        check_component( a.components[1], "horse_hoxa", 120, 9, "-",  98892, "ACAATTGCT" )
        check_component( a.components[2], "fugu_hoxa",    88, 7,  "+",  90788, "ACA--TGCT" )
        
        a = reader.next()
        assert a.score == "0.071"
        assert len( a.components ) == 3
        check_component( a.components[0], "human_unc", 9077, 8, "+", 10998, "ACAGTATT" )
        check_component( a.components[1], "horse_unc", 4555, 6, "-",  5099, "ACA--ATT" )
        check_component( a.components[2], "fugu_unc",   4000, 4, "+",  4038, "AC----TT" )

        a = reader.next()
        assert a is None

        reader.close()

    def testWriter( self ):

        val = StringIO()
        writer = maf.Writer( val, { 'scoring':'foobar' } )
        
        a = align.Alignment()
        a.score = "7009"

        a.components.append( align.Component( src="human_hoxa", start=100, size=9,  strand="+", src_size=1000257, text="ACA-TTACT" ) )
        a.components.append( align.Component( src="horse_hoxa", start=120, size=10, strand="-",   src_size=98892, text="ACAATTGCT" ) )

        check_component( a.components[0], "human_hoxa", 100, 9,  "+", 1000257, "ACA-TTACT" )
        check_component( a.components[1], "horse_hoxa", 120, 10, "-",   98892, "ACAATTGCT" )

        writer.write( a )

        assert val.getvalue() == """##maf version=1 scoring=foobar
a score=7009
s human_hoxa 100  9 + 1000257 ACA-TTACT 
s horse_hoxa 120 10 -   98892 ACAATTGCT 

"""

    def test_slice( self ):

        a = align.Alignment()
        a.score = "7009"
        a.components.append( align.Component( src="human_hoxa", start=100, size=9,  strand="+", src_size=100257, text="ACA-TTACT" ) )
        a.components.append( align.Component( src="horse_hoxa", start=120, size=10, strand="-",   src_size=98892, text="ACAATTGCT" ) )
    
        b = a.slice_by_component( 0, 101, 105 )

        check_component( b.components[0], src="human_hoxa", start=101, size=4, strand="+", src_size=100257, text="CA-TT" )
        check_component( b.components[1], src="horse_hoxa", start=121, size=5, strand="-", src_size=98892, text ="CAATT" )

def check_component( c, src, start, size, strand, src_size, text ):
    assert c.src == src
    assert c.start == start 
    assert c.size == size 
    assert c.strand == strand 
    assert c.src_size == src_size 
    assert c.text == text

if __name__ == "__main__":
    unittest.main()
