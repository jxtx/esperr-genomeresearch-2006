import rp_train
import unittest
import rp.mapping

from StringIO import StringIO

class TestCase( unittest.TestCase ):

    def test_train_from_maf( self ):
        
        out = StringIO()
        mapping = rp.mapping.alignment_mapping_from_file( open( "test_data/hm.5a.mapping.txt" ) )
        
        rp_train.run( open( "test_data/reg.maf" ), 
                      open( "test_data/ar.maf" ), 
                      out, 
                      'maf', 
                      mapping, 
                      None, 
                      2 )
                      
        expected = open( "test_data/compare/hm.5a_scoreMatrix.txt" )
        
        compare_sm( out, expected )
        
    def test_train_from_ints( self ):
        
        out = StringIO()
        
        rp_train.run( open( "test_data/reg.ints" ), 
                      open( "test_data/ar.ints" ), 
                      out, 
                      None, 
                      None, 
                      None, 
                      2 )
                      
        expected = open( "test_data/compare/hm.5a_scoreMatrix.txt" )
        
        compare_sm( out, expected )
        
def compare_sm( a, b ):
    for line1, line2 in zip( a, b ):
        assert round( float( a ), 9 ) == round( float( b ), 9 )
        
if __name__ == "__main__": unittest.main()