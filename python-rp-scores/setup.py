from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup( name = "RP",
       packages = [ 'rp' ],
       ext_modules=[ Extension( "rp.standard_model", ["rp/standard_model.pyx", "rp/standard_model_helper.c"] ),
                     Extension( "rp.mapping_helper", ["rp/mapping_helper.pyx"] ) ],
       scripts = [ 'rp_train.py', 'rp_score.py' ],
       cmdclass = {'build_ext': build_ext} )
