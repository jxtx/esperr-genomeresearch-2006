from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup( name = "RP",
       packages = [ 'rp', 'rp.models' ],
       ext_modules=[ Extension( "rp.models.standard", ["rp/models/standard.pyx", "rp/models/standard_core.c"] ),
                     Extension( "rp.mapping_helper", ["rp/mapping_helper.pyx"] ) ],
       scripts = [ 'rp_train.py', 'rp_score.py' ],
       cmdclass = {'build_ext': build_ext} )
