from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

setup( name = "RP",
       packages = [ 'rp', 'rp.models' ],
       scripts = [ 'rp_train.py', 'rp_score.py', 'rp_score_maf.py', 'rp_cv.py', 'rp_adapt.py', 'rp_adapt_mpi.py' ],
       ext_modules=[ Extension( "rp.models.standard", ["rp/models/standard.pyx", "rp/models/standard_core.c"] ),
                     Extension( "rp.mapping_helper", ["rp/mapping_helper.pyx"] ) ],
       cmdclass = {'build_ext': build_ext} )
