from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

gill_code = [ 'pst.c', 'priority_queue.c', 'pfa.c', 'pattern_match.c', 'numrec.c', 'learn_algs.c', 'file_handle.c' ]
gill_code = [ 'third_party/gill_pst/' + l for l in gill_code ]

setup( name = "RP",
       packages = [ 'rp', 'rp.models' ],
       scripts = [ 'rp_train.py', 'rp_score.py', 'rp_score_maf.py', 'rp_cv.py', 'rp_adapt.py', 'rp_adapt_mpi.py' ],
       ext_modules=[ Extension( "rp.models.standard", ["rp/models/standard.pyx", "rp/models/standard_core.c"] ),
                     Extension( "rp.models.gill_pst", ["rp/models/gill_pst.pyx"] + gill_code, include_dirs=['third_party/gill_pst/'] ), 
                     Extension( "rp.models.tree", ["rp/models/tree.pyx"] ), 
                     Extension( "rp.models.tree_pruned_1", ["rp/models/tree_pruned_1.pyx"] ), 
                     Extension( "rp.mapping_helper", ["rp/mapping_helper.pyx"] ) ],
       cmdclass = {'build_ext': build_ext} )
