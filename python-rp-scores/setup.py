from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

gill_code = [ 'pst.c', 'priority_queue.c', 'pfa.c', 'pattern_match.c', 'numrec.c', 'learn_algs.c', 'file_handle.c' ]
gill_code = [ 'third_party/gill_pst/' + l for l in gill_code ]

setup( name = "RP",
       packages = [ 'rp', 'rp.models' ],
       scripts = [ 'rp_train.py', 'rp_score.py', 'rp_score_maf.py', 'rp_score_maf_BUGGED.py',
                   'rp_score_bed_ints.py', 'rp_cv.py', 'rp_adapt.py', 
                   'rp_adapt_mpi.py', 'rp_adapt_mc.py', 'maf_to_ints.py',
                   'rp_prob_classify.py',
                   'rp_make_mapping.py' ],
       ext_modules=[ Extension( "rp.models.standard", ["rp/models/standard.pyx", "rp/models/standard_core.c"] ),
                     #Extension( "rp.models.gill_pst", ["rp/models/gill_pst.pyx"] + gill_code, include_dirs=['third_party/gill_pst/'] ), 
                     #Extension( "rp.models.tree", ["rp/models/tree.pyx"] ), 
                     #Extension( "rp.models.tree_pruned_1", ["rp/models/tree_pruned_1.pyx"] ), 
                     Extension( "rp.models.tree_pruned", ["rp/models/tree_pruned.pyx"] ), 
                     #Extension( "rp.models.simple_periodic", ["rp/models/simple_periodic.pyx", "rp/models/simple_periodic_core.c"] ),
                     #Extension( "rp.models.complex_periodic", ["rp/models/complex_periodic.pyx", "rp/models/simple_periodic_core.c"] ),
                    ],
       cmdclass = {'build_ext': build_ext} )
