from distutils.core import setup

scripts = open( "scripts.list" ).read().split()
all_packages = open( "packages.list" ).read().split()
py_packages = [ p[:-3] for p in all_packages if p.endswith( ".py" ) ]
print py_packages
packages = [ p for p in all_packages if not p.endswith( ".py" ) ]

setup( 
        name = "python-bio-tools",
        py_modules = py_packages,
        packages = packages,
        scripts = open( "scripts.list" ).read().split(),
     )
