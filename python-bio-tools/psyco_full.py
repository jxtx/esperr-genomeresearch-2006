import sys

try:
    import psyco
    psyco.full()
except:
    print >> sys.stderr, "Psyco not found, continuing without it"


