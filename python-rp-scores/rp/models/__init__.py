def get( name ):
    try:
        return getattr( __import__( "rp.models", globals(), locals(), [ name ] ), name )
    except:
        raise "Unknown model: '%s'" % name