import rp.models.standard as standard

def train( order, radix, pos_strings, neg_strings ):
    return standard.train( order, radix, pos_strings, neg_strings, 1 )

def from_file( f ):
    return standard.from_file( f )
