seq_types = type( () ), type( [] )

def flatten( *args ):
    for arg in args:
        if type( arg ) in seq_types:
            for elem in arg:
                for f in flatten( elem ): 
                    yield f
        else:
            yield arg

def cross_lists(*sets):
    """Return the cross product of the arguments"""
    wheels = map(iter, sets) 
    digits = [it.next() for it in wheels]
    while True:
        yield digits[:]
        for i in range(len(digits)-1, -1, -1):
            try:
                digits[i] = wheels[i].next()
                break
            except StopIteration:
                wheels[i] = iter(sets[i])
                digits[i] = wheels[i].next()
        else:
            break

