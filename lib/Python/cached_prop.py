def cached_property(prop):
    '''This is a simple decorator that calculates the value for a function
    taking only the self argument once, and then return the results.'''
    # Note the use of an array here rather than just saying cache = None is
    # because of a subtle problem with closures in python. You can read the
    # details at http://stackoverflow.com/questions/141642/what-limitations-have-closures-in-python-compared-to-language-x-closures,
    # but basically any outer scope variables need to be mutable containers
    cache = [None]
    def wrapper(self):
        if(cache[0] is None):
            cache[0] = prop(self)
        return cache[0]

    return property(wrapper, None, None, prop.__doc__)
