from __future__ import absolute_import
from builtins import object
from full_physics.docopt import docopt
import re

class DocOptSimple(object):
    '''The package docopt (http://docopt.org) is a nice package,
    but it has the disadvantage that getting the options etc. from it
    uses a somewhat unnatural interface (e.g., --my-opt=n comes as
    int(d["--my-opt"]) rather than OptionParser cleaner sort of
    options.my_opt). This class tries to give a simple interface that
    is sufficient for many purposes. If you try to access the attribute
    "my_val", we first look for the argument "my_val", then "<my_val>",
    then "--my_val", then "--my-val". If the value looks like a integer,
    we return it to an integer. If it looks like a float, we return it to
    a float. If we don't otherwise recognize it, we return this as a string.'''
    def __init__(self, doc, argv=None, help=True, version=None, 
                 options_first=False):
        self.args = docopt(doc, argv=argv, help=help, version=version, 
                           options_first=options_first)

    def __contains__(self, name):
        for key in (name, 
                    "<" + name + ">", 
                    "--" + name,
                    "--" + name.replace("_", "-")):
            if(key in self.args):
                return True
        return False

    def __getattr__(self, name):
        for key in (name, 
                    "<" + name + ">", 
                    "--" + name,
                    "--" + name.replace("_", "-")):
            if(key in self.args):
                return self.__find_type(key)
        raise AttributeError(name)

    def __find_type(self, key):
        '''Find the type of the value, and return in'''
        v = self.args[key]
        if(isinstance(v, str)):
            if(re.match('[+-]?\d+$', v)):
                return int(v)
            if(re.match('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$', v)):
                return float(v)
            if(re.match('[-+]?[0-9]+\.([eE][-+]?[0-9]+)?$', v)):
                return float(v)
        return v

def docopt_simple(doc, argv=None, help=True, version=None, 
                 options_first=False):
    '''The package docopt (http://docopt.org) is a nice package,
    but it has the disadvantage that getting the options etc. from it
    uses a somewhat unnatural interface. This gives a simpler interface
    sufficient for many programs. See DocOptSimple class for details
    on the interface.'''
    return DocOptSimple(doc, argv=argv, help=help, version=version, 
                        options_first=options_first)

