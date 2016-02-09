from builtins import object
import os

import ply.lex as lex
import ply.yacc as yacc

class PlyParser(object):
    """
    Base class for a lexer/parser that has the rules defined as methods
    """
    tokens = ()
    precedence = ()

    def __init__(self, **kw):
        self.debug = kw.get('debug', 0)
        self.filename = kw.get('filename', None)
        self.names = { }

        try:
            modname = os.path.split(os.path.splitext(__file__)[0])[1] + "_" + self.__class__.__name__
        except:
            modname = "parser"+"_"+self.__class__.__name__
        self.debugfile = modname + ".dbg"
        self.tabmodule = modname + "_" + "parsetab"
        self.outputdir = os.path.dirname(__file__)

        # Build the lexer and parser
        lex.lex(module=self, debug=self.debug)
        # We can't assume that we can write to this directory (e.g., we've
        # been installed as a module). So go ahead and skip writing out. Our
        # grammers are simple enough that this shouldn't be much of a problem,
        # we just regenerate the grammer each time.
        #yacc.yacc(module=self,
        #          debug=self.debug,
        #          debugfile=self.debugfile,
        #          tabmodule=self.tabmodule,
        #          outputdir=self.outputdir,
        #          method='SLR')
        yacc.yacc(module=self, debug=self.debug, write_tables=0, method='SLR')
        

    def parse(self, text=None):

        if self.filename != None and text == None:
            with open(self.filename, 'r') as input_fo:
                file_contents = input_fo.read()
            return yacc.parse(file_contents)
        elif text != None:
            return yacc.parse(text)
        else:
            return None
