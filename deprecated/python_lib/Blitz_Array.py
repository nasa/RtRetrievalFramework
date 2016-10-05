# -----------------------------------------------------------------------------
# Parser for Blitz Array Files
# -----------------------------------------------------------------------------

import os
import re
import sys
from types import ListType

import numpy

from Ply_Parser import Ply_Parser

class Full_Parser(Ply_Parser):
    """This parser is slow and probably should only be used as reference"""

    def __init__(self, filename=None, **kw):
        'Initialize the class'
        Ply_Parser.__init__(self, filename=filename, **kw)

    ##############
    # Lexing rules

    tokens = (
        'COMMENT', 'INTEGER', 'FLOAT', 'LPAREN', 'RPAREN', 'X'
        )

    # Tokens

    # A string containing ignored characters (spaces and tabs)
    t_ignore  = ' \t\n'

    t_LPAREN = r'\['
    t_RPAREN = r'\]'
    t_X      = r'[x]'

    def t_COMMENT(self, t):
        r'\#.*'
        return t

    def t_FLOAT(self, t):
        r'[-]?((\d*\.\d+)([EeDd][\+-]?\d+)?|([1-9]\d*[EeDd][\+-]?\d+)|(\d+\.))'
        return t

    def t_INTEGER(self, t):
        r'\d+'
        return t

    def t_error(self, t):
        next_space_loc = t.value.find(' ')
        raise ValueError("Illegal value '%s'... at line %d." % (t.value[0:next_space_loc+1], t.lexer.lineno))
        t.lexer.skip(1)

    ###############
    # Parsing rules

    def p_file_root(self, p):
        'file_root : file_contents'
        p[0] = tuple( p[1] )

    def p_file_contents_comments(self, p):
        'file_contents : COMMENT file_contents'
        p[0] = p[2]

    def p_file_contents_array(self, p):
        'file_contents : array file_contents'

        # Insert backwards since arrays are parsed backwards
        p[2].insert( 0, p[1] )
        p[0] = p[2]

    def p_file_contents_empty(self, p):
        'file_contents : empty'
        p[0] = []

    def p_array(self, p):
        'array : dimensions LPAREN value_list RPAREN'

        p[0] = numpy.array(p[3], dtype=float)
        p[0].shape = tuple(p[1])

    def p_dimensions(self, p):
        """dimensions : dimensions X INTEGER
                      | INTEGER
        """
        if len(p) > 2:
            p[1].append( int(p[3]) )
            p[0] = p[1]
        else:
            p[0] = [ int(p[1]) ]

    def p_value_list(self, p):
        """value_list : value_list value
                      | empty
        """
        if len(p) > 2:
            p[1].append( float(p[2]) )
            p[0] = p[1]
        else:
            p[0] = []

    def p_value(self, p):
        """
        value : INTEGER
              | FLOAT
        """
        p[0] = p[1]

    def p_empty(self, p):
        'empty :'
        pass

    def p_error(self, p):
        if p == None:
            raise ValueError("Syntax error with no info!")
        else:
            print dir(p)
            raise ValueError("Syntax error at token '%s' with value '%s' at line number %d of file %s" % (p.type, p.value, p.lineno, self.filename))

class Fast_Parser:

    def __init__(self, filename=None):
        self.filename = filename
        self.data = None
        pass

    def parse(self, filename=None):
        """Parse file and return data"""
        
        if filename != None: self.filename = filename
        if self.filename == None:
            raise Exception("File name must be specified for parsing and net set to None")

        # Use a list initially as there might be multiple definitions internally
        self.data = []

        # Read looking for [ ] marks and dim1 x dim2 ...
        # dim marks probably always fall on a seperate line,
        # but [ ] might be before or after data
        with open(self.filename) as file_obj:
            curr_data_obj = None
            # Mode flags
            in_data = False
            curr_dims = None
            curr_data = []
            for line_idx, line_txt in enumerate(file_obj):
                if not in_data:
                    if re.search('^\s*#', line_txt) or re.search('^\s+$', line_txt):
                        # comment or empty text
                        pass
                    elif line_txt.find('x') >= 0 or re.search('^\d+$', line_txt):
                        try:
                            # If something existing in curr dims then we parsed a number not
                            # a dimension
                            if curr_dims != None:
                                self.data.append(curr_dims)
                            curr_dims = [int(d) for d in re.split('\s*[x]\s*', line_txt)]
                        except ValueError:
                            raise ValueError("At line %d, error parsing dimension definition: %s" % (line_idx+1,line_txt))
                    elif line_txt.find('[') >= 0:
                        in_data = True
                        line_txt = line_txt.replace('[', '')
                    else:
                        # Single number
                        curr_dims = None
                        try:
                            self.data.append( float(line_txt) )
                        except ValueError:
                            # Everything else is just strings
                            self.data.append( line_txt )

                # Fall through if data begins on current line
                if in_data:
                    at_end = False
                    if line_txt.find(']') >= 0:
                        at_end = True
                        line_txt = line_txt.replace(']', '')
                        
                    fields = line_txt.strip().split()
                    curr_data += [float(f) for f in fields]

                    if at_end:
                        n_data = numpy.array(curr_data)
                        self.data.append( n_data.reshape(curr_dims) )

                        # Reset state variables
                        in_data = False
                        curr_dims = None
                        curr_data = []

        # If only one array then just return that single one
        # otherwise just change to a tuple so it is immutable
        if len(self.data) == 1:
            self.data = self.data[0]
        else:
            self.data = tuple(self.data)

        return self.data

    def data(self):
        """Return previously parsed data"""
        if self.data == None:
            self.parse()
        return self.data

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_file>"
        sys.exit(1)

    input_file  = sys.argv[1]

    parsed_values = Full_Parser(input_file).parse()
    print parsed_values
