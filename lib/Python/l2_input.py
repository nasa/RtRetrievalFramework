from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import six
import sys

# -----------------------------------------------------------------------------
# Parser for OCO L2 Input Files
# -----------------------------------------------------------------------------

my_list_type = None
try:
    # This is for python 2
    from types import ListType
    my_list_type = ListType
except ImportError:
    # This is for python 3
    my_list_type = list

import os
import re
import sys
import copy

import warnings
import platform

import six

from .ply_parser import PlyParser

import string
from xml.parsers import expat

class _Node(object):
    def extract_whitespace(self, value):
        frontSearch = re.findall('^[\s\t\n]+', value)
        if len(frontSearch) > 0:
            frontspace = frontSearch[0]
        else:
            frontspace = ''

        endSearch = re.findall('[\s\t\n]+$', value)
        if len(endSearch) > 0:
            endspace = endSearch[0]
        else:
            endspace = ''

        value = re.sub('^[\s\t\n]+', '', value)
        value = re.sub('[\s\t\n]+$', '', value)

        return [value, frontspace, endspace]

    def __init__(self, nodeType, leaf=None, children=None, value_list=None):
         self.value_list = value_list
         self.type = nodeType
         if children:
              self.children = children
         else:
              self.children = [ ]

         if leaf != None:
             if type(leaf) is my_list_type:
                 self.leaf = []
                 self.frontspace = []
                 self.endspace = []
                 
                 for leafVal in leaf:
                     cleanedVals = self.extract_whitespace(leafVal)
                     self.leaf.append(cleanedVals[0])
                     self.frontspace.append(cleanedVals[1])
                     self.endspace.append(cleanedVals[2])
             else:
                 cleanedVals = self.extract_whitespace(leaf)
                 self.leaf = cleanedVals[0]
                 self.frontspace = cleanedVals[1]
                 self.endspace = cleanedVals[2]
         else:
             self.leaf = leaf


class _Section(_Node):
    def __init__(self, leaf=None, children=None):
        if not type(leaf) is my_list_type:
            leaf = [leaf, leaf]

        _Node.__init__(self, 'section', leaf, children)

    def get_all_keyword_names(self):
        return [ child.leaf for child in self.get_all_keyword_objs() ]

    def get_all_keyword_objs(self):
        keywords = []
        for child in self.children:
            if child.type == 'assignment':
                keywords.append( child )

        return keywords

    def get_keyword_value(self, keyword, addWhiteSpace=False, descend=False):
        'Get a particular named keyword'

        keywords = []
        for child in self.children:
            if child.type == 'assignment' and child.leaf.lower() == keyword.lower():
                if addWhiteSpace:
                    valueList = [x.frontspace + x.leaf + x.endspace.rstrip('\n') for x in child.children]
                else:
                    valueList = [x.leaf for x in child.children]
                if len(valueList) == 1:
                    keywords.append( valueList[0] )
                else:
                    keywords.append( valueList )

            if child.type == 'section' and descend:
                child_key_list = child.get_keyword_value(keyword, addWhiteSpace, descend)

                if child_key_list != None:
                    if not type(child_key_list) is my_list_type:
                        child_key_list = [child_key_list]

                    for child_key in child_key_list:
                        keywords.append(child_key)
        if len(keywords) == 0:
            return None
        elif len(keywords) == 1:
            return keywords[0]
        else:
            return keywords

    def get_keywords_dict(self, **kwargs):
        keywordsDict = {}
        for keywordName in self.get_all_keyword_names():
            keywordsDict[keywordName] = self.get_keyword_value(keywordName, **kwargs)
        return keywordsDict

    def set_keyword_value(self, keyword, values, which=None, indent=0):
        'Set a keyword value in the section or add a new keyword assignment'

        if values == None:
            values = [ "" ]
        elif type(values) != my_list_type:
            values = [ values ]

        assignmentNodes = []
        found_count = 0
        for child in copy.copy(self.children):
            # Loop and find all occurances of the assigment.            # if they exist            
            if child.type == 'assignment' and child.leaf.lower() == keyword.lower():
                if which == None or found_count in which:
                    assignmentNodes.append(child)
                    
                found_count = found_count + 1

        # Create a new node if one does not exist and a specific index (which) was not supplied
        if len(assignmentNodes) == 0 and which == None:
            newNode = _Node('assignment', leaf=keyword)
            assignmentNodes.append( newNode )
            self.children.append( newNode )

        # Assign values to the keyword using the existing
        # children nodes if possible
        for currAssignNode in assignmentNodes:
            # Indent front space
            currAssignNode.frontspace += ' ' * int(indent)
            
            prevChildren = currAssignNode.children

            maxNumVals = max(len(values), len(prevChildren))
            currAssignNode.children = []

            lastEndSpace = ''
            for childIdx in range(maxNumVals):

                if childIdx < len(values):
                    if childIdx < len(prevChildren):
                        prevChildren[childIdx].leaf = values[childIdx]
                        space_match = re.search('\n+', prevChildren[childIdx].endspace)
                        if space_match:
                            beg = space_match.start()
                            end = space_match.end()
                            child_space = prevChildren[childIdx].endspace[beg:end]
                            lastEndSpace += child_space
                            prevChildren[childIdx].endspace = prevChildren[childIdx].endspace[:beg] + prevChildren[childIdx].endspace[end:]
                        
                        currAssignNode.children.append(prevChildren[childIdx])
                    else:
                        newChild = _Node('value', leaf=values[childIdx])
                        currAssignNode.children.append(newChild)

            # Ensure that the last node has a newline as whitespace after it
            if len(currAssignNode.children) > 0:
                lastValue = currAssignNode.children[len(currAssignNode.children)-1]
                lastValue.endspace += lastEndSpace

    def delete_keyword(self, keyword, which=None):
        'Delete keyword from a section matching a particular keyword name'

        key_count = 0
        for child in self.children:
            if child.type == 'assignment' and child.leaf == keyword:
                if which == None or key_count in which:
                    self.children.remove(child)
                key_count += 1

        return None

    def get_all_sections(self, node=None, nameTree='', namesList=None, nodesList=None):
        """Returns hierarchical names and flattended nodes in order
        of appearance"""

        if namesList == None:
            namesList = []

        if nodesList == None:
            nodesList = []
        
        if node == None:
            return self.get_all_sections(self)
        else:           
            if node.leaf[0] != 'root':
                if nameTree == '':
                    nameTree = node.leaf[0]
                else:
                    nameTree += '->' + node.leaf[0]

                namesList.append( nameTree )
                nodesList.append( node )
            
            for child in node.children:
                if child.type == 'section':
                    self.get_all_sections(child, nameTree, namesList, nodesList)

            return [namesList, nodesList]

    def get_all_section_names(self):
        'Return list of all section names'
        return self.get_all_sections()[0]

    def get_all_section_nodes(self):
        'Return list of all section nodes'
        return self.get_all_sections()[1]

    def __getitem__(self, key):
        '''This provides a simpler interface to get_section and
        get_keyword_value. We take something like
        "input/PGENameGroup/PGEName" and return
        get_section("input->PGENameGroup")[0].get_keyword_value("PGEName").
        This doesn't do anything you couldn't do by calling these
        directly, but has the advantage of making the most common
        function on L2Input easier.

        This can take either "/" or "->" as a separator.

        To handle matrix results, we first try to retrieve "key->LIST->VALUES".
        If this returns a list, then we return the matrix data (so you can
        read sounding ids by "input/GOSATFullPhysics").'''
        key = key.replace("/", "->")
        # First, try to retrieve as a list
        t = self.get_section(key + "->LIST->VALUES")
        if(len(t) > 0):
            return t[0].get_matrix_data()

        # If this fails, retrieve as a keyword value
        t = key.split('->')
        return self.get_section('->'.join(t[0:-1]))[0].get_keyword_value(t[-1])

    def sounding_ids(self):
        '''This gets the PGE name, and then uses that to read the sounding
        ids. This is used common enough that it is worth having a special
        function for this.'''
        return self["input/" + self["input/PGENameGroup/PGEName"]]

    def get_section(self, sectionName):
        'Return the one or many sections using a hierarchical name'

        if sectionName.find('->') == 0:
            sectionName = self.leaf[0] + sectionName
            
        (sectNames, sectNodes) = self.get_all_sections()

        foundSects = []
        for (name, node) in zip(sectNames, sectNodes):
            if name.upper() == sectionName.upper():
                foundSects.append(node)

        return foundSects

    def delete_section(self, section, node=None, nameTree='', which=None):
        'Delete sections with a given section name'

        if which != None:
            raise Exception('which keyword not yet implemented')

        if node == None:
            self.Delete_Section(section, self)
        else:           
            if node.leaf[0] != 'root':
                if nameTree == '':
                    nameTree = node.leaf[0]
                else:
                    nameTree += '->' + node.leaf[0]
            
            for child in copy.copy(node.children):
                if child.type == 'section':
                    childName = nameTree + '->' + child.leaf[0]
                    if type(section) == 'str' and re.search(section.upper(), childName.upper()):
                        node.children.remove(child)
                    elif section == child:
                        node.children.remove(child)
                    else:
                        self.Delete_Section(section, child, nameTree)

    def get_matrix_data(self):
        'Returns matrix data if contained in the section'

        for child in self.children:
            if child.type == 'matrix':
                matrixValues = []
                # Faster handling if we pulled the values out (e.g.
                # for _AsciiParser)
                if(child.value_list is not None):
                    mval = child.value_list[int(child.children[0][0].leaf)]
                    # Each line is a row, and each entry on a row is a column
                    res = [i.split() for i in mval.split("\n")]
                    # Now flatten this, so a 1 column matrix gets returned as a
                    # vector
                    res = [i[0] if len(i) == 1 else i for i in res
                           if len(i) != 0]
                    return res
                # Drop down to slower handling (e.g., for _XmlParser)
                for rowNode in child.children:
                    rowValues = [ x.leaf for x in rowNode ]

                    if len(rowValues) == 1:
                        matrixValues.append(rowValues[0])
                    else:
                        matrixValues.append(rowValues)

                return matrixValues
        return None

    def set_matrix_data(self, newData):
        'Sets the matrix data into the section'

        if not type(newData) is my_list_type:
            newData = [newData]

        # Create children nodes for new data
        newNodeChildren = []
        for rowData in newData:
            if not type(rowData) is my_list_type:
                rowData = [rowData]
                
            newRow = []
            for colValue in rowData:
                if(sys.version_info > (3,) and
                   isinstance(colValue, bytes)):
                    newRow.append( _Node('value', leaf=colValue.decode("utf-8") + '\t' ) )
                else:
                    newRow.append( _Node('value', leaf=str(colValue) + '\t' ) )
                
            lastValue = newRow[len(newRow)-1]
            lastValue.endspace = '\n'

            newNodeChildren.append(newRow)

        foundExistingMatrix = False
        for child in self.children:
            if child.type == 'matrix':
                # Add new matrix children to first found matrix node
                child.children = newNodeChildren
                foundExistingMatrix = True
                break

        # Append matrix data to section since a matrix object does not already exist
        if not foundExistingMatrix:
            self.children.append( _Node('matrix', children=newNodeChildren) )

    def add_child(self, newChild):
        if not type(newChild) is my_list_type:
            newChild = [newChild]

        for childObj in newChild:
            self.children.append(childObj)

    # FAO
    def del_child(self, deadChild):
        if not type(deadChild) is my_list_type:
            deadChild = [deadChild]

        for childObj in deadChild:
            self.children.remove(childObj)

#####################
# Parse ASCII L2 Input file

class _AsciiParser(PlyParser):

    def __init__(self, value_list=None, contents='', **kw):
        'Initialize the class'
        self.value_list=value_list
        PlyParser.__init__(self, **kw)

    ##############
    # Lexing rules

    lastToken = None
    frontSpace = ''

    keywords = (
        'begin', 'end',
        )

    tokens = keywords + (
        'KEYWORD', 'QUOTE_KEYWORD', 'COMMENT', 'DOUBLE_QUOTE', 'SINGLE_QUOTE', 'DATETIME',
        'RANGE', 'INT_TEXT', 'INTEGER', 'FLOAT', 'TEXT'
        )

    # Tokens

    def t_RANGE(self, t):
        r'\d+\:\d+'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_DOUBLE_QUOTE(self, t):
        r'"""(?:[^"]|\\"|"{1,2}(?!"))*"""|"(?:[^"]|\\")*"(?!")'
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_QUOTE_KEYWORD(self, t):
        r"'(?:[^'])*'(?!')\s+\=" # Nested quotes not allowed
        t.value = self.frontSpace + t.value.strip("'").replace("'", "")
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_SINGLE_QUOTE(self, t):
        r"'''(?:[^']|\\'|'{1,2}(?!'))*'''|'(?:[^']|\\')*'(?!')"
        t.value = self.frontSpace + t.value.strip("'")
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_DATETIME(self, t):
        r'\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d.\d\d\dZ'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_FLOAT(self, t):
        r'((\d*\.\d+)([EeDd][\+-]?\d+)?|([1-9]\d*[EeDd][\+-]?\d+)|(\d+\.))'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_INT_TEXT(self, t):
        r'\d+[a-zA-Z][a-zA-Z0-9]*'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_INTEGER(self, t):
        r'\d+'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_COMMENT(self, t):
        r'\#.*'
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t

    def t_KEYWORD(self, t):
        r'[a-zA-Z][a-zA-Z0-9\(\)_:]*\s+\='
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t
   
    def t_TEXT(self, t):
        r'[a-zA-Z_.\/~<>\-\[\]\*\+\(\)\?\%\$][a-zA-Z0-9_/.~<>\-:,\[\]\*\+\(\)\?\%\{\}\$]*'
        if t.value in self.keywords:
            t.type = t.value.lower()
        t.value = self.frontSpace + t.value
        self.frontSpace = ''
        self.lastToken = t
        return t
   
    def t_end_whitespace(self, t):
        r'[\s\t\n]+'
        if self.lastToken != None:
            if t.value.find('\n') >= 0:
                self.lastToken.value += t.value[0:t.value.rfind('\n')+1]
                self.frontSpace = t.value[t.value.rfind('\n')+1:]
            else:
                self.lastToken.value += t.value
        else:
            self.frontSpace = t.value
        t.lexer.lineno += t.value.count("\n")
    
    def t_error(self, t):
        next_space_loc = t.value.find(' ')
        raise ValueError("Illegal value '%s'... at line %d. Last token: %s" % (t.value[0:next_space_loc+1], t.lexer.lineno, self.lastToken))
        t.lexer.skip(1)

    ###############
    # Parsing rules

    def p_file_root(self, p):
        'file_root : section_contents'
        p[0] = _Section(leaf=['root','root'], children=p[1])
        p[0].filename = self.filename

    def p_section_contents_comments(self, p):
        'section_contents : COMMENT section_contents'
        newNode = _Node('comment', leaf=p[1])
        if type(p[2]) is my_list_type:
            p[0] = [ newNode ] + p[2]
        else:
            p[0] = [ newNode, p[2] ]

    def p_section_contents_assignment(self, p):
        """
        section_contents : assignment section_contents
                         | section_block section_contents
        """
        if type(p[2]) is my_list_type:
            p[0] = [ p[1] ] + p[2]
        else:
            p[0] = [ p[1], p[2] ]

    def p_section_contents_value_list(self, p):
        'section_contents : value_list'
        matrixData = [[]]
        rowIndex = 0
        for value in p[1]:
            matrixData[rowIndex].append(value)

            if value.endspace.count('\n') > 0:
                rowIndex += 1
                matrixData.append([])

        # Delete last row which will be empty if newline after last
        # matrix row
        if len(matrixData[rowIndex]) == 0:
            del(matrixData[rowIndex])
        p[0] = [ _Node('matrix', children=matrixData,
                       value_list=self.value_list) ]

    def p_section_contents_empty(self, p):
        'section_contents : empty'
        p[0] = _Node('empty')

    def p_section_block(self, p):
        'section_block : begin TEXT section_contents end TEXT'

        if not type(p[3]) is my_list_type:
            p[0] = _Section(leaf=[p[2], p[5]], children=[p[3]])
        else:
            p[0] = _Section(leaf=[p[2], p[5]], children=p[3])

        begin_ws = p[0].extract_whitespace(p[1])
        end_ws   = p[0].extract_whitespace(p[4])

        p[0].frontspace = [ begin_ws[1], end_ws[1] ]
        p[0].filename = self.filename
        #p[0].endspace   = [ begin_ws[2], end_ws[2] ]

    def p_assignment(self, p):
        """
        assignment : KEYWORD
                   | QUOTE_KEYWORD
                   | QUOTE_KEYWORD value_list
                   | KEYWORD value_list
        """

        keywordName = (p[1])[:p[1].rfind('=')]
        extraws = (p[1])[p[1].rfind('=')+1:]

        if len(p) == 3:           
            p[0] = _Node('assignment', leaf=keywordName, children=p[2])
        elif len(p) == 2:
            blankVal = _Node('value', leaf='')
            blankVal.endspace = extraws
            p[0] = _Node('assignment', leaf=keywordName, children=[blankVal])

    def p_value_list(self, p):
        """
        value_list : value_list value
                   | value
        """
        if len(p) == 3:
            p[0] = p[1] + [ p[2] ]
        elif len(p) == 2:
            p[0] = [ p[1] ]

    def p_value(self, p):
        """
        value : TEXT
              | DATETIME
              | RANGE
              | INT_TEXT
              | INTEGER
              | FLOAT
              | DOUBLE_QUOTE
              | SINGLE_QUOTE
        """
        p[0] = _Node('value', leaf=p[1])

    def p_empty(self, p):
        'empty :'
        pass

    def p_error(self, p):
        if p == None:
            raise ValueError("Syntax error with no info!")
        else:
            print(dir(p))
            raise ValueError("Syntax error at token '%s' with value '%s' at line number %d of file %s" % (p.type, p.value, p.lineno, self.filename))

########

class _XmlParser(object):
    "Parse XML L2 Input into Sections and Nodes"

    def __init__(self, contents='', **kw):
        # For XML Parsing
        self.nodeStack  = []

    def start_element(self,name,attributes):
        "Sets the node stack based on elements encountered in the XML file"

        if str(name).lower() == 'group' and 'name' in list(attributes.keys()):
            # Turn group elements into Sections
            group_name = attributes.pop('name')
            element = _Section(leaf=[str(group_name), str(group_name)])
        elif str(name).lower() == 'scalar' and 'name' in list(attributes.keys()):
            # Turn scalars into Assigments
            key_name = attributes.pop('name')
            element = _Node('assignment', leaf=key_name)
        elif str(name).lower() == 'vector' and 'name' in list(attributes.keys()):
            # For vector elements:
            # Create a Section(LIST) -> Section(VALUES)
            # like would be done in ASCII config
            # place a vector on the node stack where vector elements will be
            # appended 
            vec_name = attributes.pop('name')
            list_sect = _Section(leaf=['LIST', 'LIST'])

            if len(self.nodeStack) > 0:
                parent = self.nodeStack[-1]
                parent.add_child(list_sect)
            else:
                self.rootNode = list_sect

            self.nodeStack.append(list_sect)
            
            list_sect.set_keyword_value('name', vec_name)
            vals_sect = _Section(leaf=['VALUES', 'VALUES'])
            list_sect.add_child(vals_sect) 

            self.nodeStack.append(vals_sect)

            element = [] 
        elif str(name).lower() == 'element':
            # No new items for an element, just place data into values section
            return
        else:
            # Default to catch all other types of elements 
            element = _Section(leaf=[str(name), str(name)])

        # If the current terminal element is a section and there are attributes
        # to the XML node, then add these attributes to the section
        if hasattr(element, 'type') and element.type == 'section':
            for attr_key in list(attributes.keys()):
                element.set_keyword_value(attr_key, attributes[attr_key])
        elif len(list(attributes.keys())) > 0:
                raise ValueError('Element type %s can not support additional attributes' % element.type)
        
        # If the current terminal element is not a list (for vector parsing) then 
        # Add it as a child to the last parent object, and only set the root node if the
        # node stack is empty
        if hasattr(element, 'type') and len(self.nodeStack) > 0:
            parent = self.nodeStack[-1]
            parent.add_child(element)
        elif len(self.nodeStack) == 0:
            self.rootNode = element
            
        self.nodeStack.append(element)
        
    def end_element(self,name):
        if str(name).lower() == 'vector':
            # vector has ended so add the list as the matrix data of the VALUES Section
            # under the LIST Section. Then make the LIST Section the bottom of the node stack
            self.nodeStack[-2].set_matrix_data(self.nodeStack[-1])
            self.nodeStack = self.nodeStack[:-3]
        elif str(name).lower() != 'element':
            # For any other ending XML node other than element (ie vector data)
            # pop the object off the node stack
            self.nodeStack = self.nodeStack[:-1]

    def character_data(self,data):

        if data.strip():
            data = data.encode()

            element = self.nodeStack[-1]
            if type(element) is my_list_type: 
                element.append(data)
            elif element.type == 'assignment':
                if(sys.version_info > (3,) and
                   isinstance(data, bytes)):
                    element.children.append( _Node('value', leaf=data.decode("utf-8")))
                else:
                    element.children.append( _Node('value', leaf=str(data)) )
            else:
                raise ValueError('Element of type %s does not support setting character data' % element.type)
                
            return

    def parse(self,fileContents):
        # Create a SAX parser
        expatParser = expat.ParserCreate()

        # SAX event handlers
        expatParser.StartElementHandler = self.start_element
        expatParser.EndElementHandler = self.end_element
        expatParser.CharacterDataHandler = self.character_data

        # Ensure that expat does not split values between multiple CharacterDataHandler calls
        expatParser.buffer_text = True

        # Parse the XML File
        ParserStatus = expatParser.Parse(fileContents, 1)
       
        return self.rootNode


######################
# Main Interface Class

class L2InputFile(object):

    def __init__(self, file_input=None, **kw):
        'Initialize the class'

        # Parsed data
        self.rootNode = None

        # Initialize root node to a default value
        self.rootNode = _Section(leaf=['root','root'])

        self.filename = None
        if file_input != None:
            self.read(file_input)

    def __getitem__(self, key):
        return self.rootNode[key]

    def __getattr__(self, name):
        if name in dir(self.rootNode):
            return eval("self.rootNode.%s" % name)
        else:
            raise AttributeError('Attribute %s not found in rootNode or in class %s' % (name, __name__))

    def read(self, file_input):
        'Read data from disk optionally specifying an alternative filename'

        if isinstance(file_input, six.string_types):
            self.filename = file_input

            if not os.path.exists(self.filename):
                raise IOError('Can not read non-existant file: %s' % self.filename)
            
            inputFile = open(self.filename, 'r')
            fileContents = inputFile.read()
            inputFile.close()
        elif hasattr(file_input, 'read'):
            fileContents = file_input.read()
        else:
            raise IOError('Passed unknown object for reading: %s' % file_input)

        if len(fileContents.strip()) == 0:
            self.rootNode = _Section(leaf=['root','root'], children=[])
        else:
            is_xml = False
            for fileLine in fileContents.split():
                fileLine = fileLine.strip()
                if len(fileLine) > 0:
                    if fileLine.lower().find('<?xml') >= 0:
                        is_xml = True
                    break

            # The _AsciiParser can be extremely slow for large files.
            # We pull out any values to handle separately. We have a array
            # value_list that contains the contents of the value, we replace
            # this in fileContentsShort with the index number for this. So
            # _AsciiParser never actually sees the potentially long list of
            # values, which allows it to run much faster. We translate the
            # value back to the contents when we eventually look up the value.
            value_list = None
            if(not is_xml):
                fileContentsShort = ""
                value_list = []
                in_value = False
                for t in re.split('(begin|end)\s+values', fileContents, flags=re.IGNORECASE):
                    if(t.lower() == "begin"):
                        if(in_value):
                            raise RuntimeError("Confused processing %s" % self.filename)
                        in_value = True
                        fileContentsShort += "begin VALUES\n"
                        fileContentsShort += "%d\n" % len(value_list)
                        value_list.append("")
                    elif(t.lower() == "end"):
                        if(not in_value):
                            raise RuntimeError("Confused processing %s" % self.filename)
                        in_value = False
                        fileContentsShort += "end VALUES\n"
                    else:
                        if(in_value):
                            value_list[-1] = t
                        else:
                            fileContentsShort += t

            if (is_xml):
                self.rootNode = _XmlParser().parse(fileContents)
            else:
                self.rootNode = _AsciiParser(value_list=value_list,filename=self.filename).parse(fileContentsShort)


    def write(self, file_output=None, doIndent=False):
        'Write data to disk with outputFilename'

        if file_output == None:
            file_output = self.filename

        if isinstance(file_output, six.string_types):
            file_obj = open(file_output, "w")
            self.filename = file_output
        elif hasattr(file_output, 'write'):
            file_obj = file_output
        else:
            raise IOError('Passed unknown object type for writing: %s' % file_output)
        
        self.write_node(self.rootNode, file_obj, doIndent=doIndent)

        if isinstance(file_output, six.string_types):
            file_obj.close()
                
    def write_node(self, node, fileObj, level=0, doIndent=False, nextNode=None, lastNode=None):
        'Write contents of node recursively to file object'

        indent = '   ' * level

        if node == None:
            raise Exception('Can not write non-existant node!')

        if node.type == 'comment':

            lastNewline = True
            if lastNode != None and len(lastNode.children) > 0:
                lastChild = lastNode.children[len(lastNode.children)-1]
                if lastChild.type != 'empty' and lastChild.leaf != None:
                    lastNewline = lastChild.endspace.count('\n') > 0
                else:
                    lastNewline = True

            if doIndent and lastNewline:
                comment_fs = indent
            else:
                comment_fs = node.frontspace
                
            comment_es = node.endspace
            
            if comment_es == None or len(comment_es) == 0:
                comment_es = '\n'
                
            fileObj.write( comment_fs + node.leaf + comment_es )
                
        elif node.type == 'assignment':
            valuesList = []

            for child in node.children:
                if child.frontspace != None and len(child.frontspace) > 0:
                    valuesList.append(child.frontspace)
                    
                valuesList.append(child.leaf)

                if child.endspace == None or len(child.endspace) == 0:
                    valuesList.append(' ')
                else:
                    valuesList.append(child.endspace)

            # Preserve any empty values list exactly
            valuesStr = ''
            for oneValue in valuesList:
                valuesStr += oneValue

            if len(valuesStr.rstrip()) > 0:
                valuesStr = ' ' + valuesStr

            if valuesStr.count('\n') < 1 and (nextNode == None or (nextNode != None and nextNode.type != 'comment')):
                valuesStr += '\n'

            if doIndent:
                node_fs = indent
            else:
                node_fs = node.frontspace
                
            node_es = node.endspace
            if node_es == None or len(node_es) == 0:
                node_es = ' '

            if not re.search('^[a-zA-Z]', node.leaf):
                node_fs += "'"
                node_es  = "'" + node_es
            
            fileObj.write( node_fs + node.leaf + node_es + '=' + valuesStr )

        elif node.type == 'section':
            if node.leaf[0] == 'root':
                for child in node.children:
                    self.write_node(child, fileObj, level, doIndent)
            else:
                if doIndent:
                    node_fs1 = indent
                else:
                    node_fs1 = node.frontspace[0]

                node_es1 = node.endspace[0]
                if node_es1 == None or len(node_es1) == 0 or node_es1 == '':
                    node_es1 = '\n'
                
                fileObj.write( '%sbegin %s%s' % (node_fs1, node.leaf[0], node_es1) )

                childLast = None
                for child_idx in range(len(node.children)):
                    child = node.children[child_idx]
                    
                    if child_idx+1 < len(node.children):
                        childNext = node.children[child_idx+1]
                    else:
                        childNext = None
                    
                    self.write_node(child, fileObj, level+1, doIndent, nextNode=childNext, lastNode=childLast)

                    childLast = child

                if doIndent:
                    node_fs2 = indent
                else:
                    node_fs2 = node.frontspace[1]

                node_es2 = node.endspace[1]
                if node_es2 == None or len(node_es2) == 0 or node_es2 == '':
                    node_es2 = '\n\n'
                elif doIndent and node_es2.count('\n') == 1:
                    node_es2 += '\n'

                fileObj.write( '%send %s%s' % (node_fs2, node.leaf[1], node_es2) )

        elif node.type == 'matrix':
            for row in node.children:
                if doIndent:
                    fileObj.write(indent)
                    
                for colVal in row:
                    fileObj.write(colVal.frontspace + colVal.leaf + colVal.endspace)

        elif node.type == 'value':
            fileObj.write(node.leaf + node.endspace)

        elif node.type == 'empty':
            pass
                
        else:
            errMsg = 'Unknown node type: %s' % node.type
            raise Exception(errMsg)

