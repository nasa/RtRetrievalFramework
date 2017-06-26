import os
import sys
import copy
from pyIDL import idl

class IDL_Util:

    def __init__(self, sav_filename=None):
        idl_routines = []

        if sav_filename != None:
            self.sav_filename = sav_filename
        else:
            self.sav_filename = '%s/%s' % (os.path.dirname(sys.argv[0]), str(self.__class__).split('.')[-1]) + '.sav'

    def load_idl(self, recompile=False):
        self.setup_idl()

        if not recompile and os.path.exists(self.sav_filename):
            self.load_sav()
        else:
            self.compile_sav()
 
    def setup_idl(self, stdout=True):
        # Remove IDL_STARTUP from environment before loading IDL
        if 'IDL_STARTUP' in os.environ:
            del(os.environ['IDL_STARTUP'])
        
        # Load the IDL interpreter
        self.idl = idl(stdout=stdout)

        self.idl.eval('comp_routines = ROUTINE_NAMES(/PROCEDURES)')
        comp_routines = self.idl.get('comp_routines')
        self.idl.delete('comp_routines')

        if not 'WORK' in comp_routines:
            # Find where IDL utilities are located
            l2_support_path = os.getenv('L2_SUPPORT_PATH')

            if len(l2_support_path) == 0:
                raise EnvironmentError('L2_SUPPORT_PATH environmental variable must be defined using setup_env.sh')

            idl_utils_path = l2_support_path + '/IDL_Utils'

            # Set up L2 IDL environment
            self.idl.eval('.compile ' + idl_utils_path + '/work.pro')

            self.idl.work()
        
    def compile_sav(self, sav_filename=None):
        '''Compiles a .sav file of IDL code used by the software'''

        if sav_filename != None:
            self.sav_filename = sav_filename
        
        # Ensure IDL objects are fully loaded
        for r_name in self.idl_routines:
            self.idl.eval('.compile ' + r_name)
        
        self.idl.eval('Resolve_All, /CONTINUE')
        self.idl.eval('Save, /ROUTINE, FILE="' + self.sav_filename + '"')

    def load_sav(self, sav_filename=None):
        '''Loads a compiled .sav file back into the IDL session'''

        if sav_filename != None:
            self.sav_filename = sav_filename

        self.idl.restore(self.sav_filename)
