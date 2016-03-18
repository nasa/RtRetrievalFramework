#!/usr/bin/env python

from builtins import range

import os
import re
import sys
import inspect
from optparse import OptionParser

from full_physics import l2_analysis

def standalone_main():
    # Set up command line arguments
    parser = OptionParser(usage="usage: %prog [options] [l2_file_1] [l2_file_2]...")

    parser.add_option( "-n", "--name", dest="obj_names",
                       metavar="NAME",
                       type="string",
                       action="append",
                       default=None,
                       help="alternative names to be each filename passed in, specify once for each file"
                       )

    parser.add_option( "-s", "--script", dest="script",
                       metavar="FILENAME",
                       type="string",
                       default=None,
                       help="run the passed script file in the same context as would be present in the interactive shell"
                       )

    parser.add_option( "-a", "--ancillary_file", dest="ancillary_files",
                       metavar="NAME",
                       type="string",
                       action="append",
                       default=None,
                       help="additional ancillary files for querying data which match the size of the L2 files"
                       )

    parser.add_option( "-r", "--reference_file", dest="reference_file",
                       metavar="NAME",
                       type="string",
                       default=None,
                       help="a single ancillary files that is repeated for all input files"
                       )


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("Must specify at least one file for analysis!")

    if(options.ancillary_files != None and options.reference_file != None):
        parser.error("The options --ancillary_file or --reference_file can not be used together")
    elif options.ancillary_files != None:
        addl_files = options.ancillary_files
    elif options.reference_file != None:
        addl_files = [ options.reference_file for idx in range(len(args)) ]
    else:
        addl_files = None
           
    analysis_env = l2_analysis.AnalysisEnvironment(args, obj_names=options.obj_names, addl_files=addl_files)

    routine_objs = []
    for item_name, mod_item in inspect.getmembers(l2_analysis):
        if inspect.isclass(mod_item) and re.search("Routines$", item_name):
            routine_objs.append( mod_item(analysis_env=analysis_env) )
    analysis_env.add_routines(routine_objs)

    if options.script != None:
        analysis_env.run_script(options.script)
    else:
        analysis_env.launch_shell()
    
if __name__ == "__main__":
    standalone_main()


