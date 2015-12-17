import os
import sys
import glob
import logging

from Generator_Utils import *

from OCO_TextUtils import evaluate_bool_str

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

   if len(moduleSections) > 1:
       raise RuntimeError('Only one setup_absco_files block allowed')

   absco_path  = Apply_Template(moduleSections[0].Get_Keyword_Value('absco_path'), valuesDict, mapDict=mapDict)
   prefer_file = Apply_Template(moduleSections[0].Get_Keyword_Value('prefer_file'), valuesDict, mapDict=mapDict)

   if absco_path[-1] != '/':
       absco_path += '/'

   inp_file = L2_Input.Input_File(source)

   sound_info_sect = inp_file.Get_Section('SOUNDING_INFO')
   gas_sections    = inp_file.Get_Section('PARAMETER_DEFINITION->GAS')
    
   sound_info_sect[0].Set_Keyword_Value('absco_path', absco_path)

   for curr_gas in gas_sections:
       # Initialize a list of globs to try and go through to find a single filename matching
       file_glob_checklist = []

       # Use old filename for building search globs
       (old_file, old_ext) = os.path.splitext(curr_gas.Get_Keyword_Value('absco_file'))

       # Get gasname and search first of all for anything that matches starting with the gas name
       gas_name_str = curr_gas.Get_Keyword_Value('name').lower()
       file_glob_checklist.append( os.path.join(absco_path, '%s_*%s' % (gas_name_str, old_ext)) )
       
       # Build a list of globs that progressively add more wildcards for parts of existing filename
       old_parts = old_file.split('_')

       for change_idx in range(len(old_parts)-1, -1, -1):
           old_parts[change_idx] = '*'
           file_glob_checklist.append( os.path.join(absco_path, '_'.join(old_parts) + old_ext) )

       # Go through list of search globs and see if we can match the first time we find
       # a unique file given our set of defined preferences from the configuration
       new_file = None
       for srch_count, curr_search_glob in enumerate(file_glob_checklist):
           found_files = glob.glob(curr_search_glob)

           prefer_results = []
           if len(found_files) == 1:
               new_file = found_files[0]
           elif prefer_file != None:
               for curr_find in found_files:
                   for pref_idx, pref_str in enumerate(prefer_file):
                       if curr_find.find(pref_str) >= 0:
                           prefer_results.append( (curr_find, pref_idx) )

               def rank_prefer(a, b):
                   return cmp(a[1], b[1])

               if len(prefer_results) > 0:
                   prefer_results.sort(rank_prefer)
                   new_file = prefer_results[0][0]
               elif srch_count == len(file_glob_checklist)-1:
                  raise Exception('Could not find preferred result among list: %s for glob: %s to match old file: %s%s' % (found_files, curr_search_glob, old_file, old_ext))

           else:
               raise Exception('Could not determine which absco file to use from list: %s for glob: %s to match old file: %s%s' % (found_files, curr_search_glob, old_file, old_ext))

           if (new_file != None and 0 <= len(prefer_results) <= 1):
               break

       if new_file == None:
           raise Exception('Could not find absco filename in dir: %s for gas: %s to match old file: %s%s' % (absco_path, curr_gas.Get_Keyword_Value('name'), old_file, old_ext))

       curr_gas.Set_Keyword_Value('absco_file', os.path.basename(new_file))

   inp_file.Write(destination)
