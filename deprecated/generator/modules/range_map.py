import os
import logging

from Generator_Utils import *
from OCO_TextUtils import index_range_list, evaluate_bool_str
from OCO_Matrix import OCO_Matrix

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
   logger = logging.getLogger(os.path.basename(__file__))

   for currSection in moduleSections:
       if str(source) == str(destination):
           raise IOError('source and destination must be different. will not overwrite source file')

   rows          = Apply_Template(moduleSections[0].Get_Keyword_Value('rows'), valuesDict, mapDict=mapDict)
   columns       = Apply_Template(moduleSections[0].Get_Keyword_Value('columns'), valuesDict, mapDict=mapDict)
   identifier    = Apply_Template(moduleSections[0].Get_Keyword_Value('identifier'), valuesDict, mapDict=mapDict)
   initial_value = Apply_Template(moduleSections[0].Get_Keyword_Value('initial_value'), valuesDict, mapDict=mapDict)
   map_filename  = Apply_Template(moduleSections[0].Get_Keyword_Value('map_filename'), valuesDict, mapDict=mapDict)
   modify        = moduleSections[0].Get_Keyword_Value('modify')

   # Load ranges from RANGES section of module
   max_range_val = None
   range_values = {}
   for range_sect in moduleSections[0].Get_Section('->RANGES'):
      for range_spec in range_sect.Get_Matrix_Data():
         (range_name, range_str) = range_spec

         if range_str.find(',') > 0:
            curr_range = [ float(val) for val in range_str.split(',') ]
         else:
            curr_range = [ float(val) for val in range_str.split() ]

         if max_range_val == None:
            max_range_val = max(curr_range)
         else:
            max_range_val = max(max_range_val, max(curr_range))

         range_values[range_name] = curr_range

   if len(range_values) == 0:
      logger.error('No index range list supplied for operating on source: %s' % source)
      return

   # Load source for data to map agains
   data_obj = OCO_Matrix(source)

   # Set columns to all if argument not supplied,
   # Otherwise try parsing as an index range list failing that try
   # using the specified columns as label names
   if columns == None:
      columns = range(data_obj.dims[1])
   else:
      try:
         columns = index_range_list(columns)
      except ValueError, TypeError:
         columns = data_obj.find_labels(columns, match_case=False, indexes=True)

   # Set rows to all if no argument supplied, otherwise use the specified
   # index range list
   if rows == None:
      rows = range(data_obj.dims[0])
   else:
      rows = index_range_list(rows)

   # Load the modifer string for specifying how to operate on rows and col values
   # By default sum all values encountered
   if modify == None:
      modify = '<current> + <data_result>'

   # Set initial value for before we look at source data file
   if initial_value == None:
      data_result = 0
   else:
      # Use eval so string is evaluated to the type expected
      data_result = eval(intial_value)

   # Load value to check against ranges
   modify_dict = copy.copy(valuesDict)

   # Process data source with modifying expression
   for row_idx in rows:
      for col_idx in columns:
         modify_dict['data_result'] = str(data_result)
         modify_dict['current']     = str(data_obj.data[row_idx, col_idx])
         modify_expr = Apply_Template(modify, modify_dict, mapDict=mapDict)
         data_result = eval(modify_expr)

   # Find match of data_result in the RANGE section specified
   mapped_value = None
   for (curr_name, curr_values) in range_values.items():
      beg_val = curr_values[0]
      end_val = curr_values[1]

      if beg_val <= float(data_result) < end_val:
         mapped_value = curr_name
         break
      
   if mapped_value == None:
      err_msg = 'RANGE values specified but none matched for value: %s using source: %s' % (data_result, source)
      logger.error(err_msg)
      raise Exception(err_msg)

   # Finally append to a map_filename if specified or set into a destination keyword if map_filename != None:
      if type(map_filename) is str:
         dest_obj = open(map_filename, 'a')
      elif hasattr(map_filename, 'write'):
         dest_obj = map_filename
      else:
         raise TypeError('Unknown object type used for map_filename: %s' % map_filename)

      print >>dest_obj, identifier, mapped_value
         
   else:
      # Return mapped value into values dictionary for use elsewhere
      valuesDict[identifier] = mapped_value

