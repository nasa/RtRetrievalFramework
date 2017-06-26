import os
import bisect
import logging

import h5py

import Orbit_Sim

from Generator_Utils import *
from OCO_TextUtils import evaluate_bool_str

SOUNDING_ID_GROUP   = 'SoundingGeometry'
SOUNDING_ID_DATASET = 'sounding_id'

FRAME_ID_GROUP   = 'FrameHeader'
FRAME_ID_DATASET = 'frame_id'

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
   logger = logging.getLogger(os.path.basename(__file__))

   if len(moduleSections) > 1:
       raise RuntimeError('Only one value map creation allowed per FILE')

   if str(source) == str(destination):
      raise IOError('source and destination must be different. will not overwrite source file')

   list_file    = Apply_Template(moduleSections[0].Get_Keyword_Value('list_file'), valuesDict, mapDict=mapDict)
   data_file    = Apply_Template(moduleSections[0].Get_Keyword_Value('data_file'), valuesDict, mapDict=mapDict)
   data_col     = Apply_Template(moduleSections[0].Get_Keyword_Value('data_column'), valuesDict, mapDict=mapDict)
   section      = Apply_Template(moduleSections[0].Get_Keyword_Value('section'), valuesDict, mapDict=mapDict)
   static_value = Apply_Template(moduleSections[0].Get_Keyword_Value('static_value'), valuesDict, mapDict=mapDict)
   is_log_file  = evaluate_bool_str(moduleSections[0].Get_Keyword_Value('is_log_file'))
   l1b_file     = Apply_Template(moduleSections[0].Get_Keyword_Value('l1b_file'), valuesDict, mapDict=mapDict)
   modify       = moduleSections[0].Get_Keyword_Value('modify')

   max_range_val = None
   range_values = {}
   for range_sect in moduleSections[0].Get_Section('->RANGES'):
      logger.debug('Using range section')

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

   id_list = Read_Id_List_File(list_file, section, valuesDict=valuesDict, mapDict=mapDict)

   data_values = []
   if data_file != None:
      
      if len(data_file) == 0 or not os.path.exists(data_file):
         raise IOError('Could not read data_file')

      if is_log_file:
         if l1b_file == None or len(l1b_file) == 0:
            raise ValueError('Need L1B file specified for using log file as source of data')
         if not os.path.exists(l1b_file):
            raise IOError('L1B file specified does not exist: %s' % l1b_file)
         
         log_file_obj = Orbit_Sim.Log_File(data_file)
         col_index = log_file_obj.get_column_index(data_col)

         if not type(col_index) is ListType:
            col_index = [ col_index ]

         h5_obj = h5py.File(l1b_file, 'r')
         snd_id_matrix = h5_obj[SOUNDING_ID_GROUP][SOUNDING_ID_DATASET]
         frame_id_arr  = h5_obj[FRAME_ID_GROUP][FRAME_ID_DATASET]

         for curr_sounding in id_list:
            curr_frame_id = int(str(curr_sounding)[0:-1])
            frame_index = bisect.bisect_left(frame_id_arr, curr_frame_id)

            for snd_index in range(snd_id_matrix.shape[1]):
               if snd_id_matrix[frame_index, snd_index] == int(curr_sounding):
                  break

            if snd_id_matrix[frame_index, snd_index] != int(curr_sounding):
                raise ValueError('did not find correct sounding id: %d at index: %s in hdf file: %s, instead found: %d' % (curr_sounding, (frame_index, snd_index), l1b_file, snd_id_matrix[frame_index, snd_index]))
             
            curr_log_val = 0.0
            for curr_val_idx in col_index:
               curr_log_val += log_file_obj.data[frame_index, snd_index, curr_val_idx]

            data_values.append(curr_log_val)

      else:
          if data_col == None:
             data_col = 0
          else:
             data_col = int(data_col)

          logger.debug('Reading mapped values from column %d of file %s' % (data_col, data_file))
          data_fobj = open(data_file)
          for data_line in data_fobj.readlines():
               data_line = data_line.strip()
               if len(data_line) > 0 and data_line.find('#') != 0:
                   line_parts = data_line.split()

                   if len(line_parts)-1 < data_col:
                       raise IOError('data file %s does not have column %d' % (data_file, data_col))
                   data_values.append(line_parts[data_col])

   if static_value != None:
      logger.debug('Setting mapped value to static value: %s' % static_value)
      for idx in range(len(id_list) - len(data_values)):
         data_values.append(static_value)

   if len(id_list) != len(data_values):
       raise IOError('Length of id list %d from file %s does not match length of data values %d from %s' % (len(id_list), list_file,  len(data_values), data_file))

   mapValues = None
   mapSects = moduleSections[0].Get_Section('->MAP')
   if mapSects != None and len(mapSects) > 0:
      mapValues = Get_Map_Values(mapSects, valuesDict)

   logger.debug('Writing map file: %s' % destination)

   if type(destination) is str:
      dstFileObj = open(destination, 'w')
   elif hasattr(destination, 'write'):
      dstFileObj = destination
   else:
      raise Exception('Unrecognized source object: %s' % destination)


   if modify != None and len(modify) > 0:
      modifyDict = copy.copy(valuesDict)
      
   for (id_val, data_val) in zip(id_list, data_values):
      if modify != None and len(modify) > 0:
         modifyDict['original'] = str(data_val)
         modify_expr = Apply_Template(modify, modifyDict, mapDict=mapDict)
         data_val = eval(modify_expr)

      if len(range_values) > 0:
         found_range_value = False
         for (curr_name, curr_values) in range_values.items():
            beg_val = curr_values[0]
            end_val = curr_values[1]

            if float(data_val) >= beg_val and float(data_val) < end_val:
               data_val = curr_name
               found_range_value = True
               break

         if not found_range_value:
            raise LookupError('RANGE values specified but none matched for value: %s' % data_val)
      
      if mapValues != None and (str(data_val) in mapValues[DEFAULT_MAP_NAME]):
         print >>dstFileObj, id_val, mapValues[DEFAULT_MAP_NAME][str(data_val)]
      else:
         print >>dstFileObj, id_val, str(data_val)

   if type(destination) is str:
      dstFileObj.close()

