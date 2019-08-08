from builtins import range
import os
import re
import sys
import copy
import datetime

my_list_type = None
try:
    # This is for python 2
    from types import ListType
    my_list_type = ListType
except ImportError:
    # This is for python 3
    my_list_type = list

ASCII_TIMESTAMP_FORMAT  = '%Y-%m-%dT%H:%M:%S.%fZ' # ie 2006-09-14T11:45:56.001Z

def escape_underscore(string):
   string2 = r''
   for c in string:
      if (c == '_'):
         string2 = string2 + "\\"
      string2 = string2 + c
   return string2

def split_path(string):
   dirs = string.split('/')

   path = ""
   for i in range (0, len(dirs) - 1):
      path = path + dirs[i] + "/"

   file = dirs[len(dirs)-1]

   return (path, file)

def TeX_string(string):
   string2 = r''
   exponent = False
   for c in string:
      if (c == '-'):
         string2 = string2 + "$^{-"
         exponent = True
         continue
      elif (c == '^'):
         string2 = string2 + "$^{"
         exponent = True
         continue
      if (exponent == True and c.isdigit() == False):
         string2 = string2 + "}$"
         exponent = False
      string2 = string2 + c

   if (exponent == True): string2 = string2 + "}$"
   return string2

# relpath
# R.Barran 30/08/2004

def relpath(target, base=os.curdir):
    """
    Return a relative path to the target from either the current dir or an optional base dir.
    Base can be a directory specified either as absolute or relative to current dir.
    """

    if not os.path.exists(target):
        raise OSError('Target does not exist: '+target)

    if not os.path.isdir(base):
        raise OSError('Base is not a directory or does not exist: '+base)

    base_list = (os.path.abspath(base)).split(os.sep)
    target_list = (os.path.abspath(target)).split(os.sep)

    # On the windows platform the target may be on a completely different drive from the base.
    if os.name in ['nt','dos','os2'] and base_list[0] != target_list[0]:
        raise OSError('Target is on a different drive to base. Target: '+target_list[0].upper()+', base: '+base_list[0].upper())

    # Starting from the filepath root, work out how much of the filepath is
    # shared by base and target.
    for i in range(min(len(base_list), len(target_list))):
        if base_list[i] != target_list[i]: break
    else:
        # If we broke out of the loop, i is pointing to the first differing path elements.
        # If we didn't break out of the loop, i is pointing to identical path elements.
        # Increment i so that in all cases it points to the first differing path elements.
        i+=1

    rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
    return os.path.join(*rel_list)

def extract_run_names(run_dirs, name_join_str='-', return_common=False):
   # No common names if only one item listed
   if len(run_dirs) == 1:
      return [os.path.basename(run_dirs[0])]
   
   common_part = os.path.commonprefix([ curr_dir.strip('/') for curr_dir in run_dirs])

   if not os.path.exists(common_part):
      common_dir = os.path.dirname(common_part)
   else:
      common_dir = common_part

   run_names = []
   for curr_dir in run_dirs:
      name_str = curr_dir.replace(common_dir, '')
      name_str = name_str.replace('//', '/').replace('./','')
      name_str = name_str.replace('/', name_join_str)
      name_str = name_str.strip(name_join_str)
      run_names.append(name_str)

   if return_common:
      return run_names, common_dir
   else:
      return run_names

def index_range_list(range_spec, max_value=None):

    if type(range_spec) is my_list_type:
        range_list = [ int(item) for item in range_spec ]
    elif range_spec.isdigit():
        range_list = [ int(range_spec) ]
    else:
        range_list = eval('range('+range_spec+')')

    if max_value != None:
       filt_range_list = []
       for range_val in range_list:
          if range_val < max_value:
             filt_range_list.append(range_val)
       range_list = filt_range_list

    return range_list

def evaluate_bool_str(boolStr, default=False):
    if boolStr == None:
        return default
    elif boolStr.lower() == "yes" or boolStr.lower() == "true":
        return True
    else:
        return False

def convert_timestamp_to_struct(timestamp_str):
   return datetime.datetime.strptime(timestamp_str, ASCII_TIMESTAMP_FORMAT).timetuple()

def BackwardsReader(file, BLKSIZE = 4096):
   """
Read a file  line by line, backwards
Stolen from:
http://code.activestate.com/recipes/439045-read-a-text-file-backwards-yet-another-implementat/"""
   if not hasattr(file, 'seek'):
      raise IOError('file argument must be a stream object')
  
   buf = ""

   try:
       # This is not legal in python 3 for text like files (see
       # https://stackoverflow.com/questions/21533391/seeking-from-end-of-file-throwing-unsupported-exception)
       # Think this is because we don't know the size of a character in
       # something like utf-8 unless we read forward.
       #file.seek(-1, 2)
       # But can do it this way
       file.seek(0, os.SEEK_END)
       file.seek(file.tell() - 1, os.SEEK_SET)
   except IOError:
      # File is of zero size
      return
   
   lastchar = file.read(1)
   trailing_newline = (lastchar == "\n")

   while 1:
      newline_pos = buf.rfind("\n")
      pos = file.tell()
      if newline_pos != -1:
         # Found a newline
         line = buf[newline_pos+1:]
         buf = buf[:newline_pos]
         if pos or newline_pos or trailing_newline:
            line += "\n"
         yield line
      elif pos:
         # Need to fill buffer
         toread = min(BLKSIZE, pos)
         file.seek(file.tell() - toread, os.SEEK_SET)
         buf = file.read(toread) + buf
         file.seek(file.tell() - toread, os.SEEK_SET)
         if pos == toread:
            buf = "\n" + buf
      else:
         # Start-of-file
         return
