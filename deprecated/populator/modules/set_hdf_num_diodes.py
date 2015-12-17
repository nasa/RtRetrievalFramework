import os
import re
import sys
from types import ListType
import contextlib
import logging

import numpy
import h5py

from OCO_Matrix import OCO_Matrix
from OCO_TextUtils import evaluate_bool_str
from Generator_Utils import *

RADIANCE_DATASETS = ('SoundingSpectra/radiance_o2', 'SoundingSpectra/radiance_weak_co2', 'SoundingSpectra/radiance_strong_co2')

DIODE_SECTION = 'ALGORITHMS'
DIODE_KEYWORD = 'num_diodes'

# Cache the diode sizes so we dont need to load the HDF file each time
used_l1b_filename = []
band_diode_list = []

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    if len(moduleSections) > 1:
        raise AttributeError('only one set_num_diodes allowed per file')

    l1b_file = Apply_Template(moduleSections[0].Get_Keyword_Value('l1b_file'), valuesDict, mapDict=mapDict)

    if not os.path.exists(l1b_file):
        raise OSError('Could not find L1B file: %s' % l1b_file)

    if len(band_diode_list) == 0:
        used_l1b_filename.append(l1b_file)
        with contextlib.closing(h5py.File(l1b_file, 'r')) as l1b_obj:
            for curr_dataset in RADIANCE_DATASETS:
                band_diode_list.append( str(l1b_obj[curr_dataset].shape[-1]) )
    elif used_l1b_filename[0] != l1b_file:
        raise Exception('L1B file for cached diode sizes: %s does not match current L1B file specified: %s' % (used_l1b_filename, l1b_file))

    source_obj = L2_Input.Input_File(source)

    diode_sect = source_obj.Get_Section(DIODE_SECTION)

    if len(diode_sect) == 0:
        raise LookupError('Could not find section: %s in file: %s' % (DIODE_SECTION, source))
    elif len(diode_sect) > 1:
        raise LookupError('Found too many sections name: %s in file: %s' % (DIODE_SECTION, source))

    diode_sect[0].Set_Keyword_Value(DIODE_KEYWORD, ' '.join(band_diode_list))

    source_obj.Write(destination)
