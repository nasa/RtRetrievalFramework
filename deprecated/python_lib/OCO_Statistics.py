#!/usr/bin/env python

import sys
import math
import numpy

from OCO_MathUtil import *

def Residual_RMS(array1, array2):

    abs_residual = (array1 - array2)
    rel_residual = numpy.where(array2 == 0.0, abs_residual/1e-20, abs_residual/array2)

    if (len(rel_residual) == 0):
        rms = 0
    else:
        rms = numpy.sqrt(numpy.mean(rel_residual * rel_residual))

    return (abs_residual, rms)

def Ratio_RMS(array1, array2):

    residual = numpy.where(array2 == 0.0, array1/1e-20, array1/array2)
        
    if (len(residual) == 0):
        rms = 0
    else:
        rms = numpy.sqrt(numpy.mean(residual * residual))
   
    return (residual, rms)
