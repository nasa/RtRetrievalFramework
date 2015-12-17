#!/usr/bin/env python

# Load standard modules
import os
import sys
from optparse import OptionParser

import bisect
import re

# Load additional modules
import numpy
import matplotlib.pyplot as pyplot

# Load L2 modules
import L2_Version_Util
import L2_Plot_Util
import OCO_TextUtils
from OCO_Matrix import OCO_Matrix

DICT_KEY_COLUMN = 'Case'
FILL_VALUE = -999.99

#X_AXIS = ('LATITUDE', 'Latitude')


ID_ONLY_RE = r'\d{14}'
ID_POL_RE  = ID_ONLY_RE + r'[PS]?'
PLOT_FILTERS = ( (ID_POL_RE + r'_ABO2$',      'A band'),
                 (ID_POL_RE + r'_WCO2$',      'Weak CO2'),
                 (ID_POL_RE + r'_ABO2_WCO2$', 'A band + Weak CO2'),
                 (ID_POL_RE + r'_SCO2$',      'Strong CO2'),
                 (ID_POL_RE + r'_3BND$',      'Three Band'), )
LEGEND_VALUES = ( (ID_ONLY_RE + r'P', 'P polarization'),
                  (ID_ONLY_RE + r'S', 'S polarization'),
                  (ID_ONLY_RE, '(P+S)/2'),
                  )
PLOT_TYPES = ( ('XTARG_RET', 'XCO2', 1e6),
               ('RET_PSURF', 'Surface Pressure', 1e-2),
               ('RET_ALL_TOTAL_AOD', 'Aerosol Optical Depth', 1),
               )

#X_AXIS = ('Case', 'Sounding ID', '('+ID_ONLY_RE+')')
X_AXIS = ('Case', 'Index', '('+ID_ONLY_RE+')', 1)

def plot_gosat_summ(summary_files):

    if len(summary_files) > 1:
        summary_names, common_part = OCO_TextUtils.extract_run_names(summary_files, return_common=True)
    else:
        common_part = summary_files[0]
    
    plot_data = {}
    for curr_summ_file in summary_files:
        file_obj = OCO_Matrix(curr_summ_file, as_strings=True)

        plot_data = file_obj.as_dict(DICT_KEY_COLUMN, existing_dict=plot_data)

    figure_index = 0
    x_axis_index = []
    for curr_yaxis in PLOT_TYPES:
        for curr_filter in PLOT_FILTERS:
            legend_names = []
            plot_objs = []

            fig = None
            curr_axis = None

            for leg_idx, curr_legend in enumerate(LEGEND_VALUES):
                # Extract plot data
                x_data = []
                y_data = []
                for scene_name, scene_data in plot_data.items():
                    if not re.search(curr_filter[0], scene_name) or not re.search(curr_legend[0], scene_name):
                        continue

                    if not scene_data.has_key(X_AXIS[0]) or not scene_data.has_key(curr_yaxis[0]):
                        continue

                    print 'Using scene: ', scene_name, 'for', curr_filter[1], curr_legend[1]

                    x_value = scene_data[X_AXIS[0]]
                    y_value = scene_data[curr_yaxis[0]]

                    if len(X_AXIS) >= 3:
                        axis_match = re.search(X_AXIS[2], x_value)
                        x_value = axis_match.group()

                    if len(X_AXIS) >= 4:
                        if not x_value in x_axis_index:
                            x_axis_index.append(x_value)
                            
                        x_value = x_axis_index.index(x_value)                            
                        
                    try:
                        x_value = float(x_value)
                    except:
                        pass

                    y_value = float(y_value)

                    if abs(y_value - FILL_VALUE) > 1e-2:
                        y_value *= curr_yaxis[2]
                        insert_point = bisect.bisect(x_data, x_value)

                        x_data.insert(insert_point, x_value)
                        y_data.insert(insert_point, y_value)

                if len(x_data) > 0:
                    if fig == None:
                        fig = pyplot.figure(figure_index)
                        fig.clf()
                        curr_axis = pyplot.subplot(111)
                    
                    curr_plot_obj = pyplot.plot(x_data, y_data, 'x')
                    plot_objs.append( curr_plot_obj )
                    
                    legend_names.append( curr_legend[1] )

            if len(plot_objs) > 0:
                figure_index += 1
                
                curr_axis.set_xlabel(X_AXIS[1])
                curr_axis.set_ylabel(curr_yaxis[1])
                curr_axis.set_title('%s -- %s' % (curr_filter[1], common_part))

                leg = curr_axis.legend( plot_objs,
                                        legend_names,
                                        loc='lower right',
                                        borderpad = 0.1,
                                        labelspacing = 0.05,
                                        borderaxespad = 0.2,
                                        )


    pyplot.show()


def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [summary_files]...")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) <= 0:
        parser.error('summary files for plotting are required')

    summary_files = args
       
    plot_gosat_summ(summary_files)

if __name__ == "__main__":
    standalone_main()


