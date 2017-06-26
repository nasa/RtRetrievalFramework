from __future__ import absolute_import
from builtins import zip
import time
from copy import copy
import matplotlib.dates

import numpy
from matplotlib.pyplot import *
from matplotlib.dates import DateFormatter

from .decorators import call_data_pairs, add_optional_datasets
from . import utils

from .routines_base import PlotRoutinesBase, PlotMaker

DATETIME_FORMAT = "%Y-%m-%d %H:%M:%S" 

class DateLookupFormatter(Formatter):
    """Formatter for dates in PlotRoutines plots. Used to eliminate
    gaps on axis due to time gaps."""
    def __init__(self, dates, fmt=DATETIME_FORMAT, **kwargs):
        self.dates = dates
        self.fmt = fmt

    def __call__(self, x, pos=0):
        "Return the label for time x at position pos"
        ind = int(round(x))
        if ind>=len(self.dates) or ind<0: return ''
        return self.dates[ind].strftime(self.fmt)

class TimeDiffPlotRoutines(PlotRoutinesBase):
    """Defines helper functions for plotting time and diff plots"""
    
    def _plot_id_vs_value(self, sounding_ids, plot_values, get_lines=None, time_strings=None, **kwargs):
        ax = self.new_axis(**kwargs)
        error = kwargs.get("error", [])
        show_gaps = kwargs.get("show_gaps", True)

        markers = Line2D.filled_markers[:len(plot_values)]
        for idx, (snd_ids, p_val, marker) in enumerate(zip(sounding_ids, plot_values, markers)):
            # Fall back to converting sounding ids into times if time strings are not supplied 
            if time_strings != None:
                datetimes = utils.time_string_to_dt(time_strings)
            else:
                datetimes = utils.id2time(snd_ids)
            
            if show_gaps:
                ax.xaxis.set_major_formatter(DateFormatter(DATETIME_FORMAT))
                x_val = datetimes
            else:
                ax.xaxis.set_major_formatter(DateLookupFormatter(datetimes, **kwargs))
                x_val = numpy.arange(len(p_val))

            if idx < len(error) and error[idx] != None:
                lines = ax.errorbar(x_val, p_val, yerr=error[idx], marker=marker, linestyle='')
                # errorbar routine creates multiple lines which can screw up
                # legends if you use the default legend with jut names call method
                if get_lines != None: get_lines.append(lines[0]) 
            else:
                line = ax.plot(x_val, p_val, marker)
                if get_lines != None: get_lines.append(line)

        # Set some padding some padding on the xaxis
        if not show_gaps:
            ax.set_xlim(-10, len(plot_values[0])+10)
        
        ax.get_figure().autofmt_xdate()
  
        return ax

    def _plot_exceeded_diffs(self, sounding_ids, diff_values, **kwargs):
        diff_limits = kwargs.get("diff_limits", None)

        if diff_limits == None:
            fig = self._plot_id_vs_value((sounding_ids,), (diff_values,), **kwargs)
            return ax

        smaller = numpy.where(diff_values < diff_limits[0])
        larger = numpy.where(diff_values > diff_limits[1])
        limit_exceed = copy(diff_values)
        limit_exceed[:] = numpy.nan
        limit_exceed[smaller] = diff_limits[0]
        limit_exceed[larger] = diff_limits[1]

        num_exceeded = len(smaller[0])+len(larger[0])
        
        # Plot DIFF_VALUES vs sounding id time
        ax = self._plot_id_vs_value((sounding_ids,sounding_ids), (diff_values,limit_exceed), **kwargs)

        ax.set_ylim(diff_limits)

        if num_exceeded > 0:
            ax.legend( ("Within limits (%d)" % (len(diff_values)-num_exceeded),
                        "Exceeds limits (%d)" % (num_exceeded)) )
        return ax

    @call_data_pairs((3,4,5))
    def _plot__dataset__diff_corr(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs):
        """Plots two dataset where one set of values is the x-axis and the other the y-axis.
        If the values for all soundings match then they will all fall on the one-to-one line."""

        ax = self.new_axis(**kwargs)
        
        # Plot one to one line
        valid_idx_0 = numpy.where(numpy.isnan(dataset_values[0]) == False)
        valid_idx_1 = numpy.where(numpy.isnan(dataset_values[1]) == False)

        all_min = min(numpy.min(dataset_values[0][valid_idx_0]), numpy.min(dataset_values[1][valid_idx_1]))
        all_max = max(numpy.max(dataset_values[0][valid_idx_0]), numpy.max(dataset_values[1][valid_idx_1]))
        ax.plot([all_min, all_max], [all_min, all_max], '--')

        # Plot correlated values
        ax.plot(dataset_values[0], dataset_values[1], 'D')

        dataset_diff = dataset_values[0] - dataset_values[1]
        self._add_diff_stat_text(ax, dataset_diff)

        ax.set_title("%s difference correlation" % (dataset_name))
        ax.set_xlabel("%s %s" % (value_name, obj_names[0]))
        ax.set_ylabel("%s %s" % (value_name, obj_names[1]))

        return ax

    @call_data_pairs((3,4,5))
    def _plot__dataset__diff_hist(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs):
        """Plots the a histogram showing the number of the binned differences"""
        
        ax = self.new_axis(**kwargs)

        dataset_diff = dataset_values[0] - dataset_values[1]

        ax.hist(dataset_diff)
        
        self._add_diff_stat_text(ax, dataset_diff)

        ax.set_title("%s difference histogram" % dataset_name)
        ax.set_xlabel("%s difference %s - %s" % (value_name, obj_names[0], obj_names[1]))
        ax.set_ylabel("Count")

        return ax

    @call_data_pairs((3,4,5))
    @add_optional_datasets("FtsRunLog/time_string")
    def _plot__dataset__diff_time_fts(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, FtsRunLog__time_string=None, **kwargs):
        """Plots the difference for each sounding versus the sounding time"""

        time_strings = None
        if FtsRunLog__time_string:
            # Extract only 1st band's timestamps from the first file
            time_strings = FtsRunLog__time_string[0][:, 0]
            
        dataset_diff = dataset_values[0] - dataset_values[1]

        ax = self._plot_id_vs_value(sounding_ids, [dataset_diff], time_strings=time_strings, **kwargs)
        self._add_diff_stat_text(ax, dataset_diff)
        
        # Set plot labels
        ax.set_title("%s absolute difference:\n%s - %s " % (dataset_name, obj_names[0], obj_names[1]))
        ax.set_xlabel("Sounding Time")
        ax.set_ylabel("%s difference" % value_name)

        return ax

    @call_data_pairs((3,4,5))
    def _plot__dataset__diff_time(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs):
        return self._plot__dataset__diff_time_fts(dataset_name, value_name, obj_names, sounding_ids, dataset_values, None, **kwargs)


    @add_optional_datasets("FtsRunLog/time_string")
    def _plot__dataset__time(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, FtsRunLog__time_string=None, **kwargs):
        """Plots multiple sets dataset values (multiple files) on the y-axis versus the
        sounding time (sounding id) on the x-axis"""

        time_strings = None
        if FtsRunLog__time_string:
            # Extract only 1st band's timestamps from the first file
            time_strings = FtsRunLog__time_string[0][:, 0] 

        # Plot data vs sounding id time
        lines = []
        ax = self._plot_id_vs_value(sounding_ids, dataset_values, get_lines=lines, time_strings=time_strings, **kwargs)
        
        # Set plot labels
        ax.set_title("%s by time" % dataset_name)
        ax.set_xlabel("Sounding Time")
        ax.set_ylabel(value_name)

        # Use "best" location so we dont place box over data
        if len(lines) > 1:
            leg = ax.legend(obj_names, loc='best')

        return ax

    @call_data_pairs((3,4,5))
    def _plot__dataset__time_multi(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs):
        return self._plot__dataset__time(dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs)
