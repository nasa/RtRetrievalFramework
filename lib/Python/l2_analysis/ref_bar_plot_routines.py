from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import numpy
import scipy
from matplotlib.pyplot import *
import time

from . import utils

from .routines_base import PlotRoutinesBase

def ids2times(sounding_ids):
    return [ numpy.array([ time.mktime(dt.timetuple()) for dt in utils.id2time(snds) ]) for snds in sounding_ids ]

class RefBarPlotRoutines(PlotRoutinesBase):
    def _plot_bar_stats(self, obj_names, sounding_times, sounding_data, **kwargs):
        detrend = kwargs.pop("detrend", False)
        ylim_func = kwargs.pop("ylim_func", None)
        stat_func = kwargs.pop("stat_func", numpy.mean)

        if detrend:
            detrend_data = []
            detrend_slopes = []
            for time_vals, xco2_vals in zip(sounding_times, sounding_data):
                polycoeffs = scipy.polyfit(time_vals, xco2_vals, 1)
                # Calculate trend vals to remove
                yfit = scipy.polyval(polycoeffs, time_vals)
                ysub = yfit-yfit[0]

                detrend_data.append(xco2_vals - ysub)
                detrend_slopes.append(polycoeffs[0])
            sounding_data = detrend_data

        if stat_func != None:
            stat_data = [ stat_func(val) for val in sounding_data ]
        else:
            stat_data = sounding_data

        ax = self.new_axis(**kwargs)

        # Calculate where we want to put the xticks
        width = old_div(1.0,len(obj_names))
        xlocations = numpy.arange(len(stat_data))+width

        # Cycle colors for bars
        colors = [ cm.jet(1.0 * i/len(obj_names)) for i in range(len(obj_names)) ]

        # Plot stat values
        ax.bar(xlocations, stat_data, width=width, color=colors)

        # Use the object names for the xticks
        xticklocs = xlocations + old_div(width,2)
        xticks(xticklocs, obj_names)

        # Set some padding on the xaxis
        ax.set_xlim(0, xlocations[-1]+width*2)

        # Scale y axis to highlight smaller value changes better
        if ylim_func != None:
            ax.set_ylim( ylim_func(stat_data) )

        for midx, mval in enumerate(stat_data):
            val_offset = old_div((ax.get_ylim()[1]-ax.get_ylim()[0]),100)
            mloc = ( xticklocs[midx]-0.2, max(mval,0)+val_offset )
            mstr = "%7.3e" % mval

            if detrend:
                mstr += "\nw/slope:\n%.2e" % detrend_slopes[midx]
            ax.annotate(mstr, mloc)

        if stat_func != None:
            ax.set_title("%s {data_name}" % stat_func.__name__.title())
        else:
            ax.set_title("{data_name}")
            
        if detrend:
            ax.set_title(ax.get_title() + " (detrended)")

        return ax

    def plot_xco2_bar_mean(self, obj_names, sounding_id, RetrievalResults__xco2, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (numpy.min(x)-10, numpy.max(x)+10))
        
        xco2 = [ val * 1e6 for val in RetrievalResults__xco2 ]
        snd_times = ids2times(sounding_id)
        
        ax = self._plot_bar_stats(obj_names, snd_times, xco2, **kwargs)
        
        ax.set_title(ax.get_title().format(data_name="XCO2"))
        ax.set_ylabel("XCO2 (ppm)")

        return ax

    def plot_xco2_bar_std(self, obj_names, sounding_id, RetrievalResults__xco2, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (0, numpy.max(x)*1.20))
        
        xco2 = [ val * 1e6 for val in RetrievalResults__xco2 ]
        snd_times = ids2times(sounding_id)

        ax = self._plot_bar_stats(obj_names, snd_times, xco2, stat_func=numpy.std, **kwargs)
        
        ax.set_title(ax.get_title().format(data_name="XCO2"))
        ax.set_ylabel("ppm")

        return ax

    def plot_psurf_delta_bar_mean(self, obj_names, sounding_id, RetrievalResults__surface_pressure_fph, RetrievalResults__surface_pressure_apriori_fph, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (numpy.min(x)-2, numpy.max(x)+2))

        delta_psurf = [ (psurf_ap - psurf) * 1e-2 for psurf_ap,psurf in zip(RetrievalResults__surface_pressure_fph, RetrievalResults__surface_pressure_apriori_fph)  ]
        snd_times = ids2times(sounding_id)

        ax = self._plot_bar_stats(obj_names, snd_times, delta_psurf, **kwargs)
        
        ax.set_title(ax.get_title().format(data_name="Delta Surface Pressure (ap-ret)"))
        ax.set_ylabel("hPa")

        return ax

    def plot_psurf_delta_bar_std(self, obj_names, sounding_id, RetrievalResults__surface_pressure_fph, RetrievalResults__surface_pressure_apriori_fph, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (0, numpy.max(x)*1.10))

        delta_psurf = [ (psurf_ap - psurf) * 1e-2 for psurf_ap,psurf in zip(RetrievalResults__surface_pressure_fph, RetrievalResults__surface_pressure_apriori_fph)  ]
        snd_times = ids2times(sounding_id)

        ax = self._plot_bar_stats(obj_names, snd_times, delta_psurf, stat_func=numpy.std, **kwargs)
        
        ax.set_title(ax.get_title().format(data_name="Delta Surface Pressure (ap-ret)"))
        ax.set_ylabel("hPa")

        return ax
