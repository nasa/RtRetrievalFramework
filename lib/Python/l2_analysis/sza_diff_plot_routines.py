from __future__ import absolute_import
from builtins import zip
import numpy
from matplotlib.pyplot import *

from .routines_base import PlotRoutinesBase, PlotMaker

from .decorators import call_data_pairs, add_optional_datasets

class SzaDiffPlotRoutines(PlotRoutinesBase):

    def _plot_sza_vs_value(self, sza_values, plot_values, get_lines=None, **kwargs):
        ax = self.new_axis(**kwargs)
        error = kwargs.get("error", [])

        markers = Line2D.filled_markers[:len(plot_values)]
        for idx, (x_val, p_val, marker) in enumerate(zip(sza_values, plot_values, markers)):
            if idx < len(error) and error[idx] != None:
                lines = ax.errorbar(s_val, p_val, yerr=error[idx], marker=marker, linestyle='')
                # errorbar routine creates multiple lines which can screw up
                # legends if you use the default legend with jut names call method
                if get_lines != None: get_lines.append(lines[0]) 
            else:
                line = ax.plot(x_val, p_val, marker)
                if get_lines != None: get_lines.append(line)

        return ax

    @add_optional_datasets("sounding_solar_zenith")
    def _plot__dataset__sza(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, sounding_solar_zenith=None, **kwargs):

        if sounding_solar_zenith != None:
            sza = [ s[:,0] for s in sounding_solar_zenith ]
        else:
            raise Exception("No solar zenith data available")

        lines = []
        ax = self._plot_sza_vs_value(sza, dataset_values, get_lines=lines, **kwargs)
        
        # Set plot labels
        ax.set_title("%s by SZA" % dataset_name)
        ax.set_xlabel("Solar Zenith (deg)")
        ax.set_ylabel(value_name)

        # Use "best" location so we dont place box over data
        if len(lines) > 1:
            leg = ax.legend(obj_names, loc='best')

        return ax

    @add_optional_datasets("sounding_solar_zenith")
    @call_data_pairs((3,4,5))
    def _plot__dataset__diff_sza(self, dataset_name, value_name, obj_names, sounding_ids, dataset_values, sounding_solar_zenith=None, **kwargs):
        """Plots the difference for each sounding versus solar zenith"""

        if sounding_solar_zenith != None:
            sza = [ s[:,0] for s in sounding_solar_zenith ]
        else:
            raise Exception("No solar zenith data available")

        dataset_diff = dataset_values[0] - dataset_values[1]

        ax = self._plot_sza_vs_value(sza, [dataset_diff], **kwargs)
        self._add_diff_stat_text(ax, dataset_diff)
        
        # Set plot labels
        ax.set_title("%s absolute difference:\n%s - %s " % (dataset_name, obj_names[0], obj_names[1]))
        ax.set_xlabel("Solar Zenith (deg)")
        ax.set_ylabel("%s difference" % value_name)

        return ax
