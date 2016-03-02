from __future__ import absolute_import
from builtins import zip
from .decorators import call_data_pairs

import numpy
from matplotlib.pyplot import *

from .routines_base import PlotMaker
from .time_diff_plot_routines import TimeDiffPlotRoutines

class AlbedoPlotRoutines(PlotMaker, TimeDiffPlotRoutines):
    
    def _plot_albedo_ap_ret_corr(self, band_name, albedo_apriori, albedo_retrieved, sounding_id, **kwargs):
        dataset_name = "%s Albedo Retrieved vs. Apriori" % band_name
        value_name = "%s Albedo" 
        obj_names = ["Offset Apriori", "Offset Retrieved", ]

        sounding_ids = []
        dataset_values = []
        for ap_offset, ret_offset, snd_id in zip(albedo_apriori, albedo_retrieved, sounding_id):
            dataset_values.append(ap_offset)
            dataset_values.append(ret_offset)
            sounding_ids.append(snd_id)
            sounding_ids.append(snd_id)
        return self._plot__dataset__diff_corr(dataset_name, value_name, obj_names, sounding_ids, dataset_values, **kwargs)

    def plot_abo2_albedo_ap_ret_corr(self, albedo_apriori_o2_fph, albedo_o2_fph, sounding_id, **kwargs):
        return self._plot_albedo_ap_ret_corr("O2", albedo_apriori_o2_fph, albedo_o2_fph, sounding_id, **kwargs)

    def plot_wco2_albedo_ap_ret_corr(self, albedo_apriori_weak_co2_fph, albedo_weak_co2_fph, sounding_id, **kwargs):
        return self._plot_albedo_ap_ret_corr("Weak CO2", albedo_apriori_weak_co2_fph, albedo_weak_co2_fph, sounding_id, **kwargs)

    def plot_sco2_albedo_ap_ret_corr(self, albedo_apriori_strong_co2_fph, albedo_strong_co2_fph, sounding_id, **kwargs):
        return self._plot_albedo_ap_ret_corr("Strong CO2", albedo_apriori_strong_co2_fph, albedo_strong_co2_fph, sounding_id, **kwargs)
