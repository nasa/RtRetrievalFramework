from __future__ import absolute_import
from builtins import zip
from builtins import range
from copy import copy

from matplotlib import pyplot as plt
import numpy

from .decorators import call_data_pairs
from .ref_bar_plot_routines import RefBarPlotRoutines, ids2times
from .time_diff_plot_routines import TimeDiffPlotRoutines
from .filter_routines import FilterBase
from .status_routines import StatusBase

from scipy import stats
    
def r_squared(x,y, **kwargs):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return r_value**2

def least_squares_slope(x,y, sigma_filter=None, **kwargs):
    if sigma_filter != None:
        xco2_diff = x-y
        mean_diff = abs(numpy.average(xco2_diff))
        filtered = numpy.where(xco2_diff <= 2.5*mean_diff)
        x = x[filtered]
        y = y[filtered]
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return slope

class TCCONPlotRoutinesReal(RefBarPlotRoutines, TimeDiffPlotRoutines):

    def _compare_l2_tccon_data(self, l2_ids_list, l2_data_list, tccon_ids_list, tccon_data_list, comparison_func):
        "Perform some operation on each L2 set to compare with the TCCON set"
    
        if self.analysis_env == None:
            raise Exception("Must have access to analysis environment for matching L2 and TCCON soundings")

        results = []
        for l2_ids, l2_data, tccon_ids, tccon_data in zip(l2_ids_list, l2_data_list, tccon_ids_list, tccon_data_list):
            # Since we may be doing a second filtering, pass the L2 ids, TCCON we want indexes for
            corr_results = self.analysis_env.correlate_soundings([l2_ids, tccon_ids], [l2_ids], [tccon_ids])
            l2_match_indexes = corr_results['data_id_indexes']
            tccon_match_indexes = corr_results['addl_id_indexes']

            results.append( comparison_func(l2_data[l2_match_indexes], tccon_data[tccon_match_indexes]) )
        return results

    def _get_combined_values(self, obj_names, l2_ids, l2_xco2, tccon_ids, tccon_xco2, **kwargs):
        "Package TCCON XCO2 values with L2 values"
        
        tccon_ids_only = kwargs.get("tccon_ids_only", False)

        # Place TCCON values first since for multiple L2 files
        # the time diff code will always match to the first
        sounding_ids = [ tccon_ids[0] ]
        sounding_data = [ tccon_xco2[0] ]

        obj_names = [ "TCCON" ] + list(obj_names)

        # Scale L2 XCO2 to ppm
        l2_xco2_scaled = [ val * 1e6 for val in l2_xco2 ]
                    
        # Start with just L2 data
        if tccon_ids_only:
            sounding_ids += [ tccon_ids[0] for idx in range(len(l2_ids)) ]
            sounding_data += self._compare_l2_tccon_data(l2_ids, l2_xco2_scaled, tccon_ids, tccon_xco2, lambda x,y: x)
        else:
            sounding_ids += l2_ids
            sounding_data += l2_xco2_scaled

        return (obj_names, sounding_ids, sounding_data)
  
    def _plot_tccon_xco2_bar(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, comparison_func, **kwargs):
        l2_xco2_scaled = [ val * 1e6 for val in  RetrievalResults__xco2 ]
        
        xco2_comp_data = self._compare_l2_tccon_data(sounding_id, l2_xco2_scaled, SoundingHeader__sounding_id, TCCON__xco2_true, comparison_func)
        snd_times = ids2times(SoundingHeader__sounding_id)

        ax = self._plot_bar_stats(obj_names, snd_times, xco2_comp_data, **kwargs)

        return ax

    def plot_tccon_xco2_bar_mean(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (numpy.min(x)-1, numpy.max(x)+1))
        
        ax = self._plot_tccon_xco2_bar(obj_names,
                                       sounding_id,
                                       RetrievalResults__xco2,
                                       SoundingHeader__sounding_id,
                                       TCCON__xco2_true,
                                       lambda x,y: x-y,
                                       **kwargs)
        ax.set_title(ax.get_title().format(data_name="L2 - TCCON XCO2"))
        ax.set_ylabel("ppm")
        return ax
         
    def plot_tccon_xco2_bar_std(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (0, numpy.max(x)*1.10))
        
        ax = self._plot_tccon_xco2_bar(obj_names,
                                       sounding_id,
                                       RetrievalResults__xco2,
                                       SoundingHeader__sounding_id,
                                       TCCON__xco2_true,
                                       lambda x,y: x-y,
                                       stat_func=numpy.std, **kwargs)
        ax.set_title(ax.get_title().format(data_name="L2 - TCCON XCO2"))
        ax.set_ylabel("ppm")
        return ax

    def plot_tccon_xco2_bar_rsquared(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (0, numpy.max(x)*1.10))

        ax = self._plot_tccon_xco2_bar(obj_names,
                                       sounding_id,
                                       RetrievalResults__xco2,
                                       SoundingHeader__sounding_id,
                                       TCCON__xco2_true,
                                       r_squared,
                                       stat_func=None, **kwargs)
        ax.set_title(ax.get_title().format(data_name="R^2 L2 vs. TCCON XCO2"))
        return ax


    def plot_tccon_xco2_bar_lfit_slope(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):
        kwargs["ylim_func"] = kwargs.get("ylim_func", lambda x: (0, numpy.max(x)*1.10))
        sigma_filter = kwargs.get("sigma_filter", None)

        def slope_filtered(x,y):
            return least_squares_slope(x,y,sigma_filter)
        
        ax = self._plot_tccon_xco2_bar(obj_names,
                                       sounding_id,
                                       RetrievalResults__xco2,
                                       SoundingHeader__sounding_id,
                                       TCCON__xco2_true,
                                       slope_filtered,
                                       stat_func=None, **kwargs)
        ax.set_title(ax.get_title().format(data_name="Least Squares Slope L2 vs. TCCON XCO2"))
        return ax


    def plot_tccon_xco2_time(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):

        (obj_names, sounding_ids, sounding_data) = self._get_combined_values(obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs)

        ax = self._plot__dataset__time_multi("XCO2", "ppm", obj_names, sounding_ids, sounding_data, **kwargs)

        return ax

    def plot_tccon_xco2_diff_corr(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):

        (obj_names, sounding_ids, sounding_data) = self._get_combined_values(obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, tccon_ids_only=True)

        return self._plot__dataset__diff_corr("XCO2", "ppm", obj_names, sounding_ids, sounding_data, **kwargs)
        
    def get_tccon_l2_combined_data(self, obj_names, sounding_id, RetrievalResults__xco2, SoundingHeader__sounding_id, TCCON__xco2_true, **kwargs):
        """Returns TCCON and L2 XCO2 data packaged together and returned as namedtuple.
        Optional keyword 'tccon_ids_only = True' will return only L2 data for the same sounding ids"""

        snd_times = ids2times(sounding_ids)

        from collections import namedtuple
        Results = namedtuple('CombinedData', 'obj_names sounding_ids sounding_times sounding_data')
        return Results(obj_names, sounding_ids, snd_times, sounding_data)

class TCCONPlotRoutinesGuard(TCCONPlotRoutinesReal):
    def __new__(cls, analysis_env=None, **kwargs):
        """Make sure that these routines make sense to in the context of data loaded into environment"""

        # Only allow the object to be returned if the analysis envirionment contains
        # additional objects where the TCCON group is defined in the datasets
        if analysis_env != None and len(analysis_env.addl_objs) > 0 and 'TCCON' in list(analysis_env.addl_objs[0].keys()):
            return super(TCCONPlotRoutinesGuard, cls).__new__(cls)
        else:
            return None

class TCCONStatusFilterRoutines(FilterBase, StatusBase):
    
    def filter_tccon_site(self, sounding_id, SoundingHeader__sounding_id, TCCON__site_name, **kwargs):
        """Returns sounding ids which match a given TCCON site name"""

        # Make sure to filter ONLY on ids both in L2 and TCCON datasets
        common_ids = self.analysis_env.correlate_soundings(sounding_id + SoundingHeader__sounding_id)['comparison_ids']
        return self._filter_on_value([common_ids], [TCCON__site_name[0]], kwargs.get('by', None))

    def tccon_sites(self, TCCON__site_name, **kwargs):
        """Returns a count of each site name in the TCCON reference data file"""
        return self._count_items(TCCON__site_name)

