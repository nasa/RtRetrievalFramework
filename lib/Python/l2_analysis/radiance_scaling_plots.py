from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
from matplotlib import pyplot as plt
import numpy

from .time_diff_plot_routines import TimeDiffPlotRoutines
from .data_access_routines import DataAccessRoutines
from .filter_routines import FilterRoutines

class RadianceScalingRoutines(DataAccessRoutines, TimeDiffPlotRoutines):

    def plot_radiance_scaling(self, obj_names, sounding_id, **kwargs):
        if len(obj_names) > 1:
            raise NotImplementedError("More than one input file not yet supported")
        ax = kwargs.pop('axis', self.new_axis(**kwargs))
        band_list = ["A-Band", "WC-Band", "SC-Band"] 
        datasets = ["/RetrievalResults/radiance_scaling_o2", "/RetrievalResults/radiance_scaling_weak_co2", "/RetrievalResults/radiance_scaling_strong_co2"]
        for band_name, ds_name in zip(band_list, datasets):
            # Pick off only the constant term from the first object
            rad_scaling = [self.analysis_env.get_object_data(ds_name, ids=sounding_id, **kwargs)[0][:,0]]
            self._plot__dataset__time("Radiance Scaling", "scale factor (unitless)", obj_names, sounding_id, rad_scaling, axis=ax, **kwargs)
        plt.legend(band_list, 0)
        return ax

    def plot_radiance_scaling_by_pos(self, obj_names, sounding_id, **kwargs):
        fig = plt.figure()
        filt = FilterRoutines()
        positions = kwargs.get("positions", list(range(1, 9))) 
        res = []
        # Remove the ids argument if passed, sounding_id would already be filtered
        _ = kwargs.pop("ids", None)
        subplot_idx = 1
        for pos in positions:
            pos_ids = filt.filter_snd_pos(sounding_id, by=pos)
            axis = fig.add_subplot(2, old_div(len(positions),2), subplot_idx)
            subplot_idx += 1
            self.plot_radiance_scaling(obj_names, [pos_ids], axis=axis, **kwargs)
            axis.set_title("Radiance Scaling Pos: %d" % pos)
            res.append(axis)
        return res
