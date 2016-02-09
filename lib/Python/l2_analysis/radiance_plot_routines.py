from __future__ import absolute_import
from builtins import zip
from builtins import range
import numpy
from matplotlib.pyplot import *

from .routines_base import PlotMaker
from .time_diff_plot_routines import TimeDiffPlotRoutines

class RadiancePlotRoutines(PlotMaker, TimeDiffPlotRoutines):
    def __init__(self, **kwargs):
        TimeDiffPlotRoutines.__init__(self, **kwargs)

        def mean_radiance(rad, colors, band_idx):
            res = []
            for snd_idx in range(colors.shape[0]):
                beg_idx = band_idx > 0 and numpy.sum(colors[snd_idx,:band_idx]) or 0
                end_idx = band_idx < colors.shape[0] and numpy.sum(colors[snd_idx,:band_idx+1]) or numpy.sum(colors[snd_idx,:])
                res.append( numpy.mean(rad[snd_idx,beg_idx:end_idx]) )
            return numpy.array(res)

        for band_idx, band_name in enumerate(("abo2", "wco2", "sco2")):
            def mean_data_rad_band(band_idx):
                return lambda rad, colors: [ mean_radiance(r,c,band_idx) for r,c in zip(rad, colors) ]

            self._create_dataset_routines("%s_rad_mean" % band_name,
                                          "Mean /SpectralParameters/modeled_radiance %s" % band_name.upper(), 
                                          "Mean %s Modeled Radiance" % band_name.upper(), 
                                          source_datasets=["/SpectralParameters/modeled_radiance","num_colors_per_band"],
                                          translate_func=mean_data_rad_band(band_idx))

            self._create_dataset_routines("%s_rad_mean_uncert" % band_name,
                                           "Mean /SpectralParameters/measured_radiance_uncert %s" % band_name.upper(), 
                                          "Mean %s Measured Radiance Uncertainty" % band_name.upper(), 
                                          source_datasets=["/SpectralParameters/measured_radiance_uncert","num_colors_per_band"],
                                          translate_func=mean_data_rad_band(band_idx))
            

