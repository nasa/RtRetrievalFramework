from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
import numpy
from matplotlib.pyplot import *

import scipy.stats as stats

from .routines_base import PlotRoutinesBase

BAND_NAMES = ("ABO2", "WCO2", "SCO2")
BAND_COLORS = ("b", "g", "r")
CHI2_DATASETS = ("reduced_chi_squared_o2_fph", "reduced_chi_squared_weak_co2_fph", "reduced_chi_squared_strong_co2_fph")

class FrankenburgPlotRoutinesReal(PlotRoutinesBase):

    def _chi2(self, **kwargs):
        chi2 = []
        for ds_name in CHI2_DATASETS:
            try:
                chi2.append( self.analysis_env.get_object_data(ds_name, **kwargs) )
            except ValueError:
                # Did not find dataset name
                chi2.append(None)

        return chi2

    def _plot_spectral_fit_wavelength(self, obj_names, sounding_id, num_colors_per_band, 
                                      SpectralParameters__measured_radiance, SpectralParameters__measured_radiance_uncert, SpectralParameters__modeled_radiance,
                                      SpectralParameters__wavelength,
                                      **kwargs):

        return self._plot_spectral_fit(obj_names, sounding_id, num_colors_per_band,
                                       SpectralParameters__measured_radiance, SpectralParameters__measured_radiance_uncert, SpectralParameters__modeled_radiance,
                                       SpectralParameters__wavelength,
                                       self._chi2(**kwargs),
                                       radiance_units="ph s$^{-1}$ m$^{-2}$ sr$^{-1}$ micron$^{-1}$",
                                       grid_units="$\mu$m (microns)",
                                       **kwargs)

    def _plot_spectral_fit_wavenumber(self, obj_names, sounding_id, num_colors_per_band,
                                      SpectralParameters__measured_radiance, SpectralParameters__measured_radiance_uncert, SpectralParameters__modeled_radiance,
                                      SpectralParameters__wavenumber,
                                      **kwargs):


        return self._plot_spectral_fit(obj_names, sounding_id, num_colors_per_band,
                                      SpectralParameters__measured_radiance, SpectralParameters__measured_radiance_uncert, SpectralParameters__modeled_radiance,
                                      SpectralParameters__wavenumber,
                                      self._chi2(**kwargs),
                                      radiance_units="W cm$^{-2}$ sr$^{-1}$ cm$^{-1}$",
                                      grid_units="wavenumber / cm$^{-1}$",
                                      **kwargs)

    def _plot_fit(self, x, y_meas, y_mod, x_label, y_label, axis, color='k', fontsize=10, leg_fontsize=8):
        axis.plot(x,y_meas, color+'-', x, y_mod, 'k--')
        axis.set_ylabel(y_label, fontsize=fontsize)
        axis.set_xlabel(x_label, fontsize=fontsize)
        axis.grid(True)

        leg = axis.legend( (r'measured', r'modeled'), loc=0, fontsize=leg_fontsize)
        leg.get_frame().set_alpha(0.5)

    def _plot_residual(self, x, y_meas, y_mod, meas_sigma, chi2, y_label, axis, color='k', fontsize=10, leg_fontsize=8):
        axis.plot(x, y_meas - y_mod, color+'-')
        axis.plot(x, meas_sigma,'k-', x, -meas_sigma, 'k-', alpha=0.5)
        axis.set_xlim(x[0], x[-1])
        axis.set_ylabel(y_label, fontsize=fontsize)
        axis.grid(True)

        leg1 = axis.legend(['meas-mod', 'meas err'], loc=0, fontsize=leg_fontsize)
        leg1.get_frame().set_alpha(0.5)

        # Remove labels for x axis of residual plot
        for label in axis.get_xticklabels():
            label.set_visible(False)

        stats_txt = '$\chi^2$=%.2f, SNR=%d' % (chi2, int(old_div(numpy.max(y_meas),numpy.mean(meas_sigma))))
        axis.set_xlabel(stats_txt)

    def _plot_histogram(self, y_meas, y_mod, meas_sigma, axis):
        axis.xaxis.set_major_formatter(NullFormatter())
        axis.yaxis.set_major_formatter(NullFormatter())

        nn, nbins, npatches = axis.hist(y_meas - y_mod, 30, normed=1, orientation='horizontal', facecolor='green', alpha=0.75)
        bincenters = 0.5*(nbins[1:]+nbins[:-1])
        y = mlab.normpdf(bincenters, 0,  np.mean(meas_sigma))
        l = axis.plot(y, bincenters, 'r-', linewidth=1)

    def _plot_qq(self, y_meas, y_mod, meas_sigma, left, width, axis):
        axis.xaxis.set_major_formatter(NullFormatter())
        axis.yaxis.set_major_formatter(NullFormatter())

        (osm, osr), (m, b, r) = stats.probplot(old_div((y_meas - y_mod), meas_sigma),  dist='norm')

        osmf = osm.take([0, -1])  # endpoints
        osrf = m * osmf + b       # fit line
        axis.plot(osm, osr, 'k.', osmf, osrf, 'r-')
        figtext(left+width+0.002, 0.58, 'slope=' + str(m)[0:4])

    def _plot_spectral_fit(self, obj_names, sounding_id, num_colors_per_band,
                           measured_radiance, measured_radiance_uncert, modeled_radiance,
                           grid_points, chi2,
                           radiance_units, grid_units,
                           **kwargs):
        """Creates one spectral residual plot per sounding per band.
        The plots show the residual differences, an over plot of the
        measured and modeled radiances, a histogram of residual
        differences and a QQ plot of the deviation of the residual
        differences from a gaussian.

        In the case of aggregated L2 files this routine could generate
        an enormous amount of plots which could in an interactive environment
        be very annoying. To guard against accidental execution on large
        data sets, this routine will only plot more than one sounding
        if the keyword argument all_soundings = True."""

        if len(sounding_id[0]) > 1 and not kwargs.get("all_soundings", False):
            raise Exception("In order to avoid crowded the screen with way too many plots in an interactive environment you must supply the keyword argument all_soundings=True before more than one sounding will be plotted.")

        left, width = 0.09, 0.785
        rect1 = [left, 0.73, width, 0.2]
        rect2 = [left, 0.1, width, 0.58]
        rect_hist = [left+width+0.005, 0.73, 0.1, 0.2]
        rect_qq = [left+width+0.005, 0.62, 0.1, 0.1]

        for obj_idx in range(len(obj_names)):
            for snd_idx, sounding_id in enumerate(sounding_id[obj_idx]):
                # Extract necessary input from the file
                num_colors = num_colors_per_band[obj_idx][snd_idx]

                meas = measured_radiance[obj_idx][snd_idx]
                meas_sigma = measured_radiance_uncert[obj_idx][snd_idx]
                mod = modeled_radiance[obj_idx][snd_idx]

                points = grid_points[obj_idx][snd_idx]

                band_plots = []
                for i, (band_name, color) in enumerate(zip(BAND_NAMES, BAND_COLORS)):
                    # Skip empty bands
                    if num_colors[i] == 0:
                        continue

                    band_chi2 = chi2[i][obj_idx][snd_idx]

                    idx_beg = sum(num_colors[:i])
                    idx_end = min(num_colors[i] + sum(num_colors[:i]) - 1, len(points) - 1)

                    band_points = points[idx_beg:idx_end]
                    band_meas = meas[idx_beg:idx_end]
                    band_mod = mod[idx_beg:idx_end]
                    band_sigma = meas_sigma[idx_beg:idx_end]

                    # Plot retrieval radiances overplot and residual
                    fig1 = figure()
                    band_plots.append(fig1)

                    ax_res = fig1.add_axes(rect1)  #left, bottom, width, height
                    ax_res.set_title("%s %s Spectral Fit" % (sounding_id, band_name))
                    ax_fit = fig1.add_axes(rect2, sharex=ax_res)

                    self._plot_fit(band_points, band_meas, band_mod, grid_units, 'Radiance / (%s)' % radiance_units, ax_fit, color) 
                    self._plot_residual(band_points, band_meas, band_mod, band_sigma, band_chi2, '$\Delta$ radiance', ax_res, color)

                    ax_hist = fig1.add_axes(rect_hist)
                    ax_hist.set_ylim(ax_res.get_ylim())
                    self._plot_histogram(band_meas, band_mod, band_sigma, ax_hist)

                    ax_qq = fig1.add_axes(rect_qq)
                    self._plot_qq(band_meas, band_mod, band_sigma, left, width, ax_qq)

                yield band_plots

class FrankenburgPlotRoutinesGuard(FrankenburgPlotRoutinesReal):

    def __new__(cls, analysis_env=None, **kwargs):
        """Make sure that these routines make sense to in the context of data loaded into environment"""

        if analysis_env != None:
            new_obj = super(FrankenburgPlotRoutinesGuard, cls).__new__(cls)
            if analysis_env.data_objs[0].get("SpectralParameters/wavelength", None) != None:
                new_obj.plot_spectral_fits = new_obj._plot_spectral_fit_wavelength
                return new_obj
            elif analysis_env.data_objs[0].get("SpectralParameters/wavenumber", None) != None:
                new_obj.plot_spectral_fits = new_obj._plot_spectral_fit_wavenumber
                return new_obj

        return None
