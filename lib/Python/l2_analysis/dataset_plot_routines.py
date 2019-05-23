from __future__ import absolute_import
import numpy
from matplotlib.pyplot import *

from .routines_base import PlotMaker
from .time_diff_plot_routines import TimeDiffPlotRoutines
from .sza_diff_plot_routines import SzaDiffPlotRoutines

BAD_VALUE =-9.99999000e+05

class DatasetPlotRoutines(PlotMaker, TimeDiffPlotRoutines, SzaDiffPlotRoutines):
    def __init__(self, **kwargs):
        TimeDiffPlotRoutines.__init__(self, **kwargs)

        self._create_dataset_routines(
            "xco2", "RetrievalResults/xco2", "XCO2 (ppm)", 
            translate_func=lambda ov: [ x * 1e6 for x in ov ], 
            diff_limits=(-2, 2))

        self._create_dataset_routines(
            "psurf", "RetrievalResults/surface_pressure_fph", "Surface Pressure Retrieved (mb)",
            translate_func=lambda ov: [ x * 1e-2 for x in ov ])

        self._create_dataset_routines(
            "psurf_ap", "RetrievalResults/surface_pressure_apriori_fph", "Surface Pressure A priori (mb)",
            translate_func=lambda ov: [ x * 1e-2 for x in ov ])

        self._create_dataset_routines(
            "temp_offset", "RetrievalResults/temperature_offset_fph", "Temperature Offset (k)")

        self._create_dataset_routines(
            "h2o_scale_factor", "/RetrievalResults/h2o_scale_factor", "H2O Scale Factor")

        self._create_dataset_routines(
            "abo2_chi2", "reduced_chi_squared_o2_fph", "O2 reduced chi^2")
        self._create_dataset_routines(
            "wco2_chi2", "reduced_chi_squared_weak_co2_fph", "Weak CO2 reduced chi^2")
        self._create_dataset_routines(
            "sco2_chi2", "reduced_chi_squared_strong_co2_fph", "Strong CO2 reduced chi^2")

        self._create_dataset_routines(
            "abo2_rms", "/SpectralParameters/relative_residual_mean_square_o2", "O2 RMS")
        self._create_dataset_routines(
            "wco2_rms", "/SpectralParameters/relative_residual_mean_square_weak_co2", "Weak CO2 RMS")
        self._create_dataset_routines(
            "sco2_rms", "/SpectralParameters/relative_residual_mean_square_strong_co2", "Strong CO2 RMS")

        self._create_dataset_routines(
            "abo2_albedo", "albedo_o2", "O2 Albedo")
        self._create_dataset_routines(
            "wco2_albedo", "albedo_weak_co2", "Weak CO2 Albedo")
        self._create_dataset_routines(
            "sco2_albedo", "albedo_strong_co2", "Strong CO2 Albedo")                    

        def rs_val_trans(orig_val):
            return [ rs[:,0] for rs in orig_val ]

        self._create_dataset_routines(
            "abo2_radiance_scaling", "RetrievalResults/radiance_scaling_o2", "O2 Radiance Scaling",
            translate_func=rs_val_trans)
        self._create_dataset_routines(
            "wco2_radiance_scaling", "RetrievalResults/radiance_scaling_weak_co2", "Weak CO2 Radiance Scaling",
            translate_func=rs_val_trans)
        self._create_dataset_routines(
            "sco2_radiance_scaling", "RetrievalResults/radiance_scaling_strong_co2", "Strong CO2 Radiance Scaling",
            translate_func=rs_val_trans)

        def rs_slope_trans(orig_val):
            return [ rs[:,1] for rs in orig_val ]

        self._create_dataset_routines(
            "abo2_radiance_scaling_slope", "RetrievalResults/radiance_scaling_o2", "O2 Radiance Scaling Slope",
            translate_func=rs_slope_trans)
        self._create_dataset_routines(
            "wco2_radiance_scaling_slope", "RetrievalResults/radiance_scaling_weak_co2", "Weak CO2 Radiance Scaling Slope",
            translate_func=rs_slope_trans)
        self._create_dataset_routines(
            "sco2_radiance_scaling_slope", "RetrievalResults/radiance_scaling_strong_co2", "Strong CO2 Radiance Scaling Slope",
            translate_func=rs_slope_trans)

        # Set any BAD values to nan so they are filtered out of the plot
        def filter_bad_values(value_sets):
            return [numpy.where(values == BAD_VALUE, numpy.nan, values) for values in value_sets]

        self._create_dataset_routines(
            "windspeed", "RetrievalResults/wind_speed", "Windspeed", translate_func=filter_bad_values)   

        aer_types_ds = "/Metadata/AerosolTypes"
        if self.analysis_env != None and self.analysis_env.data_objs[0].get(aer_types_ds, False):
            aer_types = self.analysis_env.data_objs[0].get_sounding_data(aer_types_ds)[0, :]
            for aer_idx, aer_name in enumerate(aer_types):
                if "retrieved_aerosol_aod_by_type" in list(self.analysis_env.data_objs[0]["/RetrievalResults"].keys()):
                    # So that scoping works out right we have to have a
                    # function return a function so aer_idx binds to the
                    # value for the current index, not the last one
                    def aer_filter_wrap(aer_idx):
                        def filter_aer_type(file_aer_data):
                            return [ aer_data[:, aer_idx] for aer_data in file_aer_data ]
                        return filter_aer_type

                    routine_name = "aerosol_%s" % aer_name.lower().decode('utf-8')
                    display_name = "Aerosol %s" % aer_name.decode('utf-8')
                    dataset_name = "RetrievalResults/retrieved_aerosol_aod_by_type"
                    self._create_dataset_routines(routine_name, display_name, display_name, translate_func=aer_filter_wrap(aer_idx), source_datasets=[dataset_name])
                else:
                    routine_name = "aerosol_%s" % aer_name.lower().decode('utf-8')
                    display_name = "Aerosol %s" % aer_name.decode('utf-8')
                    dataset_name = "RetrievalResults/aerosol_%s_aod" % aer_name.lower().decode('utf-8')
                    self._create_dataset_routines(routine_name, display_name, display_name, source_datasets=[dataset_name])
