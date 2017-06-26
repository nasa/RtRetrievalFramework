from __future__ import absolute_import
from matplotlib import pyplot as plt
import numpy

from .routines_base import PlotMaker
from .time_diff_plot_routines import TimeDiffPlotRoutines

class L1bPlotRoutinesReal(PlotMaker, TimeDiffPlotRoutines):
    def __init__(self, analysis_env=None, **kwargs):
        self.analysis_env = analysis_env
        self._create_dataset_routines("latitude", "FootprintGeometry/footprint_latitude", "Latitude (deg)")
        self._create_dataset_routines("longitude", "FootprintGeometry/footprint_longitude", "Longitude (deg)")
        self._create_dataset_routines("solar_zenith", "FootprintGeometry/footprint_solar_zenith", "Solar Zenith (deg)")
        self._create_dataset_routines("zenith", "FootprintGeometry/footprint_zenith", "Sounding Zenith (deg)")
        self._create_dataset_routines("solar_azimuth", "FootprintGeometry/footprint_solar_azimuth", "Solar Azimuth (deg)")
        self._create_dataset_routines("azimuth", "FootprintGeometry/footprint_azimuth", "Sounding Azimuth (deg)")

class L1bPlotRoutinesGuard(L1bPlotRoutinesReal):
    def __new__(cls, analysis_env=None, **kwargs):
        """Make sure that these routines make sense to in the context of data loaded into environment"""

        # Only allow the object to be returned if the analysis envirionment contains
        # additional objects where the L1B group is defined in the datasets
        if analysis_env != None and len(analysis_env.addl_objs) > 0 and 'FootprintGeometry' in list(analysis_env.addl_objs[0].keys()):
            return super(L1bPlotRoutinesGuard, cls).__new__(cls)
        else:
            return None
