from __future__ import absolute_import
###################################################################
# Classes can be placed in this module to be picked up and used by
# analyze_l2.py
# Any new class must be named like: *Routines. Also the class must
# be listed below to beand be imported into this module's context.
#
# Classes should inherit from one of the base_routines classes.
#
# All methods of the class will be picked up as analysis routines
# as long as its name does not start with the _ character.
###################################################################

# Import objects into analysis context
from .analysis_environment import AnalysisEnvironment
from .status_routines import StatusRoutines
from .custom_plot_routines import CustomPlotRoutines
from .dataset_plot_routines import DatasetPlotRoutines
from .filter_routines import FilterRoutines
from .frankenburg_plot_routines import FrankenburgPlotRoutinesGuard as FrankenburgPlotRoutines
from .radiance_plot_routines import RadiancePlotRoutines
from .ref_bar_plot_routines import RefBarPlotRoutines
from .status_routines import StatusRoutines
from .tccon_plot_routines import TCCONPlotRoutinesGuard as TCCONPlotRoutines, TCCONStatusFilterRoutines
from .albedo_plot_routines import AlbedoPlotRoutines
from .l1b_plot_routines import L1bPlotRoutinesGuard as L1bPlotRoutines
from .data_access_routines import DataAccessRoutines
from .radiance_scaling_plots import RadianceScalingRoutines
