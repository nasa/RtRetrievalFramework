from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import range
from collections import defaultdict

from .routines_base import RoutinesBase

OUTCOME_MAPPING_DICT = {0: "error", 1: "converged", 2: "converged", 3: "max_iter", 4: "max_div"}

class StatusBase(RoutinesBase):
    def _count_items(self, data_values, obj_names=None, mapping_func=lambda x: x):
        if obj_names == None:
            data_names = list(range(len(data_values)))
        else:
            data_names = obj_names
        
        tot_counts = {}
        for o_name, o_values in zip(data_names, data_values):
            tot_counts[o_name] = defaultdict(int)
            mapped_values = list(map(mapping_func, o_values))
            for count_val in mapped_values:
                if hasattr(count_val, "strip"):
                    tot_counts[o_name][count_val.strip()] += 1
                else:
                    tot_counts[o_name][count_val] += 1
            # Convert to normal dict, so can be pretty printed better
            tot_counts[o_name] = dict(tot_counts[o_name])

        if obj_names == None:
            return list(tot_counts.values())
        else:
            return tot_counts

class StatusRoutines(StatusBase):
    def outcome(self, obj_names, outcome, **kwargs):
        """Returns a count of outcome strings per set of data"""
        mapping_func = lambda x: OUTCOME_MAPPING_DICT[x]
        counts = self._count_items(outcome, obj_names, mapping_func)
        return counts

    def surface_type(self, obj_names, surface_type, **kwargs):
        """Returns a count of surface type strings per set of data"""
        return self._count_items(surface_type, obj_names)

    def cloud_flag(self, obj_names, cloud_flag, **kwargs):
        """Returns counts of cloud flag values"""
        return self._count_items(cloud_flag, obj_names)
