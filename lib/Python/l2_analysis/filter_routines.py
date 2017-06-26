from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import filter
from builtins import str
import sys
import re
import types

from .routines_base import RoutinesBase
from .status_routines import OUTCOME_MAPPING_DICT 

class FilterBase(RoutinesBase):
    def _filter_on_value(self, sounding_ids, data_values, filter_comparison=None, mapping_func=lambda x: x):
        if filter_comparison == None:
            print("Filter value not supplied returning all ids", file=sys.stderr)
            return sounding_ids

        if not isinstance(filter_comparison, types.FunctionType):
            def regexp_compare(val):
                return re.search(filter_comparison, str(val)) != None 
            comparison_func = regexp_compare
        else:
            comparison_func = filter_comparison

        ret_ids = []
        for obj_snd_ids, obj_values in zip(sounding_ids, data_values):
            obj_filtered_ids = []
            mapped_values = list(map(mapping_func, obj_values))
            for curr_id, curr_value in zip(obj_snd_ids, mapped_values):
                if hasattr(curr_value, "strip"):
                    curr_value = curr_value.strip()

                if comparison_func(curr_value):
                    obj_filtered_ids.append(curr_id)
            obj_filtered_ids.sort()
            ret_ids.append( tuple(obj_filtered_ids) )

        return tuple(ret_ids)
        
class FilterRoutines(FilterBase):
    def filter_surface_type(self, sounding_id, RetrievalResults__surface_type, **kwargs):
        """Returns sounding ids which match a given surface type supplied by the
        keyword argument 'by' which is used as a regular expression"""
        comparison_str = None
        # Make this filter easier to use as Lambertian is in both modes, add searching beginning
        # of string to make matching easier without user having to know regexp
        if ('by' in kwargs) != None:
            comparison_str = "^" + kwargs['by']
        return self._filter_on_value(sounding_id, RetrievalResults__surface_type, comparison_str)

    def filter_outcome(self, sounding_id, outcome, **kwargs):
        """Returns sounding ids which match a given outcome string supplied by the
        keyword argument 'by' which is used as a regular expression"""
        mapping_func = lambda x: OUTCOME_MAPPING_DICT[x]
        return self._filter_on_value(sounding_id, outcome, kwargs.get('by', None), mapping_func)
    
    def filter_solar_zenith(self, sounding_id, SoundingGeometry__sounding_solar_zenith, **kwargs):
        """Returns sounding ids which match the comparison function given in the 'by' keyword applied the 
        L1B solar zenith angle"""
        return self._filter_on_value(sounding_id, SoundingGeometry__sounding_solar_zenith, kwargs.get('by', None))

    def filter_cloud_free(self, sounding_id, cloud_flag, **kwargs):
        """Returns sounding ids which match are marked as cloud free by the Cloud file"""
        comparison = lambda x: x == 0
        return self._filter_on_value(sounding_id, cloud_flag, comparison)

    def filter_snd_pos(self, sounding_id, **kwargs):
        pos_index = kwargs.get('by', None)
        def is_snd_pos(pos):
            def is_pos_helper(snd_id):
                if str(snd_id)[-1] == str(pos):
                    return True
                else:
                    return False
            return is_pos_helper

        ret_ids = []
        for obj_ids in sounding_id:
            ret_ids += list(filter(is_snd_pos(pos_index), obj_ids))
        return ret_ids
