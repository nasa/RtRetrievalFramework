from __future__ import absolute_import
from builtins import zip
import numpy
import re

from .routines_base import RoutinesBase

class DataAccessRoutines(RoutinesBase):
    def get_data(self, _data_name=None):
        return self.analysis_env.get_object_data(_data_name)

    def _get_sv_values(self, sv_elem_name_re, sv_data):
        sv_names = self.analysis_env.get_object_data("state_vector_names")

        sv_filtered = []
        for obj_sv_names, obj_sv_data in zip(sv_names, sv_data):
            obj_filtered = []
            for snd_names, snd_data in zip(obj_sv_names, obj_sv_data):
                snd_filtered = []
                for name, value in zip(snd_names, snd_data):
                    if re.search(sv_elem_name_re, name):
                        snd_filtered.append(value)
                obj_filtered.append(snd_filtered)
            sv_filtered.append( numpy.array(obj_filtered) )

        return sv_filtered

    def get_sv_names(self, _sv_elem_name_re):
        sv_names = self.analysis_env.get_object_data("state_vector_names")
        return self._get_sv_values(_sv_elem_name_re, sv_names)

    def get_sv_result(self, _sv_elem_name_re):
        sv_results = self.analysis_env.get_object_data("state_vector_result")
        return self._get_sv_values(_sv_elem_name_re, sv_results)

    def get_sv_apriori(self, _sv_elem_name_re):
        sv_apriori = self.analysis_env.get_object_data("state_vector_apriori")
        return self._get_sv_values(_sv_elem_name_re, sv_apriori)

