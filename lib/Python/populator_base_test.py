from __future__ import print_function
from __future__ import absolute_import
from nose.tools import *
import os
import shutil
import filecmp
from .populator_base import *
import logging

test_data = os.path.dirname(__file__) + "/../../unit_test_data/"
expected_dir = test_data + "expected/create_run_scripts/"

class FakePopulator(PopulatorBase):
    def __init__(self):
        self.processing_dir = "create_run_scripts_test"
        PopulatorBase.__init__(self)
        os.environ["L2_SUPPORT_PATH"] = "/l2_support_fake_path"
        os.environ["L2_SUPPORT_UTIL_PATH"] = "/l2_support_fake_path/utils"
        os.environ["PATH"] = "/fake_bin_path"
        os.environ["PYTHONPATH"] = "/fake_python_path"
        os.environ["LD_LIBRARY_PATH"] = "/fake_lib_path"
        os.environ["LUA_PATH"] = "/l2_lua_fake_path"
        os.environ["abscodir"] = "/fake_absco_path"
        os.environ["merradir"] = "/fake_merra_path"
        self.l2_config_filename = "/fake_path/input/gosat/config/config.lua"
        self.l2_binary_filename = "/fake_path/l2_fp"
        self.abscoversion = "v4.2.0_unscaled"
        self.aggregate = False
        self.target_cluster = 'pbs_pro'
        self.array_idx_var_name = 'PBS_ARRAY_INDEX'
        self.job_array_option = '-J'
        self.depend_option = "-W depend=afterany"
        self.group_size = 1

    @property
    def variable_exports_list(self):
        return {'met_file': "/fake_path/ecmwf.h5",
                'imap_file': "/fake_path/imap.h5",
                'spectrum_file': "/fake_path/spectrum.h5"}
    @property
    def spectrum_file(self):
        return self.variable_exports_list["spectrum_file"]
    @property
    def additional_var_list(self):
        return { }
def test_create_run_scripts():
    '''Test the creation of the torque scripts'''
    shutil.rmtree("create_run_scripts_test", True)
    os.mkdir("create_run_scripts_test")
    original_env = os.environ.copy()
    try:
        pb = FakePopulator()

        # Create default SCF scripts
        pb.create_run_scripts(test_data + "gosat_tccon.config",
                              "input/GOSATFullPhysics/SoundingIds")

        # Create pleiades script
        pb.target_cluster = 'pleiades'
        pb.create_run_scripts(test_data + "gosat_tccon.config",
                              "input/GOSATFullPhysics/SoundingIds")
    finally:
        os.environ.clear()
        os.environ.update(original_env)
    for f in ["sounding_id.list", "aggregate.sh",
              "launch_jobs.sh", "l2_fp_job.sh", "pleiades_job.sh"]:
        print("Comparing file " + f)
        res = filecmp.cmp("create_run_scripts_test/" + f, expected_dir + f)
        if(not res):
            print("The comparison failed. You can run 'diff create_run_scripts_test/" + f + " " + expected_dir + f + "' for more info if needed")
        assert res
    

def test_read_id_list_file():
    '''Test the reading the ID list'''
    pb = PopulatorBase()
    id_list = pb.read_id_list_file(test_data + "gosat_tccon.config",
                                "input/GOSATFullPhysics/SoundingIds")
    assert id_list == ['20090827005603', '20090925130238', 
                       '20100223034944', '20100831023103', '20100914193918']
    
def test_read_map_values_file():
    # Need a test here, but don't actually know what this does so I couldn't
    # create one.
    pass

def test_get_config_keyword_value():
    pb = PopulatorBase()
    val = pb.get_config_keyword_value(test_data + "gosat_tccon.config",
                                      "input/InputProductFiles/L1BFile")
    assert (val ==
 "/data/smyth/Level2/test/tccon_small_set/acos_L1bB2900_tccon_5_good_qual.h5")


def test_read_id_list_file_large():
    '''Test the reading the ID list for a large file. Historically this has
    been really slow, so we have a test in place here to check the speed of 
    this.'''
    pb = PopulatorBase()
    id_list = pb.read_id_list_file(test_data + "large_sounding_ids.list")
    assert len(id_list) == 199158
    
def test_read_id_list_config_large():
    '''Test the reading the ID list for a large file. Historically this has
    been really slow, so we have a test in place here to check the speed of 
    this. 

    This checks a file that uses our old ASCII format, rather than the 
    simpler list of soundings'''
    pb = PopulatorBase()
    id_list = pb.read_id_list_file(test_data + "large.config",
                                   "input/OCOFullPhysics/SoundingIds")
    assert len(id_list) == 199158
    

