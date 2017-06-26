from __future__ import absolute_import
from nose.tools import *
import os
from .l2_input import *
test_data = os.path.dirname(__file__) + "/../../unit_test_data/"

def test_get_all_section_names():
    '''Make sure we can get all the section names from a XML file.'''
    config = L2InputFile(test_data + "acos_Fph_090627_37_Test_v050050_Fph2900_r11_110523154640i.config")
    assert config.get_section("input")[0].get_all_section_names() == \
        ['input', 'input->PrimaryExecutable', 'input->Geometry', 
         'input->LogMetadata', 'input->LogMetadata->LIST', 
         'input->LogMetadata->LIST->VALUES', 'input->InputProductFiles', 
         'input->SCFIdentification', 'input->JobIdentification', 
         'input->GOSATFullPhysics', 'input->GOSATFullPhysics->LIST', 
         'input->GOSATFullPhysics->LIST->VALUES', 
         'input->RecordedAuxiliaryInputFiles', 
         'input->PGENameGroup', 'input->StaticFileIdentificationFiles', 
         'input->ProductPathGroup', 'input->GOSATIdentification', 
         'input->DynamicAuxiliaryInputFiles', 'input->MonitorGroup']

def test_get_config_type():
    '''Make sure we can read the configuration type from an XML file.'''
    config = L2InputFile(test_data + "acos_Fph_090627_37_Test_v050050_Fph2900_r11_110523154640i.config")
    assert config.get_section("input->PGENameGroup")[0].get_keyword_value("PGEName") == 'GOSATFullPhysics'
    assert config["input/PGENameGroup/PGEName"] == "GOSATFullPhysics"
    print(config["input/GOSATFullPhysics"][0])
    assert config["input/GOSATFullPhysics"][0] == '20090627210055'
    assert config.sounding_ids()[0] == "20090627210055"

def test_ascii_file():
    '''Make sure we can read an ascii file'''
    config = L2InputFile(test_data + "gosat.config")
    assert config.get_section("input")[0].get_all_section_names() == ['input', 'input->PGENameGroup', 'input->InputProductFiles', 'input->AlternativeSettings', 'input->GOSATFullPhysics', 'input->GOSATFullPhysics->LIST', 'input->GOSATFullPhysics->LIST->VALUES']
    assert config.get_section("input->PGENameGroup")[0].get_keyword_value("PGEName") == "GOSATFullPhysics"
    assert config["input/PGENameGroup/PGEName"] == "GOSATFullPhysics"
    assert config["input/GOSATFullPhysics"][0] == "20100223034944"
    assert config.sounding_ids()[0] == "20100223034944"
    
    


