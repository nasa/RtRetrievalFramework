#include "heritage_file.h"
#include "unit_test_support.h"
#include "fp_exception.h"
using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(heritage_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  HeritageFile c(test_data_dir() + "heritage_file_test.run");
  c.parse_file(test_data_dir() + "in/l1b/soundinginfo.dat");
  BOOST_CHECK(c.has_value("CONTROL/run_mode"));
  BOOST_CHECK(c.has_value("control/RUN_MODE"));
  BOOST_CHECK_EQUAL(c.file_value("SOUNDING_INFO/spectrum_file"),
  		    test_data_dir() + "in/l1b/spec/spectra.dat");
  BOOST_CHECK_EQUAL(c.directory_base(), test_data_dir());
  BOOST_CHECK(!c.has_value("Badkeyword"));
  BOOST_CHECK_EQUAL(c.value<int>("ALGORITHMS/points_sun"), 10000);
  BOOST_CHECK_EQUAL(c.value<int>("algorithms/points_sun"), 10000);
  BOOST_CHECK(c.value<bool>("ALGORITHMS/h2o_correction"));
  BOOST_CHECK(c.value<bool>("algorithms/h2o_correction"));
  BOOST_CHECK(!c.value<bool>("ALGORITHMS/inst_doppler_shift"));
  Time texpect = Time::parse_time("2006-09-14T12:27:22.001Z");
  BOOST_CHECK_CLOSE(c.value<Time>("frame_time_stamp").unix_time(), 
		    texpect.unix_time(), 1e-6);
  BOOST_CHECK(c.has_value("CONTROL/fake_empty_for_test"));
  BOOST_CHECK_EQUAL(c.value<std::string>("CONTROL/fake_empty_for_test"), "");
  std::vector<int> vi = 
    c.value<std::vector<int> >("CONTROL/fake_empty_for_test");
  BOOST_CHECK_EQUAL((int) vi.size(), 0);
  std::vector<double> vd = 
    c.value<std::vector<double> >("CONTROL/fake_empty_for_test");
  BOOST_CHECK_EQUAL((double) vd.size(), 0);
  std::vector<std::string> vs = 
    c.value<std::vector<std::string> >("CONTROL/fake_empty_for_test");
  BOOST_CHECK_EQUAL((int) vs.size(), 0);
}

BOOST_AUTO_TEST_CASE(intvec_test)
{
  HeritageFile c(test_data_dir() + "intvec_test.txt");
  std::vector<int> invec1_expect, invec2_expect;
  invec1_expect.push_back(1);
  invec1_expect.push_back(2);
  invec1_expect.push_back(3);
  invec1_expect.push_back(4);
  invec1_expect.push_back(5);
  invec2_expect.push_back(4);
  invec2_expect.push_back(5);
  invec2_expect.push_back(6);
  invec2_expect.push_back(7);
  invec2_expect.push_back(8);
  invec2_expect.push_back(15);
  invec2_expect.push_back(40);
  invec2_expect.push_back(41);
  invec2_expect.push_back(42);
  invec2_expect.push_back(43);
  invec2_expect.push_back(44);
  invec2_expect.push_back(45);
  BOOST_CHECK(c.value<std::vector<int> >("invec1") == invec1_expect);
  BOOST_CHECK(c.value<std::vector<int> >("invec2") == invec2_expect);
}

BOOST_AUTO_TEST_CASE(multiple_section_test)
{
  HeritageFile c(test_data_dir() + "heritage_file_test.run");
  BOOST_CHECK_EQUAL(c.number_block("SOUNDING_INFO"), 1);
  BOOST_CHECK_EQUAL(c.number_block("FAKESEC"), 0);
  BOOST_CHECK_EQUAL(c.number_block("PARAMETER_DEFINITION/GAS"), 3);
  BOOST_CHECK_EQUAL(c.value<std::string>("PARAMETER_DEFINITION/GAS/name", 0), 
		    "CO2");
  BOOST_CHECK_EQUAL(c.value<std::string>("PARAMETER_DEFINITION/GAS/name", 1), 
		    "H2O");
  BOOST_CHECK_EQUAL(c.value<std::string>("PARAMETER_DEFINITION/GAS/name", 2), 
		    "O2");
  BOOST_CHECK_EQUAL(c.value<std::string>("PARAMETER_DEFINITION/GAS/O2/name"), 
		    "O2");
}

BOOST_AUTO_TEST_CASE(bad_value)
{
  HeritageFile c(test_data_dir() + "heritage_file_test.run");
  CHECK_THROW_EXCEPTION(c.value<bool>("WINDOW_INFO/spectral_window_file"),
"Error reading file. The value for WINDOW_INFO/spectral_window_file is not\n"
"'true' or 'false'\n"
"File: " + test_data_dir() + "heritage_file_test.run\n"
"Line: 76"
			);
  CHECK_THROW_EXCEPTION(c.value<Time>("WINDOW_INFO/spectral_window_file"),
"Error reading file. The value for WINDOW_INFO/spectral_window_file is not\n"
"a datetime\n"
"File: " + test_data_dir() + "heritage_file_test.run\n"
"Line: 76"
			);
  CHECK_THROW_EXCEPTION(c.value<int>("WINDOW_INFO/spectral_window_file"),
"Error reading file. The value for WINDOW_INFO/spectral_window_file is not\n"
"the expected type\n"
"File: " + test_data_dir() + "heritage_file_test.run\n"
"Line: 76"
			);
  BOOST_CHECK_EQUAL(
    c.value<std::string>("WINDOW_INFO/spectral_window_file"),
    "oco_l2.win");
  BOOST_CHECK_EQUAL(
    c.value<std::string>("window_info/spectral_window_file"),
    "oco_l2.win");
  CHECK_THROW_EXCEPTION(HeritageFile c(test_data_dir() + "bad_no_header.dat"), 
"Error reading file. End of header seen before header started\n"
"File: " + test_data_dir() + "bad_no_header.dat\n"
"Line: 1"
			);
  CHECK_THROW_EXCEPTION(HeritageFile c(test_data_dir() + "bad_short_matrix_1.dat"), 
"Error reading file. Not enough lines of data\n"
"File: " + test_data_dir() + "bad_short_matrix_1.dat\n"
"Line: 17"
 			);
  CHECK_THROW_EXCEPTION(HeritageFile c(test_data_dir() + "bad_short_matrix_2.dat"), 
"Error reading file. Trouble processing a row of data\n"
"File: " + test_data_dir() + "bad_short_matrix_2.dat\n"
"Line: 13"
			);
}

BOOST_AUTO_TEST_CASE(value_vector)
{
  HeritageFile c(test_data_dir() + "in/l1b/soundinginfo.dat");
  std::vector<double> r = c.value<std::vector<double> >("sounding_latitude");
  BOOST_CHECK_EQUAL((int) r.size(), 3);
  BOOST_CHECK_CLOSE(r[0], 77.1828918457, 1e-8);
  BOOST_CHECK_CLOSE(r[1], 77.1828918457, 1e-8);
  BOOST_CHECK_CLOSE(r[2], 77.1828918457, 1e-8);
}

BOOST_AUTO_TEST_CASE(read_matrix_file)
{
  HeritageFile c(test_data_dir() + "old_ascii/solar_cont_v1.dat");
  BOOST_CHECK_EQUAL(c.value<int>("solar_model_version"), 2);
  BOOST_CHECK_EQUAL(c.value<std::string>("File_ID"), "Solar Parameters");
  std::vector<std::string> r = c.value<std::vector<std::string> >("Labels");
  BOOST_CHECK_EQUAL((int) r.size(), 2);
  BOOST_CHECK_EQUAL(r[0], "PARAMETER");
  BOOST_CHECK_EQUAL(r[1], "SOLAR");
  blitz::Array<double, 2> dat_expected(6, 2);
  dat_expected =  
    1, 8.83596E21,
    2, -9.48206E20,
    3, -1.517E22,  
    4, 1.74114E22, 
    5, -7.73485E21,
    6, 1.2313E21;
  BOOST_CHECK(max(abs(c.data() - dat_expected)) < 1e-8);
}

BOOST_AUTO_TEST_CASE(read_column)
{
  HeritageFile c(test_data_dir() + "old_ascii/aerosol_015_log_20.dat");
  BOOST_CHECK_THROW(c.data("bad_col"), Exception);
  blitz::Array<double, 1> dexpect(20);
  dexpect = -34.538776, -34.538776, -19.288269, -14.815826, -14.391203, 
    -14.675235, -14.947882, -15.110210, -15.116857, -15.031474, -14.926888,
    -14.836818, -14.766849, -14.710517, -14.657920, -14.599772, -14.529411,
    -14.443572, -14.342169, -14.227380;
  BOOST_CHECK_MATRIX_CLOSE(c.data("water"), dexpect);
  dexpect= -34.538776, -34.538776, -18.882804, -14.413419, -13.992814,
    -14.291637, -14.590424, -14.794397, -14.848530, -14.810188, -14.757100,
    -14.731176, -14.745469, -14.800121, -14.890932, -15.012790, -15.160861,
    -15.330923, -15.519402, -15.723299;
  BOOST_CHECK_MATRIX_CLOSE(c.data("ice"), dexpect);
}
BOOST_AUTO_TEST_SUITE_END()
