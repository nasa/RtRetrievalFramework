#include "unit_test_support.h"
#include "output_hdf.h"
#include <cerrno>

using namespace FullPhysics;
using namespace blitz;

class CalcData {
public:
  CalcData(int D) : d(D) {}
  int d_val() const {return d;}
  int d_val_fake_error() const {throw Exception("fake error");}
private:
  int d;
};

BOOST_FIXTURE_TEST_SUITE(output_hdf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFileGenerating> 
    hf(new HdfFileGenerating("output_hdf.h5"));
  add_file_to_cleanup("output_hdf.h5");
  OutputHdf h(hf, 20, 112, 5, 3);
  boost::shared_ptr<CalcData> c1(new CalcData(1));
  boost::shared_ptr<CalcData> c2(new CalcData(2));
  h.register_data_source("/Test/C1", &CalcData::d_val, c1);
  h.register_data_source("/Test/C2", &CalcData::d_val, c2);
  // Try each of the hardcoded shapes
  Array<int, 1> iarr;
  Array<double, 1> darr;
  iarr.resize(20);
  darr.resize(20);
  iarr = 1;
  darr = 2.0;
  h.register_data_source("/Test/ilev", iarr);
  h.register_data_source("/Test/dlev", darr);
  iarr.resize(112);
  darr.resize(112);
  iarr = 1;
  darr = 2.0;
  h.register_data_source("/Test/istate", iarr);
  h.register_data_source("/Test/dstate", darr);
  iarr.resize(5);
  darr.resize(5);
  iarr = 1;
  darr = 2.0;
  h.register_data_source("/Test/iaer", iarr);
  h.register_data_source("/Test/daer", darr);
  // And one unrecognized size
  iarr.resize(12);
  darr.resize(12);
  iarr = 1;
  darr = 2.0;
  h.register_data_source("/Test/iunr", iarr);
  h.register_data_source("/Test/dunr", darr);

  // Same for 2d arrays.
  Array<int, 2> iarr2(112, 112);
  Array<double, 2> darr2(112, 112);
  iarr2 = 1;
  darr2 = 2.0;
  h.register_data_source("/Test/istatestate", iarr2);
  h.register_data_source("/Test/dstatestate", darr2);
  iarr2.resize(10, 10);
  darr2.resize(10, 10);
  iarr2 = 1;
  darr2 = 2.0;
  h.register_data_source("/Test/iunkunk", iarr2);
  h.register_data_source("/Test/dunkunk", darr2);
  h.write();
  hf->close();
  HdfFile hread("output_hdf.h5");
  BOOST_CHECK_EQUAL(hread.read_field<int>("/Test/C1"), 1);
  BOOST_CHECK_EQUAL(hread.read_field<int>("/Test/C2"), 2);
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/C1/Shape"),
		    "Retrieval_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/daer/Shape"),
		    "Retrieval_Aerosol_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/iaer/Shape"),
		    "Retrieval_Aerosol_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/dlev/Shape"),
		    "Retrieval_Level_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/ilev/Shape"),
		    "Retrieval_Level_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/dstate/Shape"),
		    "Retrieval_StateVectorElement_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/istate/Shape"),
		    "Retrieval_StateVectorElement_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/dstatestate/Shape"),
		    "Retrieval_StateVectorElement_StateVectorElement_Array");
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>("/Test/istatestate/Shape"),
		    "Retrieval_StateVectorElement_StateVectorElement_Array");
}

BOOST_AUTO_TEST_CASE(test_error)
{
  boost::shared_ptr<HdfFileGenerating> 
    hf1(new HdfFileGenerating("output_hdf.h5"));
  boost::shared_ptr<HdfFileGenerating> 
    hf2(new HdfFileGenerating("output_hdf.h5.error"));
  OutputHdf h(hf1, 20, 112, 5, 3);
  OutputHdf herr(hf2, 20, 112, 5, 3);
  add_file_to_cleanup("output_hdf.h5");
  add_file_to_cleanup("output_hdf.h5.error");
  boost::shared_ptr<CalcData> c1(new CalcData(1));
  boost::shared_ptr<CalcData> c2(new CalcData(2));
  h.register_data_source("/Test/C1", &CalcData::d_val_fake_error, c1);
  h.register_data_source("/Test/C2", &CalcData::d_val, c2);
  herr.register_data_source("/Test/C1", &CalcData::d_val_fake_error, c1);
  herr.register_data_source("/Test/C2", &CalcData::d_val, c2);
  BOOST_CHECK_THROW(h.write(), Exception);
  herr.write_best_attempt();
  hf1->close();
  hf2->close();
  // Normal file should not exist
  BOOST_CHECK(fopen("output_hdf.h5", "r") ==0);
  BOOST_CHECK_EQUAL(errno, ENOENT);
  // Error file should have C2, but not C1
  HdfFile hread("output_hdf.h5.error");
  BOOST_CHECK_EQUAL(hread.read_field<int>("/Test/C2"), 2);
  BOOST_CHECK_THROW(hread.read_field<int>("/Test/C1"), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
