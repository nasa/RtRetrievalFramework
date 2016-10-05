#include "unit_test_support.h"
#include "output_hdf_iteration.h"
#include <cerrno>

using namespace FullPhysics;
using namespace blitz;

class CalcData {
public:
  int d_val() const {return d;}
  Array<int, 2> a_val() const {return a;}
  int d;
  Array<int, 2> a;
};

BOOST_FIXTURE_TEST_SUITE(output_hdf_iteration, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFileGenerating> 
    hf(new HdfFileGenerating("output_hdf_iteration.h5"));
  add_file_to_cleanup("output_hdf_iteration.h5");
  OutputHdfIteration h(hf);
  boost::shared_ptr<CalcData> c1(new CalcData);
  h.register_data_source("/Test/d_val", &CalcData::d_val, c1);
  h.register_data_source("/Test/a_val", &CalcData::a_val, c1);
  c1->d = 1;
  c1->a.resize(2, 3);
  c1->a = 
    1,2,3,
    4,5,6;
  h.write();
  c1->d = 2;
  c1->a =
    7,8,9,
    10,11,12;
  h.write();
  h.close();
  hf->close();
  HdfFile hread("output_hdf_iteration.h5");
  Array<int, 2> d_expect(1, 2);
  d_expect = 1, 2;
  Array<int, 4> a_expect(1, 2, 2, 3);
  a_expect = 
    1,2,3,
    4,5,6,

    7,8,9,
    10,11,12;
  Array<int, 2> d_read(hread.read_field<int, 2>("/Iteration/Test/d_val"));
  Array<int, 4> a_read(hread.read_field<int, 4>("/Iteration/Test/a_val"));
  BOOST_CHECK_MATRIX_CLOSE(d_read, d_expect);
  BOOST_CHECK_MATRIX_CLOSE(a_read, a_expect);
}

BOOST_AUTO_TEST_SUITE_END()
