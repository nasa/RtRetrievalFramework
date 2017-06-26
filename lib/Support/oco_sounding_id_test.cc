#include "unit_test_support.h"
#include "oco_sounding_id.h"

/// Place holder for better sounding id reading code
// TESTING
namespace FullPhysics {
class TestMe {
public:

  template<class T, int D> blitz::Array<T, D> read_sounding_field(const HdfFile& file, const std::string& dataset_name) const
  {
    std::string shape_name = file.read_attribute<std::string>(dataset_name + "/Shape");
    std::vector<std::string> dim_names = file.read_attribute<std::vector<std::string> >("/Shapes/" + shape_name + "/Dimensions");
    
    return read_sounding_field<T, D>(file, dim_names);
  }
  template<class T, int D> blitz::Array<T, D> read_sounding_field(const HdfFile& file, std::vector<std::string>& dim_names)  const
  {
    std::string frame_dim = "Exposure";
    std::string sounding_dim = "Polarization";
    
    blitz::Array<T, D> out_data;
    return out_data;
  }
 
};
}
// TESTING

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oco_sounding_id, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "/oco2_L1bScND_80008a_111017225030d_spliced.h5");

  OcoSoundingId sid_1(h, "2010090900133575");
  BOOST_CHECK_EQUAL(sid_1.frame_number(), 0);
  BOOST_CHECK_EQUAL(sid_1.sounding_id(), (int64_t) 2010090900133575);
  BOOST_CHECK_EQUAL(sid_1.sounding_number(), 2);

  OcoSoundingId sid_2(h, "2010090900134374");
  BOOST_CHECK_EQUAL(sid_2.frame_number(), 6);
  BOOST_CHECK_EQUAL(sid_2.sounding_id(), (int64_t) 2010090900134374);
  BOOST_CHECK_EQUAL(sid_2.sounding_number(), 1);

}

BOOST_AUTO_TEST_CASE(read_sounding_field)
{
  HdfFile h(test_data_dir() + "l1b.h5");
  TestMe x;
  x.read_sounding_field<double, 3>(h, "FootprintGeometry/footprint_altitude");
}

BOOST_AUTO_TEST_SUITE_END()
