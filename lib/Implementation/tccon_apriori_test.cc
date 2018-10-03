#include "tccon_apriori.h"
#include "level_1b_acos.h"
#include "unit_test_support.h"
#include "acos_met_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(tccon_apriori, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> sfile
    (new HdfFile(test_data_dir() + "in/sounding_id.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*sfile, "20091009203401", AcosSoundingId::S_SOUNDING));
  boost::shared_ptr<AcosMetFile> ecmwf(new AcosMetFile(test_data_dir() + "in/ecmwf.h5", 
						   sid, true));
  boost::shared_ptr<Level1bAcos> l1b(new Level1bAcos(sfile, sid));

  Array<double, 1> pecmwf = ecmwf->pressure_levels();
  Array<double, 1> tecmwf = ecmwf->temperature();

  // Test two interfaces for using class
  TcconApriori a1(ecmwf, l1b);
  TcconApriori a2(l1b, pecmwf, tecmwf, ecmwf->surface_pressure());

  std::cerr.precision(12);

  Array<double, 1> expected(91);
  expected = 0.000378671522655, 0.000378711870109, 0.00037874809863, 
    0.000378792438354, 0.000378845386017, 0.00037890728977, 
    0.000378978354167, 0.000379058650405, 0.000379148129854, 
    0.000379246639393, 0.000379353936988, 0.000379469707467, 
    0.000379593576961, 0.000379725126337, 0.000379863903451, 
    0.000380009433196, 0.000380161227178, 0.000380318790869, 
    0.00038048163038, 0.000380649257613, 0.000380821194298, 
    0.000380996976093, 0.000381176154971, 0.000381358300781, 
    0.000381543001973, 0.000381729866813, 0.000381918524902, 
    0.000382108626389, 0.000382299841854, 0.000382491857979, 
    0.000382684376347, 0.000382877087539, 0.000383069657247, 
    0.000383261721923, 0.000383452996431, 0.000383643717078, 
    0.000383834586047, 0.000384026409253, 0.000384219931221, 
    0.000384415710805, 0.000384614055847, 0.000384816221345, 
    0.000385025838296, 0.000385244728044, 0.000385471209863, 
    0.000385704674926, 0.000385945251792, 0.000386186664605, 
    0.000386393354745, 0.000386520533465, 0.000386754292694, 
    0.000387281426903, 0.000387947550126, 0.000388716998469, 
    0.000389370709282, 0.000387542685961, 0.000387291513479, 
    0.000387005191999, 0.000386679528566, 0.000386310224352, 
    0.000385893214589, 0.000385425255632, 0.000384906891049, 
    0.000384343638268, 0.000383744360869, 0.000383121291565, 
    0.000382489342201, 0.000381863809791, 0.000381259163382, 
    0.000380689370405, 0.00038016732226, 0.000379704587461, 
    0.00037931030603, 0.000378990128847, 0.000378746362752, 
    0.000378578297716, 0.000378457365817, 0.000378805616047, 
    0.000379626521834, 0.000380808635624, 0.000382226004244, 
    0.000383760919092, 0.000385316181574, 0.000386815791525, 
    0.000388202019272, 0.000389437717114, 0.000390531040575, 
    0.000391485630329, 0.000392272374519, 0.000392879866269, 
    0.000393322562575;

  for(int i = 0; i < pecmwf.rows(); ++i) {
    BOOST_CHECK_CLOSE(a1.co2_vmr(pecmwf(i)), expected(i), 1e-6);
    BOOST_CHECK_CLOSE(a2.co2_vmr(pecmwf(i)), expected(i), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(error_case)
{
  // Don't normally run this. This was for looking at a specific error
  // case that has been fixed. Leave this in place in case we need to
  // revisit this, but eventually we can just delete this code.
  return;
  // Specific error case that depends on real data. Won't run this normally
  boost::shared_ptr<HdfFile> sfile
    (new HdfFile("/acos/product/Production/v050050/L1b20900/r01/091128/acos_L1b_091128_23_Production_v050050_L1b20900_r01_110910225735.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*sfile, "20091128134921", AcosSoundingId::S_SOUNDING));
  boost::shared_ptr<AcosMetFile> ecmwf(new AcosMetFile("/acos/product/Production/v050050/Ecm20900/r01/091128/acos_Ecm_091128_23_Production_v050050_Ecm20900_r01_110910235754.h5", 
						   sid, true));
  boost::shared_ptr<Level1bAcos> l1b(new Level1bAcos(sfile, sid));
  TcconApriori a(ecmwf, l1b);

  Array<double, 1> pecmwf = ecmwf->pressure_levels();
  Array<double, 1> tecmwf = ecmwf->temperature();

  for(int i = 0; i < pecmwf.rows(); ++i)
    std::cerr << pecmwf(i) << "  " << a.co2_vmr(pecmwf(i)) << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
