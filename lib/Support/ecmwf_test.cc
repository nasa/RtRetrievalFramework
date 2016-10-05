#include "unit_test_support.h"
#include "acos_sounding_id.h"
#include "oco_sounding_id.h"
#include "acos_ecmwf.h"
#include "oco_ecmwf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ecmwf_acos, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "20091009203401";
  HdfFile sfile(test_data_dir() + "in/sounding_id.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(sfile, sid);
  AcosEcmwf e(test_data_dir() + "in/ecmwf.h5", sidv[0], sidv.size() > 1);
  BOOST_CHECK_CLOSE(e.surface_pressure(), 99682.828125, 1e-6);
  Array<double, 1> press(20);
  press = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    100000, 105000;
  Array<double, 1> texpect(20);
  texpect =
    243.90097378304642461, 214.55495206293338128, 218.03548301856420721,
    222.54439025512246531, 218.33300494385892421, 221.37919630741211563,
    227.3962815102570687,
    233.51051500723281151, 239.38876092919977623, 244.52113267396558172,
    248.71471551251292453, 251.98636348202509794, 254.54141229243495559,
    256.65775880671048981, 
    258.52334065260964735, 260.15648388783131395, 261.74845655310156189,
    261.7317782935307946, 256.31781429833779384, 257.36699412179075352;
  BOOST_CHECK_MATRIX_CLOSE(e.temperature(press), texpect);
  Array<double, 1> vmr_expect(20);
  vmr_expect = 
    6.52915085338131434138e-06, 4.74802211086569229482e-06, 4.92840307027326174783e-06,
    7.14068695346565009954e-06, 3.06465709946443932554e-05, 7.57551402530849389647e-05,
    1.16493281320304099065e-04, 1.35521593385681674623e-04, 1.83875983083302410115e-04,
    2.64376164629756744452e-04, 3.68792169005680763125e-04, 4.65498603629904477916e-04,
    5.48674415316505658877e-04, 5.73115442420909929414e-04, 6.15337784112228577613e-04,
    7.40936275682209898041e-04, 9.08908284989747017324e-04, 1.30631957917109933071e-03,
    1.41811951670399731713e-03, 1.55141541137329805507e-03;
  BOOST_CHECK_MATRIX_CLOSE(e.h2o_vmr(press), vmr_expect);

  ArrayAd<double, 1> press2(20, 20);
  press2.value() = press;
  press2.jacobian() = 0;
  for(int i = 0; i < press2.rows(); ++i)
    press2.jacobian()(i,i) = 1;
  BOOST_CHECK_MATRIX_CLOSE(e.temperature(press2).value(), texpect);
  BOOST_CHECK_MATRIX_CLOSE(e.h2o_vmr(press2).value(), vmr_expect);
  Array<double, 2> jac = e.temperature(press2).jacobian();
  Array<double, 1> t0 = e.temperature(press2).value();
  for(int i = 0; i < 20; ++i) {
    double eps = 0.1;
    press(i) += eps;
    Array<double, 1> t1 = e.temperature(press);
    // Changing the first pressure gives a larger temperature
    // differece. The values were inspected, and looked fine
    if(i == 0)
      BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-3);
    else
      BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-7);
    press(i) -= eps;
  }

  jac = e.h2o_vmr(press2).jacobian();
  t0 = e.h2o_vmr(press2).value();
  for(int i = 0; i < 20; ++i) {
    double eps = 0.1;
    press(i) += eps;
    Array<double, 1> t1 = e.h2o_vmr(press);
    BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-12);
    press(i) -= eps;
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ecmwf_oco, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "2010090900133574";
  HdfFile sfile(test_data_dir() + "/oco2_ECMWFND_80008a_111018214952d_spliced.h5");
  boost::shared_ptr<HdfSoundingId> sido(new OcoSoundingId(sfile, sid));
  OcoEcmwf e(test_data_dir() + "/oco2_ECMWFND_80008a_111018214952d_spliced.h5", sido);
  BOOST_CHECK_CLOSE(e.surface_pressure(), 95796, 1e-6);
  BOOST_CHECK_CLOSE(e.windspeed(), 9.5288243675028176938, 1e-6);
  Array<double, 1> press(20);
  press = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    100000, 105000;
  Array<double, 1> texpect(20);
  texpect =
    249.99932014712481987, 204.71691671141664415, 205.14682181545998674,
    206.67900774692736832, 211.72013798983311972, 218.03120044161263991,
    222.97529589441276698, 227.78345303636547214, 233.27363833074147692,
    238.46076429718232248, 242.14518070680597361, 245.05993612140508731,
    248.07973362637457626, 251.15477955996385617, 253.2607227236587164,
    253.89281922755003507, 249.47060311981078939, 253.32026131834302873,
    257.52964972931016518, 261.61106752656741037;
  BOOST_CHECK_MATRIX_CLOSE(e.temperature(press), texpect);
  Array<double, 1> vmr_expect(20);
  vmr_expect = 
    6.1846336386044464047e-06, 2.3696526616669067966e-06, 1.7912447195779596004e-06,
    4.4493971613491540648e-06, 4.2663811750629721735e-06, 2.310533211457202921e-05,
    4.3066677729604225546e-05, 8.0495354151999834675e-05, 0.00012825322759659445052,
    0.00011939929116084474232, 0.00013114229383684339217, 0.00019820301866710661053,
    0.00050859154716482578729, 0.00060071598235635829781, 0.00049173906754511310067,
    0.00043043330031180778751, 0.00078522606393224120801, 0.001024082783212755662,
    0.0014987653090951654714, 0.0021981331635712191876;
  BOOST_CHECK_MATRIX_CLOSE(e.h2o_vmr(press), vmr_expect);

  ArrayAd<double, 1> press2(20, 20);
  press2.value() = press;
  press2.jacobian() = 0;
  for(int i = 0; i < press2.rows(); ++i)
    press2.jacobian()(i,i) = 1;
  BOOST_CHECK_MATRIX_CLOSE(e.temperature(press2).value(), texpect);
  BOOST_CHECK_MATRIX_CLOSE(e.h2o_vmr(press2).value(), vmr_expect);
  Array<double, 2> jac = e.temperature(press2).jacobian();
  Array<double, 1> t0 = e.temperature(press2).value();
  for(int i = 0; i < 20; ++i) {
    double eps = 0.1;
    press(i) += eps;
    Array<double, 1> t1 = e.temperature(press);
    // Changing the first pressure gives a larger temperature
    // differece. The values were inspected, and looked fine
    if(i == 0)
      BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-3);
    else
      BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-7);
    press(i) -= eps;
  }

  jac = e.h2o_vmr(press2).jacobian();
  t0 = e.h2o_vmr(press2).value();
  for(int i = 0; i < 20; ++i) {
    double eps = 0.1;
    press(i) += eps;
    Array<double, 1> t1 = e.h2o_vmr(press);
    BOOST_CHECK_MATRIX_CLOSE_TOL((t1 - t0) / eps, jac(Range::all(), i), 1e-10);
    press(i) -= eps;
  }

}

BOOST_AUTO_TEST_SUITE_END()
