#include "aerosol_optical.h"
#include "pressure.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "aerosol_property_hdf.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_optical, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(config_file)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  Array<double, 1> od_expect(4);
  od_expect =  0.0300339738374, 0.0300339738374, 
    0.0325219546595, 0.036014594462;
  for(int aer_idx = 0; aer_idx < od_expect.rows(); aer_idx++)
    BOOST_CHECK_CLOSE(a->aerosol_optical_depth(aer_idx), od_expect(aer_idx), 1e-4);
  double total_expect = sum(od_expect);
  BOOST_CHECK_CLOSE(a->aerosol_optical_depth_total(), total_expect, 1e-4);
}

BOOST_AUTO_TEST_CASE(layer_parameters)
{
  firstIndex i1; secondIndex i2;
  // Expected values were gotten by running the old Fortran code and
  // extracting out the answer from that.
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  Array<double, 1> od_expect(18);
  Array<double, 1> ssa_expect(18);
  od_expect = 2.74194178225e-11, 1.57679009526e-05,
    0.0046653557663, 0.00940062333069, 0.00904058828474,
    0.00542317357956, 0.00514668378488, 0.00563371951229,
    0.00658779554283, 0.00765016197018, 0.00857484138191,
    0.00926322684202, 0.00973244118942, 0.0100665759631,
    0.0103748368244, 0.0107641227421, 0.0113251652964, 0.00416396423;
  ssa_expect = 2.62411725803e-11, 1.57624098856e-05, 0.00466213741271,
    0.00938712429319, 0.00899661403373, 0.0053483240547, 0.00501076457573,
    0.00542243887556, 0.00629774362343, 0.00728802039782, 0.008153167367,
    0.00879544773679, 0.0092289066228, 0.0095323323439,
    0.00980895573824, 0.0101600695552, 0.010671683777, 0.00391749916961;
  BOOST_CHECK_MATRIX_CLOSE(sum(a->optical_depth_each_layer(12929.94).value(),i2),
			   od_expect);
  BOOST_CHECK_MATRIX_CLOSE(a->ssa_each_layer(12929.94).value(), ssa_expect);
  od_expect = 2.74196241536e-11, 1.57678789317e-05, 0.00466534958506,
    0.00940061234245, 0.00904058419521, 0.00542318125538, 
    0.00514670462721, 0.00563375538688, 0.0065878464798, 
    0.00765022642059, 0.00857491693238, 0.00926331104332, 
    0.00973253221078, 0.0100666729585, 0.0103749400361,
    0.0107642334269, 0.0113252855649, 0.00416400976705;
  ssa_expect = 2.62412156853e-11, 1.57623889568e-05, 0.00466213132914, 
    0.0093871125113, 0.00899660480613, 0.00534832180726, 0.00501076684348, 
    0.00542244559708, 0.00629775440118, 0.00728803463588, 0.00815318440425, 
    0.00879546698692, 0.00922892768541, 0.00953235506478, 
    0.00980898022095, 0.0101600961393, 0.010671713, 0.00391751034701;
  BOOST_CHECK_MATRIX_CLOSE(sum(a->optical_depth_each_layer(12930.30).value(),i2),
			   od_expect);
  BOOST_CHECK_MATRIX_CLOSE(a->ssa_each_layer(12930.30).value(), ssa_expect);
  
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(aerosol_jac, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 2> od = a->optical_depth_each_layer(wn);
  Array<double, 2> od0(od.shape());
  od0 = od.value();
  Array<double, 3> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(od0.shape());
    jacfd = (a->optical_depth_each_layer(wn).value() - od0) / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), Range::all(), i), jacfd, 
				 5e-9);
  }
}

BOOST_AUTO_TEST_CASE(ssa_jac)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  double wn = 4820.0;
  ArrayAd<double, 1> ssa = a->ssa_each_layer(wn);
  Array<double, 1> ssa0(ssa.shape());
  ssa0 = ssa.value();
  Array<double, 2> jac = ssa.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(ssa0.shape());
    jacfd = (a->ssa_each_layer(wn).value() - ssa0) / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 
				 5e-9);
  }
}

BOOST_AUTO_TEST_SUITE_END()


