#include <boost/lexical_cast.hpp>

#include "atmosphere_fixture.h"
#include "acos_sounding_id.h"
#include "unit_test_support.h"
#include "absorber_vmr_shape.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absorber_vmr_shape, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> vmr_base(19);
  vmr_base = 6.52839575e-06, 4.74787247e-06, 4.92904187e-06, 7.15214080e-06, 3.08558918e-05,
      7.61965519e-05, 1.16499294e-04, 1.36102208e-04, 1.84245383e-04, 2.64503160e-04, 3.69431451e-04,
      4.66057518e-04, 5.48717085e-04, 5.73098800e-04, 6.15544205e-04, 7.41377218e-04,
      9.09569234e-04, 1.30658765e-03, 1.46738835e-03;

  HdfFile shape_file(test_data_dir() + "in/shape_scaling/IGRA_EOF_20220121.h5");

  // Files have 20 levels, just use first 19
  Array<double, 2> shape_prof(19, 3);
  for(int shape_idx = 0; shape_idx < shape_prof.cols(); shape_idx++) {
    shape_prof(Range::all(), shape_idx) = shape_file.read_field<double, 1>("Gas/H2O/EOF/shape_" + boost::lexical_cast<std::string>(shape_idx + 1))(Range(0,18));
  }

  Array<double, 1> shape_scaling(3);
  shape_scaling = 0.1, 0.2, 0.3;

  Array<bool, 1> scaling_flag(shape_scaling.rows());
  scaling_flag = true;

  AbsorberVmrShape absorber_shape(vmr_base, shape_prof, shape_scaling, config_pressure, scaling_flag, "H2O");

  StateVector sv;
  sv.add_observer(absorber_shape);

  for(int lev_idx = 0; lev_idx < vmr_base.rows(); lev_idx++) {
    AutoDerivative<double> pres_lev = config_pressure->pressure_grid()(lev_idx).value;
    double vmr_calc = absorber_shape.volume_mixing_ratio(pres_lev).value();
    double vmr_expect = vmr_base(lev_idx) +
      shape_scaling(0) * shape_prof(lev_idx, 0) +
      shape_scaling(1) * shape_prof(lev_idx, 1) +
      shape_scaling(2) * shape_prof(lev_idx, 2);
    BOOST_CHECK_CLOSE(vmr_calc, vmr_expect, 1e-8);
  }
}

BOOST_AUTO_TEST_SUITE_END()
