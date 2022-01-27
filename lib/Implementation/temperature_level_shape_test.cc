#include <boost/lexical_cast.hpp>

#include "temperature_met.h"
#include "atmosphere_fixture.h"
#include "acos_sounding_id.h"
#include "unit_test_support.h"
#include "temperature_level_shape.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_level_shape, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> temp_base(19);
  temp_base = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;

  HdfFile shape_file(test_data_dir() + "in/shape_scaling/IGRA_EOF_20220121.h5");

  // Files have 20 levels, just use first 19
  Array<double, 2> shape_prof(19, 3);
  for(int shape_idx = 0; shape_idx < shape_prof.cols(); shape_idx++) {
    shape_prof(Range::all(), shape_idx) = shape_file.read_field<double, 1>("Temperature/EOF/shape_" + boost::lexical_cast<std::string>(shape_idx + 1))(Range(0,19));
  }

  Array<double, 1> shape_scaling(3);
  shape_scaling = 0.1, 0.2, 0.3;

  Array<bool, 1> scaling_flag(shape_scaling.rows());
  scaling_flag = true;

  TemperatureLevelShape t_shape(temp_base, shape_prof, shape_scaling, config_pressure, scaling_flag);

  StateVector sv;
  sv.add_observer(t_shape);

  for(int lev_idx = 0; lev_idx < temp_base.rows(); lev_idx++) {
    AutoDerivativeWithUnit<double> pres_lev = config_pressure->pressure_grid()(lev_idx);
    double temp_calc = t_shape.temperature(pres_lev).convert(units::K).value.value();
    double temp_expect = temp_base(lev_idx) +
      shape_scaling(0) * shape_prof(lev_idx, 0) +
      shape_scaling(1) * shape_prof(lev_idx, 1) +
      shape_scaling(2) * shape_prof(lev_idx, 2);
    BOOST_CHECK_CLOSE(temp_calc, temp_expect, 1e-8);
  }
}

BOOST_AUTO_TEST_SUITE_END()

