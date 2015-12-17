#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "dispersion_fit.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(dispersion_fit_gosat, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  DoubleWithUnit solar_line_location(12985.16325e0, units::inv_cm);
  DoubleWithUnit solar_line_width(0.2*0.2, units::inv_cm);
  DoubleWithUnit search_width(2.83, units::inv_cm);
  // this takes into account the offset of the GOSAT band 1 ILS.
  DoubleWithUnit ils_offset(0.203e0, units::inv_cm);

  ArrayWithUnit<double, 1> offset_scaling; 
  offset_scaling.value.resize(3);
  offset_scaling.value = 1, 0.48, 0.38;
  offset_scaling.units = Unit("cm^-1");

  // Check that we get a solution similar to old code
  Array<double, 2> disp_expt(3,2);
  disp_expt =
    1.28695626e+04, 1.99493000e-01,
    5.74981805e+03, 1.99493000e-01,
    4.74980179e+03, 1.99493000e-01;

  Array<double, 2> disp_initial(disp_expt.shape());
  disp_initial = 
    12869.9, 0.199493,
    5749.98, 0.199493,
    4749.93, 0.199493;
  DispersionFit disp_fit = DispersionFit(config_level_1b);
  Array<double, 2> disp_sol = disp_fit.fit(disp_initial, solar_line_location, solar_line_width, search_width, ils_offset, offset_scaling);

  BOOST_CHECK_MATRIX_CLOSE_TOL(disp_expt, disp_sol, 1e-4);

  // Check that we can perturb slighly and get back a scaling equal to
  // the offset scaling
  Array<double, 2> disp_pert(disp_sol.shape());
  disp_pert = disp_sol;
  disp_pert(Range::all(), 0) = disp_pert(Range::all(), 0) + 0.1;
  Array<double, 2> disp_sol_pert = disp_fit.fit(disp_pert, solar_line_location, solar_line_width, search_width, ils_offset, offset_scaling);
  
  Array<double, 2> disp_expt_pert_diff(disp_expt.shape());
  disp_expt_pert_diff =
    -0.1,    0.0,
    -0.048,  0.0,
    -0.038,  0.0;
  BOOST_CHECK_MATRIX_CLOSE(disp_expt_pert_diff, disp_sol_pert - disp_pert);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(dispersion_fit_oco, ConfigurationOco2Fixture)

BOOST_AUTO_TEST_CASE(oco)
{
  DoubleWithUnit solar_line_location(0.77010968653012513, units::micron);
  DoubleWithUnit solar_line_width(pow(1.6e-5,2), units::micron);
  DoubleWithUnit search_width(1.682e-4, units::micron);
  DoubleWithUnit ils_offset(0, units::micron);

  ArrayWithUnit<double, 1> offset_scaling; 
  offset_scaling.value.resize(3);
  offset_scaling.value = 1, 0.48, 0.38;
  offset_scaling.units = Unit("cm^-1");

  Array<double, 2> disp_expt(3,5);
  disp_expt = 
    0.757453, 1.7355e-05, -2.78491e-09,            0,          0, 
    1.58995, 3.62686e-05, -5.76115e-09, -4.33245e-14,          0,
    2.04122, 4.70061e-05, -8.33508e-09,  8.14472e-13, -1.797e-16;

  Array<double, 2> disp_initial(disp_expt.shape());
  disp_initial = disp_expt;
  disp_initial(Range::all(), 0) += 1e-4 * offset_scaling.value(Range::all());;

  DispersionFit disp_fit = DispersionFit(config_level_1b);
  Array<double, 2> disp_sol = disp_fit.fit(disp_initial, solar_line_location, solar_line_width, search_width, ils_offset, offset_scaling);

  BOOST_CHECK_MATRIX_CLOSE_TOL(disp_expt(Range::all(), 0), disp_sol(Range::all(), 0), 2e-6);

  Array<double, 1> diff(disp_sol(Range::all(), 0)-disp_expt(Range::all(), 0)) ;
  /*std::cout << "Dispersion solution: " << disp_sol << std::endl
    << "Diff: " << diff << std::endl;*/

}

BOOST_AUTO_TEST_SUITE_END()

