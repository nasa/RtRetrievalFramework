#include "closest_point.h"
#include "nonuniform_spectrum_sampling.h"
#include "uniform_spectrum_sampling.h"
#include "lsi_rt.h"
#include "fp_logger.h"
#include "lidort_fixture.h"
#include "unit_test_support.h"
#include "hdf_file.h"

//#include <fstream>


using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(nonuniform_lsi_rt, LidortLowHighLambertianFixture)


void refl_compare( 
  const Array<double, 1> &grid_1, const Array<double, 1> &refl_1,
  const Array<double, 1> &grid_2, const Array<double, 1> &refl_2)
{
  BOOST_REQUIRE_EQUAL(grid_1.rows(), refl_1.rows());
  BOOST_REQUIRE_EQUAL(grid_2.rows(), refl_2.rows());
  Array<int, 1> closest_indx = closest_point(grid_1, grid_2);
  Array<double, 1> temp_grid_2(closest_indx.rows());
  Array<double, 1> temp_refl_2(closest_indx.rows());
  for(int i=0; i<closest_indx.rows(); i++) {
    temp_grid_2(i) = grid_1(closest_indx(i));
    temp_refl_2(i) = refl_1(closest_indx(i));
  }
  BOOST_CHECK_EQUAL( sum(abs(grid_2-temp_grid_2)),  0 );
  double temp = max(abs((refl_2-temp_refl_2)/(refl_2+temp_refl_2)*2.0));
  double tol = 0.001;
  std::cout << std::endl << "Checking if the following true: " << temp << " < " << tol << std::endl;
  BOOST_CHECK(temp < tol);
}


void jac_compare( 
  const Array<double, 1> &grid_1, const Array<double, 2> &jac_1,
  const Array<double, 1> &grid_2, const Array<double, 2> &jac_2)
{
  BOOST_REQUIRE_EQUAL(grid_1.rows(), jac_1.rows());
  BOOST_REQUIRE_EQUAL(grid_2.rows(), jac_2.rows());
  BOOST_REQUIRE_EQUAL(jac_1.cols(), jac_2.cols());
  Array<int, 1> closest_indx = closest_point(grid_1, grid_2);
  Array<double, 1> temp_grid_2(closest_indx.rows());
  Array<double, 2> temp_jac_2(jac_2.shape());
  Range all = Range::all();
  for(int i=0; i<closest_indx.rows(); i++) {
    temp_grid_2(i) = grid_1(closest_indx(i));
    temp_jac_2(i,all) = jac_1(closest_indx(i),all);
  }
  BOOST_CHECK_EQUAL( sum(abs(grid_2-temp_grid_2)),  0 );
  double tol = 0.001;
  for(int i_c=0; i_c<jac_2.cols(); i_c++) {
    double col_diff_norm_L1 = sum(abs(temp_jac_2(all,i_c) - jac_2(all,i_c)));
    double col_norm_L1 = sum(abs(temp_jac_2(all,i_c) + jac_2(all,i_c))) / 2.0;
    if(col_norm_L1 != 0) {
      double temp = col_diff_norm_L1/col_norm_L1;
      std::cout << std::endl << "Checking if the following true for column indexed " << i_c << ": " << temp << " < " << tol << std::endl;
      BOOST_CHECK(temp < 0.001);
    }
  }
}


BOOST_AUTO_TEST_CASE(stokes)
{
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();	     // Have log output show up.


  HeritageFile c1(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_abo2_oco__absco_v3.1.0__wn.dat");
  HeritageFile c2(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_wco2_oco__absco_v3.1.0__wn.dat");
  HeritageFile c3(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_sco2_oco__absco_v3.1.0__wn.dat");
  boost::shared_ptr<UniformSpectrumSampling> 
    unif_rt_grid(new UniformSpectrumSampling(12928.0-1, 13213.0+1, 0.01,
					     6145.0-1,  6307.0+1, 0.01,
					     4790.0-1,  4917.0+1, 0.01));
  NonuniformSpectrumSampling nonunif_rt_grid(c1,c2,c3, unif_rt_grid);

  HdfFile config(test_data_dir() + "l2_fixed_level_static_input.h5");
  LsiRt rt(low_rt, high_rt, config);

  std::cout << nonunif_rt_grid;
  std::cout << unif_rt_grid;

  for(int i = 0 ; i < 3; ++i) {
    SpectralDomain sd1 = unif_rt_grid->spectral_domain(i, lowres_grid(i), 
						       ils_half_width(i));
    SpectralDomain sd2 = nonunif_rt_grid.spectral_domain(i, lowres_grid(i), 
							 ils_half_width(i));
    Array<double, 1> refl_unif_1 = 
      rt.reflectance(sd1, i, true).spectral_range().data();
    Array<double, 1> refl_nonunif_1 = 
      rt.reflectance(sd2, 0, true).spectral_range().data();
    refl_compare(sd1.wavenumber(), refl_unif_1,
		 sd2.wavenumber(), refl_nonunif_1);
  }
}

BOOST_AUTO_TEST_CASE(stokes_and_jacobian)
{
  is_really_long_test();     // Skip unless we are running long tests.
  turn_on_logger();	     // Have log output show up.



  HeritageFile c1(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_abo2_oco__absco_v3.1.0__wn.dat");
  HeritageFile c2(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_wco2_oco__absco_v3.1.0__wn.dat");
  HeritageFile c3(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_sco2_oco__absco_v3.1.0__wn.dat");
  boost::shared_ptr<UniformSpectrumSampling> 
    unif_rt_grid(new UniformSpectrumSampling(12928.0-1, 13213.0+1, 0.01,
					     6145.0-1,  6307.0+1, 0.01,
					     4790.0-1,  4917.0+1, 0.01));
  NonuniformSpectrumSampling nonunif_rt_grid(c1,c2,c3, unif_rt_grid);

  HdfFile config(test_data_dir() + "l2_fixed_level_static_input.h5");
  LsiRt rt(low_rt, high_rt, config);

  std::cout << nonunif_rt_grid;
  std::cout << unif_rt_grid;

// std::ofstream outf;
// outf.precision( 15 );
// std::string sid;

  for(int i = 0 ; i < 3; ++i) {
    SpectralDomain sd1 = unif_rt_grid->spectral_domain(i, lowres_grid(i), 
						       ils_half_width(i));
    SpectralDomain sd2 = nonunif_rt_grid.spectral_domain(i, lowres_grid(i), 
							 ils_half_width(i));
    ArrayAd<double, 1> refl_unif_1 = 
      rt.reflectance(sd1, i).spectral_range().data_ad();
    ArrayAd<double, 1> refl_nonunif_1 = 
      rt.reflectance(sd2, i).spectral_range().data_ad();
    refl_compare(sd1.wavenumber(), refl_unif_1.value(),
		 sd2.wavenumber(), refl_nonunif_1.value());
    jac_compare(sd1.wavenumber(), refl_unif_1.jacobian(),
		sd2.wavenumber(), refl_nonunif_1.jacobian());
  }
}



BOOST_AUTO_TEST_SUITE_END()
