#include "nonuniform_spectrum_sampling.h"
#include "uniform_spectrum_sampling.h"
#include "unit_test_support.h"
#include "fp_exception.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(nonuniform_spectrum_sampling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  HeritageFile c1(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid_abo2.dat");
  HeritageFile c2(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid_wco2.dat");
  HeritageFile c3(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid_sco2.dat");
  HeritageFile c3_srtd(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid_sco2_sorted.dat");
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec(new UniformSpectrumSampling(31.8,34.5,0.05,
						  1.5, 10.5, 0.05,
						  21.7, 26.3, 0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec1(new UniformSpectrumSampling(31.8,34.5,0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec2(new UniformSpectrumSampling(1.5, 10.5, 0.05));
  boost::shared_ptr<SpectrumSampling> 
    interpolated_spec3(new UniformSpectrumSampling(21.7, 26.3, 0.05));
  NonuniformSpectrumSampling nonunif_rt_grid(c1,c2,c3, interpolated_spec);
  NonuniformSpectrumSampling nonunif_rt_grid_1(c1, interpolated_spec1);
  NonuniformSpectrumSampling nonunif_rt_grid_2(c2, interpolated_spec2);
  NonuniformSpectrumSampling nonunif_rt_grid_3_srtd(c3_srtd, interpolated_spec3);

  if(false) {			// Can turn on debugging messages if desired.
    std::cout << nonunif_rt_grid << "\n"
	      << nonunif_rt_grid_1 << "\n"
	      << nonunif_rt_grid_2 << "\n"
	      << nonunif_rt_grid_3_srtd << "\n";
  }

  SpectralDomain dummy;		// Not used by UniformSpectrumSampling
  DoubleWithUnit dummy2;
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
		    spectral_domain(0, dummy, dummy2).wavenumber().size(), 7);
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
		    spectral_domain(1, dummy, dummy2).wavenumber().size(), 10);
  BOOST_CHECK_EQUAL(nonunif_rt_grid.
		    spectral_domain(2, dummy, dummy2).wavenumber().size(), 9);

  BOOST_CHECK_EQUAL(nonunif_rt_grid_1.
		    spectral_domain(0, dummy, dummy2).wavenumber().size(), 7);
  BOOST_CHECK_EQUAL(nonunif_rt_grid_2.
		    spectral_domain(0, dummy, dummy2).wavenumber().size(), 10);
  BOOST_CHECK_EQUAL(nonunif_rt_grid_3_srtd.
		    spectral_domain(0, dummy, dummy2).wavenumber().size(), 9);

  BOOST_CHECK_EQUAL(sum(abs(nonunif_rt_grid.
		    spectral_domain(2, dummy, dummy2).wavenumber()-
		    nonunif_rt_grid_3_srtd.
		    spectral_domain(0, dummy, dummy2).wavenumber())), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
