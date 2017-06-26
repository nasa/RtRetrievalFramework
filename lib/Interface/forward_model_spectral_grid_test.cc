#include "forward_model_spectral_grid.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include <fstream>

#include "nonuniform_spectrum_sampling.h"
#include "uniform_spectrum_sampling.h"

#define OUTPUT_EXPECTED false

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(forward_model_spectral_grid, ConfigurationFixture)

// Use to print out expected values
void print_grid(std::ostream& out, std::vector<SpectralDomain>& grid)
{
  out << std::setprecision(20);
  for(int spec_idx = 0; spec_idx < (int) grid.size(); spec_idx++) {
    out << "# Spectrometer: " << spec_idx << std::endl
	<< grid[spec_idx].data() << std::endl;
  }
}

void print_pixel_list(std::ostream& out, std::vector<std::vector<int> >& plist)
{
  for(int spec_idx = 0; spec_idx < (int) plist.size(); spec_idx++) {
    out << "# Spectrometer: " << spec_idx << std::endl
	<< plist[spec_idx].size() << std::endl;
    for(int pix = 0; pix < (int) plist[spec_idx].size(); pix++) {
      out << plist[spec_idx][pix] << " ";
    }
    out << std::endl;
  }
}

void check_expected_grid(std::string filename, std::vector<SpectralDomain>& calc_grids)
{
    IfstreamCs expected_inputs(filename);

    for(int spec_idx = 0; spec_idx < (int) calc_grids.size(); spec_idx++) {
        blitz::Array<double, 1> expt_grid;
        expected_inputs >> expt_grid;

        BOOST_CHECK_MATRIX_CLOSE(expt_grid, calc_grids[spec_idx].data());
    }
}

void check_expected_pixel_list(std::string filename, std::vector<std::vector<int> >& calc_lists)
{
    IfstreamCs expected_inputs(filename);

    for(int spec_idx = 0; spec_idx < (int) calc_lists.size(); spec_idx++) {
        int num_pix;
        expected_inputs >> num_pix;

        std::vector<int> plist;
        for(int pix = 0; pix << num_pix; pix++) {
            int pix_value;
            expected_inputs >> pix_value;

            BOOST_CHECK_EQUAL(pix_value, calc_lists[spec_idx][pix]);
        }
    }
 }

BOOST_AUTO_TEST_CASE(low_resolution_grid)
{
    std::string expected_file = test_data_dir() + "expected/forward_model_spectral_grid/low_resolution_grid";

    ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, config_spectrum_sampling);

    std::vector<SpectralDomain> lgrid_calc;
    for (int spec_idx = 0; spec_idx < fg.number_spectrometer(); spec_idx++) {
        lgrid_calc.push_back(fg.low_resolution_grid(spec_idx).data());
    }

    check_expected_grid(expected_file, lgrid_calc);

    if (OUTPUT_EXPECTED) {
        std::ofstream save(expected_file.c_str());
        print_grid(save, lgrid_calc);
    }
}

BOOST_AUTO_TEST_CASE(high_resolution_grid)
{
    std::string expected_file = test_data_dir() + "expected/forward_model_spectral_grid/high_resolution_grid";

    ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, config_spectrum_sampling);

    std::vector<SpectralDomain> hgrid_calc;
    for (int spec_idx = 0; spec_idx < fg.number_spectrometer(); spec_idx++) {
        hgrid_calc.push_back(fg.high_resolution_grid(spec_idx).data());
    }
    check_expected_grid(expected_file, hgrid_calc);

    if (OUTPUT_EXPECTED) {
        std::ofstream save(expected_file.c_str());
        print_grid(save, hgrid_calc);
    }

}

BOOST_AUTO_TEST_CASE(high_resolution_interpolated_grid)
{
    std::string expected_file = test_data_dir() + "expected/forward_model_spectral_grid/high_resolution_interpolated_grid";

    ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, config_spectrum_sampling);

    std::vector<SpectralDomain> hgrid_interp_calc;
    for (int spec_idx = 0; spec_idx < fg.number_spectrometer(); spec_idx++) {
        hgrid_interp_calc.push_back(fg.high_resolution_interpolated_grid(spec_idx).data());
    }

    check_expected_grid(expected_file, hgrid_interp_calc);

    if (OUTPUT_EXPECTED) {
        std::ofstream save(expected_file.c_str());
        print_grid(save, hgrid_interp_calc);
    }
}

BOOST_AUTO_TEST_CASE(pixel_grid)
{
    std::string expected_file = test_data_dir() + "expected/forward_model_spectral_grid/pixel_grid";

    ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, config_spectrum_sampling);

    std::vector<std::vector<int> > plist_expt;
    for (int spec_idx = 0; spec_idx < fg.number_spectrometer(); spec_idx++) {
        plist_expt.push_back(fg.pixel_list(spec_idx));
    }

    check_expected_pixel_list(expected_file, plist_expt);

    if (OUTPUT_EXPECTED) {
        std::ofstream save(expected_file.c_str());
        print_pixel_list(save, plist_expt);
    }
}

BOOST_AUTO_TEST_CASE(interpolate_spectrum)
{
    std::string expected_file = test_data_dir() + "expected/forward_model_spectral_grid/interpolate_spectrum";

    std::ofstream save;
    IfstreamCs expected_inputs(expected_file);
    if (OUTPUT_EXPECTED) {
        save.open(expected_file.c_str());
    }

    // Use non uniform sampling so that interpolated_spectrum will actually do something other than
    // return the high res spectrum
    HeritageFile c1(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_abo2_oco__absco_v3.1.0__wn.dat");
    HeritageFile c2(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_wco2_oco__absco_v3.1.0__wn.dat");
    HeritageFile c3(test_data_dir() + "nonunif_rt_grid/nonunif_rt_grid__gosat_sco2_oco__absco_v3.1.0__wn.dat");
    boost::shared_ptr<SpectrumSampling> interpolated_spec(
            new UniformSpectrumSampling(12950.4, 13189.7, 0.05, 
                                        6166.29, 6285.75, 0.05, 
                                        4810.14, 4896.74, 0.05));

    boost::shared_ptr<SpectrumSampling> non_uni_spec_samp( new NonuniformSpectrumSampling(c1,c2,c3, interpolated_spec) );

    ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, non_uni_spec_samp);

    for (int spec_idx = 0; spec_idx < fg.number_spectrometer(); spec_idx++) {
        SpectralDomain lgrid = fg.low_resolution_grid(spec_idx);

        // Create a new grid of some arbitrary size with the end points of the low resolution grid
        Array<double, 1> new_grid(1000);
        Array<double, 1> grid_data(new_grid.rows());

        double delta = (lgrid.data()(lgrid.data().rows()-1) - lgrid.data()(0)) / new_grid.rows();
        for(int gidx = 0; gidx < new_grid.rows(); gidx++) {
            new_grid(gidx) = lgrid.data()(0) + delta * gidx;
            grid_data(gidx) = gidx;
        }

        SpectralDomain spec_dom(new_grid, units::inv_cm);
        SpectralRange spec_range(grid_data, Unit("W / cm^2 / sr / cm^-1"));

        Spectrum spec_in(spec_dom, spec_range);

        Spectrum spec_out = fg.interpolate_spectrum(spec_in, spec_idx);

        if (OUTPUT_EXPECTED) {
            save << std::setprecision(20);
            save << "# Spectrometer: " << spec_idx << ", SpectralDomain" << std::endl;
            save << spec_out.spectral_domain().data() << std::endl;
            save << "# Spectrometer: " << spec_idx << ", SpectralRange" << std::endl;
            save << spec_out.spectral_range().data() << std::endl;
        } else {
            Array<double, 1> expected_spec_domain, expected_spec_range;
            expected_inputs >> expected_spec_domain >> expected_spec_range;
 
            BOOST_CHECK_MATRIX_CLOSE(expected_spec_domain, spec_out.spectral_domain().data());
            BOOST_CHECK_MATRIX_CLOSE(expected_spec_range, spec_out.spectral_range().data());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
