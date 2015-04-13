#include "forward_model_spectral_grid.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include <fstream>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(forward_model_spectral_grid, ConfigurationFixture)

// Use to print out expected values
void print_grid(std::ostream& out, std::vector<SpectralDomain>& grid)
{
    out << std::setprecision(20);
    for(int spec_idx = 0; spec_idx < grid.size(); spec_idx++) {
        out << "# Spectrometer: " << spec_idx << std::endl
            << grid[spec_idx].data() << std::endl;
    }
}

void print_pixel_list(std::ostream& out, std::vector<std::vector<int> >& plist)
{
    for(int spec_idx = 0; spec_idx < plist.size(); spec_idx++) {
        out << "# Spectrometer: " << spec_idx << std::endl
            << plist[spec_idx].size() << std::endl;
        for(int pix = 0; pix < plist[spec_idx].size(); pix++) {
            out << plist[spec_idx][pix] << " ";
        }
        out << std::endl;
    }
}

void check_expected_grid(std::string filename, std::vector<SpectralDomain>& calc_grids)
{
    IfstreamCs expected_inputs(filename);

    for(int spec_idx = 0; spec_idx < calc_grids.size(); spec_idx++) {
        blitz::Array<double, 1> expt_grid;
        expected_inputs >> expt_grid;

        BOOST_CHECK_MATRIX_CLOSE(expt_grid, calc_grids[spec_idx].data());
    }
}

void check_expected_pixel_list(std::string filename, std::vector<std::vector<int> >& calc_lists)
{
    IfstreamCs expected_inputs(filename);

    for(int spec_idx = 0; spec_idx < calc_lists.size(); spec_idx++) {
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

    //std::ofstream save(expected_file.c_str());
    //print_grid(save, lgrid_calc);
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

    //std::ofstream save(expected_file.c_str());
    //print_grid(save, hgrid_calc);

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

    //std::ofstream save(expected_file.c_str());
    //print_grid(save, hgrid_interp_calc);
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

    //std::ofstream save(expected_file.c_str());
    //print_pixel_list(save, plist_expt);
}

BOOST_AUTO_TEST_SUITE_END()
