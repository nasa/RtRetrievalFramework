#include "gas_vmr_apriori.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"

#include "oco_sounding_id.h"
#include "level_1b_oco.h"
#include "oco_met_file.h"
#include "pressure_sigma.h"
#include "temperature_met.h"
#include "altitude_hydrostatic.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(gas_vmr_apriori, GlobalFixture)

BOOST_AUTO_TEST_CASE(apriori_calc)
{
    std::vector<std::string> test_sounding_id = { "2014090915251774", "2014120112331638", "2015020202201332", "2015020301255031", "2015020400304333" };
    //std::vector<std::string> test_sounding_id = { "2014090915251774" };

    std::ofstream trop_press_file("tropopause_pressure");

    for (auto curr_sid : test_sounding_id) {
        // Match setup used in TcconApriori test
        boost::shared_ptr<HdfFile> sfile(new HdfFile(test_data_dir() + "in/oco2_l1b_apriori_check.h5"));
        boost::shared_ptr<OcoSoundingId> sid (new OcoSoundingId(*sfile, curr_sid));
        boost::shared_ptr<OcoMetFile> ecmwf(new OcoMetFile(test_data_dir() + "in/oco2_ecmwf_apriori_check.h5", sid));
        boost::shared_ptr<Level1bOco> l1b(new Level1bOco(sfile, sid));
        HdfFile hdf_static_input(test_data_dir() + "../input/oco/input/l2_oco_static_input.h5");

        // Create pressure/temp and hence altitude just like done in the production software
        blitz::Array<double, 1> sigma_a = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_a");
        blitz::Array<double, 1> sigma_b = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_b");
        double surface_pressure = ecmwf->surface_pressure();
        boost::shared_ptr<PressureSigma> press
          (new PressureSigma(sigma_a, sigma_b, surface_pressure, false));

        boost::shared_ptr<TemperatureMet> temp
          (new TemperatureMet(ecmwf, press, 0, false));

        boost::shared_ptr<Altitude> alt(new AltitudeHydrostatic(press, temp, l1b->latitude(0), l1b->altitude(0)));

        GasVmrApriori gas_vmr(ecmwf, l1b, alt, hdf_static_input, "/Reference_Atmosphere/", "CO2");

        IfstreamCs expt_vmr_input(test_data_dir() + "expected/gas_vmr_apriori/expt_result");
        Array<double, 1> expt_res;
        expt_vmr_input >> expt_res;

        Array<double, 1> calc_ap_vmr = gas_vmr.apriori_vmr();
        //BOOST_CHECK_MATRIX_CLOSE_TOL(expt_res, calc_ap_vmr, 2e-7);
        trop_press_file << curr_sid << "\t" << gas_vmr.tropopause_pressure() << std::endl;
    }
}

/*  From Brendan's presentation on Reichler method
    Expected tropopause pressure:
    2014090915251774 - 114.6 hPa
    2014120112331638 - 126.6 hPa
    2015020202201332 - 96.0 hPa
    2015020301255031 - 95.3 hPa
    2015020400304333 - 316.7 hPa
*/

BOOST_AUTO_TEST_SUITE_END()
