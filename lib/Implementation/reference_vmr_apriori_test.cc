#include "reference_vmr_apriori.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

class RefVmrFixture: public GlobalFixture {
public:
    RefVmrFixture()
    {
        // Altitude grid from modeled data, aka meterological data
        Array<double, 2> model_data;
        IfstreamCs model_data_input(test_data_dir() + "expected/reference_vmr_apriori/model_data.dat");
        model_data_input >> model_data;
        
        Array<double, 1> model_pressure( model_data(Range::all(), 0) * 100 ); // Convert hPa to Pa
        Array<double, 1> model_temperature( model_data(Range::all(), 1) );
        Array<double, 1> model_altitude( model_data(Range::all(), 2) );

        // Reference VMR and altitude grid for a specific date and time
        // For us its 2005 at 35.0 deg latitude
        IfstreamCs ref_vmr_input(test_data_dir() + "expected/reference_vmr_apriori/gnd_summer.vmr.dat");
        Array<double, 2> ref_vmr_data;
        ref_vmr_input >> ref_vmr_data;
        
        Array<double, 1> ref_altitude(ref_vmr_data(Range::all(), 0));

        ref_vmr.resize(ref_vmr_data.rows(), ref_vmr_data.cols() - 1);
        ref_vmr(Range::all(), Range::all()) = ref_vmr_data(Range::all(), Range(1,toEnd));

        double ref_latitude = 35.0;
        double ref_topo_alt = 15.0;
        Time ref_time = Time::parse_time("2005-01-01T00:00:00");

        double obs_latitude = 45.945;
        Time obs_time = Time::parse_time("2004-07-21T20:41:11.484375");

        ref_ap.reset(new ReferenceVmrApriori(model_pressure, model_altitude, model_temperature, 
                                             ref_altitude, ref_latitude, ref_time, ref_topo_alt,
                                             obs_latitude, obs_time));

        IfstreamCs gas_names_input(test_data_dir() + "expected/reference_vmr_apriori/gas_names.txt");
        std::string line;
        while (std::getline(gas_names_input, line)) {
           gas_names.push_back(line);
        }

        BOOST_CHECK_EQUAL(gas_names.size(), ref_vmr.cols());

    }
    virtual ~RefVmrFixture() {}

    Array<double, 2> ref_vmr;
    std::vector<std::string> gas_names;

    boost::shared_ptr<ReferenceVmrApriori> ref_ap;
};

BOOST_FIXTURE_TEST_SUITE(reference_vmr_apriori, RefVmrFixture)

BOOST_AUTO_TEST_CASE(apriori_calc)
{
    // Check tropopause altitude calculation
    BOOST_CHECK_CLOSE(14.392430521241026, ref_ap->model_tropopause_altitude().value, 1e-6);

    // Check effective altitude calculation
    IfstreamCs eff_alt_input(test_data_dir() + "expected/reference_vmr_apriori/effective_altitude.dat");
    Array<double, 1> eff_alt_expt;
    eff_alt_input >> eff_alt_expt;
 
    Array<double, 1> eff_alt_calc = ref_ap->effective_altitude();
    BOOST_CHECK_MATRIX_CLOSE_TOL(eff_alt_calc, eff_alt_expt, 5e-6);

    IfstreamCs model_grid_input(test_data_dir() + "expected/reference_vmr_apriori/model_grid_vmr.dat");
    Array<double, 2> expt_model_grid_vmr;
    model_grid_input >> expt_model_grid_vmr;

    IfstreamCs latitude_grad_input(test_data_dir() + "expected/reference_vmr_apriori/latitude_grad_vmr.dat");
    Array<double, 2> expt_lat_grad_vmr;
    latitude_grad_input >> expt_lat_grad_vmr;

    IfstreamCs secular_trend_input(test_data_dir() + "expected/reference_vmr_apriori/secular_trend_vmr.dat");
    Array<double, 2> expt_secular_trend_vmr;
    secular_trend_input >> expt_secular_trend_vmr;

    IfstreamCs seasonal_cycle_input(test_data_dir() + "expected/reference_vmr_apriori/seasonal_cycle_vmr.dat");
    Array<double, 2> expt_seasonal_cycle_vmr;
    seasonal_cycle_input >> expt_seasonal_cycle_vmr;

    for(int gas_idx = 0; gas_idx < ref_vmr.cols(); gas_idx++) {
        // Resample to model grid, resampling here not implemented 100% the same way as from expected
        // inputs.
        Array<double, 1> mod_grid_vmr = ref_ap->resample_to_model_grid(ref_vmr(Range::all(), gas_idx));
        BOOST_CHECK_MATRIX_CLOSE_TOL(mod_grid_vmr, expt_model_grid_vmr(Range::all(), gas_idx), 1e-4);

        // Apply latitude gradient, use expected from last step due to resampling differences so we
        // can still compare to our expected value
        Array<double, 1> lat_grad_vmr = ref_ap->apply_latitude_gradient(expt_model_grid_vmr(Range::all(), gas_idx), gas_names[gas_idx]);
        BOOST_CHECK_MATRIX_CLOSE_TOL(lat_grad_vmr, expt_lat_grad_vmr(Range::all(), gas_idx), 1e-10);

        // Secular trend
        Array<double, 1> secular_trend_vmr = ref_ap->apply_secular_trend(lat_grad_vmr, gas_names[gas_idx]);
        BOOST_CHECK_MATRIX_CLOSE_TOL(secular_trend_vmr, expt_secular_trend_vmr(Range::all(), gas_idx), 1e-10);

        // Seasonal cycle
        Array<double, 1> seasonal_cycle_vmr = ref_ap->apply_seasonal_cycle(secular_trend_vmr, gas_names[gas_idx]);
        BOOST_CHECK_MATRIX_CLOSE_TOL(seasonal_cycle_vmr, expt_seasonal_cycle_vmr(Range::all(), gas_idx), 1e-10);

        // All steps together, once again differences in resampling contributes to higher level of differences here
        Array<double, 1> apriori_vmr = ref_ap->apriori_vmr(ref_vmr(Range::all(), gas_idx), gas_names[gas_idx]);
        BOOST_CHECK_MATRIX_CLOSE_TOL(apriori_vmr, expt_seasonal_cycle_vmr(Range::all(), gas_idx), 6e-5);
    }

}

BOOST_AUTO_TEST_SUITE_END()
