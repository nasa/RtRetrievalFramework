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

        double model_latitude = 45.945;
        double ref_latitude = 35.0;
        double ref_topo_alt = 15.0;
        //Time reference_time("2005-01-01T00:00:00");

        ref_ap.reset(new ReferenceVmrApriori(model_altitude, model_temperature, model_latitude,
                                             ref_altitude, ref_latitude, ref_topo_alt));
    }
    virtual ~RefVmrFixture() {}

    Array<double, 2> ref_vmr;

    boost::shared_ptr<ReferenceVmrApriori> ref_ap;
};

BOOST_FIXTURE_TEST_SUITE(reference_vmr_apriori, RefVmrFixture)

BOOST_AUTO_TEST_CASE(model_grid_vmr)
{
    // Check tropopause altitude calculation
    BOOST_CHECK_CLOSE(14.392430521241026, ref_ap->model_tropopause_altitude(), 1e-6);

    // Check effective altitude calculation
    IfstreamCs eff_alt_input(test_data_dir() + "expected/reference_vmr_apriori/effective_altitude.dat");
    Array<double, 1> eff_alt;
    eff_alt_input >> eff_alt;
 
    BOOST_CHECK_MATRIX_CLOSE_TOL(ref_ap->effective_altitude(), eff_alt, 1e-4);

    // Check resampling of vmr
    IfstreamCs model_grid_input(test_data_dir() + "expected/reference_vmr_apriori/model_grid_vmr.dat");
    Array<double, 2> ggg_model_grid_vmr;
    model_grid_input >> ggg_model_grid_vmr;

    for(int gas_idx = 0; gas_idx << ggg_model_grid_vmr.rows(); gas_idx++) {
        Array<double, 1> mod_grid_vmr = ref_ap->resample_to_model_grid(ref_vmr(Range::all(), gas_idx));
        BOOST_CHECK_MATRIX_CLOSE(mod_grid_vmr, ggg_model_grid_vmr(Range::all(), gas_idx));
    }
}

BOOST_AUTO_TEST_SUITE_END()
