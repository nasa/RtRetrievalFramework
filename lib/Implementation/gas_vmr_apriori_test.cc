#include "gas_vmr_apriori.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"

#include "level_1b_acos.h"
#include "acos_ecmwf.h"
#include "pressure_sigma.h"
#include "temperature_level_offset.h"
#include "altitude_hydrostatic.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(gas_vmr_apriori, GlobalFixture)

BOOST_AUTO_TEST_CASE(apriori_calc)
{
    // Match setup used in TcconApriori test
    boost::shared_ptr<HdfFile> sfile(new HdfFile(test_data_dir() + "in/sounding_id.h5"));
    boost::shared_ptr<AcosSoundingId> sid (new AcosSoundingId(*sfile, "20091009203401", AcosSoundingId::S_SOUNDING));
    boost::shared_ptr<AcosEcmwf> ecmwf(new AcosEcmwf(test_data_dir() + "in/ecmwf.h5", sid, true));
    boost::shared_ptr<Level1bAcos> l1b(new Level1bAcos(sfile, sid));
    
    // Compute model altitudes
    blitz::Array<double, 1> press_profile;
    blitz::Array<double, 1> temp_profile;
    double surface_pressure = ecmwf->surface_pressure();

    ecmwf->temperature_grid(press_profile, temp_profile);
    Array<double, 1> b(press_profile.shape());
    b = 0;
    boost::shared_ptr<PressureSigma> press
      (new PressureSigma(press_profile, b, surface_pressure, false));
    boost::shared_ptr<TemperatureLevelOffset> temp
      (new TemperatureLevelOffset(press, temp_profile, 0, false));

    boost::shared_ptr<Altitude> alt(new AltitudeHydrostatic(press, temp, l1b->latitude(0), l1b->altitude(0)));
    HdfFile hdf_static_input(test_data_dir() + "../input/gosat/input/l2_gosat_static_input.h5");

    GasVmrApriori gas_vmr(ecmwf, l1b, alt, hdf_static_input, "/Reference_Atmosphere/", "CO2");

    std::cerr << "ap_vmr = " << gas_vmr.apriori_vmr() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
