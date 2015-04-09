#include "ground_lambertian.h"
#include "ground_coxmunk.h"
#include "ground_breon.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

class GroundFixture : public GlobalFixture {
public:
    boost::shared_ptr<GroundLambertian> lambertian;
    boost::shared_ptr<GroundCoxmunk> coxmunk;
    boost::shared_ptr<GroundBreonVeg> breon_veg;
    boost::shared_ptr<GroundBreonSoil> breon_soil;

    GroundFixture() {
        // Lambertian
        Array<double, 2> params(3, 2);
        params(0, 0) = 0.5;
        params(0, 1) = 1e-3;

        params(1, 0) = 0.5;
        params(1, 1) = 1e-3;

        params(2, 0) = 0.5;
        params(2, 1) = 1e-3;

        Array<bool, 2> flags(3, 2);
        flags = true;

        blitz::Array<double, 1> ref_points(3);
        ref_points(0) = 0.77;
        ref_points(1) = 1.615;
        ref_points(2) = 2.06;
        ArrayWithUnit<double, 1> ref_points_w_unit(ref_points, units::micron);

        std::vector<std::string> band_name;
        band_name.push_back("ABO2");
        band_name.push_back("WCO2");
        band_name.push_back("SCO2");

        lambertian.reset(new GroundLambertian(params, flags, ref_points_w_unit, band_name));

        // Coxmunk
        Array<double, 1> refr_index(3);
        refr_index(0) = 1.331;
        refr_index(1) = 1.332;
        refr_index(2) = 1.334;

        coxmunk.reset(new GroundCoxmunk(7.1, false, refr_index));

        // Breon vegetative
        double amplitude = 0.1;
        double asymmetry = 0.3;
        double geometric = 1.5;
        breon_veg.reset(new GroundBreonVeg(amplitude, asymmetry, geometric, true, true, true));
        breon_soil.reset(new GroundBreonSoil(amplitude, asymmetry, geometric, true, true, true));

  }
};
