#include "ground_lambertian.h"
#include "ground_brdf.h"
#include "ground_brdf_weight.h"
#include "ground_coxmunk.h"
#include "ground_coxmunk_scaled.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

class GroundFixture : public GlobalFixture {
public:
    boost::shared_ptr<GroundLambertian> lambertian;
    boost::shared_ptr<GroundCoxmunk> coxmunk;
    boost::shared_ptr<GroundBrdfWeight> brdf_weight;
    boost::shared_ptr<GroundBrdfVeg> brdf_veg;
    boost::shared_ptr<GroundBrdfSoil> brdf_soil;

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

        // BrdfWeight
        std::string brdf_name("CoxMunk");
        brdf_weight.reset(new GroundBrdfWeight(params, flags, ref_points_w_unit, band_name, brdf_name));

        // BRDF (Rahman + Breon)
        Array<double, 2> brdf_coeffs(3, 7);
        brdf_coeffs = 
            // rahman_factor, hotspot, asymmetry, anisotropy,
            // breon_factor, weight_intercept, weight_slope
            0.7, 1.0, 1.3, 1.6, 1.9, 0.1, 0.4, 
            0.8, 1.1, 1.4, 1.7, 2.0, 0.2, 0.5, 
	    0.9, 1.2, 1.5, 1.8, 2.1, 0.3, 0.6;
        Array<bool, 2> flag(brdf_coeffs.shape());
        flag = true;
            
        brdf_veg.reset(new GroundBrdfVeg(brdf_coeffs, flag, ref_points_w_unit, band_name));
        brdf_soil.reset(new GroundBrdfSoil(brdf_coeffs, flag, ref_points_w_unit, band_name));

  }
};

