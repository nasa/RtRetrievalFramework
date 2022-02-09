#include "linear_interpolate.h"
#include "old_constant.h"
#include "ground_brdf_output.h"
#include<iostream>
using namespace std;

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_brdf_output_create(const boost::shared_ptr<Ground>& ground, const boost::shared_ptr<Level1b>& l1b, const std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundBrdfOutput(boost::dynamic_pointer_cast<GroundBrdf>(ground), boost::dynamic_pointer_cast<Level1b>(l1b), hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundBrdf>&, const boost::shared_ptr<Level1b>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_brdf_output_create)
]
REGISTER_LUA_END()
#endif


class BrdfOutputHelper {
public:

    BrdfOutputHelper(const boost::shared_ptr<GroundBrdf>& Brdf,
                      const boost::shared_ptr<Level1b>& L1b)
    : brdf(Brdf), l1b(L1b) { }

    double rahman_factor(int spec_idx)
    {
        return brdf->rahman_factor(spec_idx).value();
    }

    double hotspot_parameter(int spec_idx)
    {
        return brdf->hotspot_parameter(spec_idx).value();
    }

    double asymmetry_parameter(int spec_idx)
    {
        return brdf->asymmetry_parameter(spec_idx).value();
    }

    double anisotropy_parameter(int spec_idx)
    {
        return brdf->anisotropy_parameter(spec_idx).value();
    }

    double breon_factor(int spec_idx)
    {
        return brdf->breon_factor(spec_idx).value();
    }

    double uncert_value(int spec_idx, int coeff_idx)
    {
        Array<double, 2> cov = brdf->brdf_covariance(spec_idx);
        if(cov.rows() > 0 and cov.cols() > 0) {
            return (cov(coeff_idx, coeff_idx) < 0 ? 0.0 : sqrt(cov(coeff_idx, coeff_idx)));
        } else {
            return 0.0;
        }
    }

    double weight_intercept(int spec_idx)
    {
        return brdf->weight_intercept(spec_idx).value();
    }

    double weight_slope(int spec_idx)
    {
        return brdf->weight_slope(spec_idx).value();
    }

    double weight_coeff(int spec_idx, int weight_idx)
    {
        return brdf->weight_coeff(spec_idx, weight_idx).value();
    }

    double reflectance_intercept(int spec_idx)
    {
        return weight_intercept(spec_idx) * kernel_amplitude(spec_idx);
    }

    double reflectance_slope(int spec_idx)
    {
        return weight_slope(spec_idx) * kernel_amplitude(spec_idx);
    }

    double reflectance_weight(int spec_idx, int weight_idx)
    {
        return weight_coeff(spec_idx, weight_idx) * kernel_amplitude(spec_idx);
    }

    double rahman_factor_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::RAHMAN_KERNEL_FACTOR_INDEX);
    }

    double hotspot_parameter_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::RAHMAN_OVERALL_AMPLITUDE_INDEX);
    }

    double asymmetry_parameter_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::RAHMAN_ASYMMETRY_FACTOR_INDEX);
    }

    double anisotropy_parameter_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::RAHMAN_GEOMETRIC_FACTOR_INDEX);
    }

    double breon_factor_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::BREON_KERNEL_FACTOR_INDEX);
    }

    double weight_intercept_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::BRDF_WEIGHT_INTERCEPT_INDEX);
    }

    double weight_slope_uncert(int spec_idx)
    {
        return uncert_value(spec_idx, GroundBrdf::BRDF_WEIGHT_SLOPE_INDEX);
    }

    double weight_coeff_uncert(int spec_idx, int weight_idx)
    {
        return uncert_value(spec_idx, 5 + weight_idx);
    }

    double reflectance_intercept_uncert(int spec_idx)
    {
        return weight_intercept_uncert(spec_idx) * kernel_amplitude(spec_idx);
    }

    double reflectance_slope_uncert(int spec_idx)
    {
        return weight_slope_uncert(spec_idx) * kernel_amplitude(spec_idx);
    }

    double reflectance_weight_uncert(int spec_idx, int weight_idx)
    {
        return weight_coeff_uncert(spec_idx, weight_idx) * kernel_amplitude(spec_idx);
    }

    double kernel_amplitude(int spec_idx) const
    {
        Unit angle_unit("deg");
        return brdf->kernel_value(spec_idx,
                l1b->solar_zenith(spec_idx).convert(angle_unit).value,
                l1b->sounding_zenith(spec_idx).convert(angle_unit).value,
                l1b->relative_azimuth(spec_idx).convert(angle_unit).value);
    }
private:
    boost::shared_ptr<GroundBrdf> brdf;
    boost::shared_ptr<Level1b> l1b;
};

void GroundBrdfOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<GroundBrdf> brdf_freeze =
    boost::dynamic_pointer_cast<GroundBrdf>(brdf->clone());
  boost::shared_ptr<BrdfOutputHelper> helper(new BrdfOutputHelper(brdf_freeze, l1b));

  for(int spec_idx = 0; spec_idx < brdf->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::rahman_factor, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::hotspot_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_hotspot_parameter_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::asymmetry_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::anisotropy_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_anisotropy_parameter_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::breon_factor, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_apriori_" + band_name, f); }


      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_intercept, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_slope, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_apriori_" + band_name, f); }

      for (int i = 2; i < brdf->number_weight_parameters(); ++i) {
          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_coeff, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_apriori_" + band_name, f); }
      }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_intercept, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_slope, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_slope_apriori_" + band_name, f); }

      for (int i = 2; i < brdf->number_weight_parameters(); ++i) {
          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_weight, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_apriori_" + band_name, f); }
      }
  }
}

void GroundBrdfOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<BrdfOutputHelper> helper(new BrdfOutputHelper(brdf, l1b));

  for(int spec_idx = 0; spec_idx < brdf->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::rahman_factor, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::rahman_factor_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::hotspot_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_hotspot_parameter_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::hotspot_parameter_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_hotspot_parameter_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::asymmetry_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::asymmetry_parameter_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::anisotropy_parameter, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_anisotropy_parameter_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::anisotropy_parameter_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_anisotropy_parameter_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::breon_factor, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::breon_factor_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_uncert_" + band_name, f); }


      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_intercept, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_intercept_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_slope, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_slope_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_uncert_" + band_name, f); }

      for (int i = 2; i < brdf->number_weight_parameters(); ++i) {
          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_coeff, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_" + band_name, f); }

          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::weight_coeff_uncert, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_uncert_" + band_name, f); }
      }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_intercept, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_intercept_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_slope, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_slope_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_slope_uncert, helper, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_reflectance_slope_uncert_" + band_name, f); }

      for (int i = 2; i < brdf->number_weight_parameters(); ++i) {

          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_weight, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_" + band_name, f); }

          { boost::function<double ()> f = boost::bind(&BrdfOutputHelper::reflectance_weight_uncert, helper, spec_idx, i);
            out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_uncert_" + band_name, f); }
      }
  }

  surface_type = "BRDF " + brdf->breon_type();
  out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str());
}
