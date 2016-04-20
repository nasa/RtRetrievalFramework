#include "ground_brdf_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_brdf_output_create(const boost::shared_ptr<Ground>& ground, const std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundBrdfOutput(boost::dynamic_pointer_cast<GroundBrdf>(ground), hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundBrdf>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_brdf_output_create)
]
REGISTER_LUA_END()
#endif

double weight_intercept(boost::shared_ptr<GroundBrdf>& Brdf, int spec_idx)
{
    return Brdf->weight_intercept(spec_idx).value();
}

double weight_slope(boost::shared_ptr<GroundBrdf>& Brdf, int spec_idx)
{
    return Brdf->weight_slope(spec_idx).value();
}

double rahman_factor(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return Brdf->rahman_factor(spec_idx).value();
}

double overall_amplitude(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx) 
{
    return Brdf->overall_amplitude(spec_idx).value();
}

double asymmetry_parameter(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return Brdf->asymmetry_parameter(spec_idx).value();
}

double geometric_factor(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return Brdf->geometric_factor(spec_idx).value();
}

double breon_factor(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return Brdf->breon_factor(spec_idx).value();
}

double uncert_value(boost::shared_ptr<GroundBrdf>& Brdf, int spec_idx, int coeff_idx)
{
    Array<double, 2> cov = Brdf->brdf_covariance(spec_idx);
    if(cov.rows() > 0 and cov.cols() > 0) {
        return (cov(coeff_idx, coeff_idx) < 0 ? 0.0 : sqrt(cov(coeff_idx, coeff_idx)));
    } else {
        return 0.0;
    }
}

double weight_intercept_uncert(boost::shared_ptr<GroundBrdf>& Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::BRDF_WEIGHT_INTERCEPT_INDEX);
}

double weight_slope_uncert(boost::shared_ptr<GroundBrdf>& Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::BRDF_WEIGHT_SLOPE_INDEX);
}

double rahman_factor_uncert(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::RAHMAN_KERNEL_FACTOR_INDEX);
}

double overall_amplitude_uncert(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx) 
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::RAHMAN_OVERALL_AMPLITUDE_INDEX);
}

double asymmetry_parameter_uncert(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::RAHMAN_ASYMMETRY_FACTOR_INDEX);
}

double geometric_factor_uncert(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::RAHMAN_GEOMETRIC_FACTOR_INDEX);
}

double breon_factor_uncert(boost::shared_ptr<GroundBrdf> Brdf, int spec_idx)
{
    return uncert_value(Brdf, spec_idx, GroundBrdf::BREON_KERNEL_FACTOR_INDEX);
}

void GroundBrdfOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<GroundBrdf> brdf_freeze = 
    boost::dynamic_pointer_cast<GroundBrdf>(brdf->clone());

  for(int spec_idx = 0; spec_idx < brdf->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&weight_intercept, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_intercept_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&weight_slope, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&rahman_factor, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&overall_amplitude, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_overall_amplitude_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&asymmetry_parameter, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&geometric_factor, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_geometric_factor_apriori_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&breon_factor, brdf_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_apriori_" + band_name, f); }
  }
}

void GroundBrdfOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  for(int spec_idx = 0; spec_idx < brdf->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&weight_intercept, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_intercept_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&weight_intercept_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_intercept_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&weight_slope, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&weight_slope_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_weight_slope_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&rahman_factor, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&rahman_factor_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_rahman_factor_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&overall_amplitude, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_overall_amplitude_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&overall_amplitude_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_overall_amplitude_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&asymmetry_parameter, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&asymmetry_parameter_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_asymmetry_parameter_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&geometric_factor, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_geometric_factor_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&geometric_factor_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_geometric_factor_uncert_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&breon_factor, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&breon_factor_uncert, brdf, spec_idx);
        out->register_data_source("/RetrievalResults/brdf_breon_factor_uncert_" + band_name, f); }

  }

  surface_type = "BRDF " + brdf->breon_type();
  out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str());
    
}
