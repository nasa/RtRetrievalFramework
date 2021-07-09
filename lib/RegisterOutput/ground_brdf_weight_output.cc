#include "ground_brdf_weight_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_brdf_weight_output_create(const boost::shared_ptr<Ground>& ground, const std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundBrdfWeightOutput(boost::dynamic_pointer_cast<GroundBrdfWeight>(ground), hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfWeightOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundBrdfWeight>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_brdf_weight_output_create)
]
REGISTER_LUA_END()
#endif

double weight_coeff_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i)
{
    return Brdf_weight->weight_coefficients(spec_idx).value()(i);
}

double weight_intercept_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_coeff_(Brdf_weight, spec_idx, 0);
}

double weight_slope_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_coeff_(Brdf_weight, spec_idx, 1);
}

double weight_coeff_uncert_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i)
{
    Array<double, 2> cov = Brdf_weight->weight_covariance(spec_idx);
    if (cov.rows() > i and cov.cols() > i) {
      return (cov(i, i) < 0 ? 0.0 : sqrt(cov(i, i)));
    } else {
      return 0.0;
    }
}

double weight_intercept_uncert_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_coeff_uncert_(Brdf_weight, spec_idx, 0);
}

double weight_slope_uncert_(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_coeff_uncert_(Brdf_weight, spec_idx, 1);
}

double GroundBrdfWeightOutput::weight_intercept(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_intercept_(Brdf_weight, spec_idx);
}

double GroundBrdfWeightOutput::weight_slope(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_slope_(Brdf_weight, spec_idx);
}

double GroundBrdfWeightOutput::weight_coeff(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i)
{
    return weight_coeff_(Brdf_weight, spec_idx, i);
}

double GroundBrdfWeightOutput::weight_intercept_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_intercept_uncert_(Brdf_weight, spec_idx);
}

double GroundBrdfWeightOutput::weight_slope_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx)
{
    return weight_slope_uncert_(Brdf_weight, spec_idx);
}

double GroundBrdfWeightOutput::weight_coeff_uncert(boost::shared_ptr<GroundBrdfWeight>& Brdf_weight, int spec_idx, int i)
{
    return weight_coeff_uncert_(Brdf_weight, spec_idx, i);
}

void GroundBrdfWeightOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    boost::shared_ptr<GroundBrdfWeight> weight_freeze =
      boost::dynamic_pointer_cast<GroundBrdfWeight>(brdf_weight->clone());

    for(int spec_idx = 0; spec_idx < brdf_weight->number_spectrometer(); spec_idx++) {
        std::string band_name = hdf_band_names[spec_idx];

        { boost::function<double ()> f = boost::bind(&weight_intercept_, weight_freeze, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_apriori_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&weight_slope_, weight_freeze, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_slope_apriori_" + band_name, f); }

        for (int i = 2; i < brdf_weight->number_params(); ++i) {
            { boost::function<double ()> f = boost::bind(&weight_coeff_, weight_freeze, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_apriori_" + band_name, f); }
        }
    }
}

void GroundBrdfWeightOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    for(int spec_idx = 0; spec_idx < brdf_weight->number_spectrometer(); spec_idx++) {
        std::string band_name = hdf_band_names[spec_idx];

        { boost::function<double ()> f = boost::bind(&weight_intercept_, brdf_weight, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&weight_intercept_uncert_, brdf_weight, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_uncert_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&weight_slope_, brdf_weight, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_slope_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&weight_slope_uncert_, brdf_weight, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_weight_slope_uncert_" + band_name, f); }

        for (int i = 2; i < brdf_weight->number_params(); ++i) {
            { boost::function<double ()> f = boost::bind(&weight_coeff_, brdf_weight, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_" + band_name, f); }

            { boost::function<double ()> f = boost::bind(&weight_coeff_uncert_, brdf_weight, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_weight_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_uncert_" + band_name, f); }
        }
    }

    out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str());
}
