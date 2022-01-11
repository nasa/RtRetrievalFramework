#include "ground_coxmunk_scaled_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_cm_scaled_output_create
                                          (boost::shared_ptr<Level1b>& l1b,
                                           boost::shared_ptr<Ground>& coxmunk,
                                           boost::shared_ptr<Ground>& brdf_weight,
                                           std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundCoxmunkScaledOutput(boost::dynamic_pointer_cast<Level1b>(l1b),
                                       boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk),
                                       boost::dynamic_pointer_cast<GroundBrdfWeight>(brdf_weight),
                                       hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundCoxmunkScaledOutput, RegisterOutputBase)
.def(luabind::constructor<boost::shared_ptr<Level1b>&,
                          boost::shared_ptr<GroundCoxmunk>&,
                          boost::shared_ptr<GroundBrdfWeight>&,
                          std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_cm_scaled_output_create)
]
REGISTER_LUA_END()
#endif

GroundCoxmunkScaledOutput::GroundCoxmunkScaledOutput(
        const boost::shared_ptr<Level1b>& L1b,
        const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
        const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
        const std::vector<std::string>& Hdf_band_names)
: l1b(L1b), coxmunk(Coxmunk), brdf_weight(Brdf_weight), hdf_band_names(Hdf_band_names),
  surface_type("Coxmunk-Scaled")
{
    coxmunk_output.reset(new GroundCoxmunkOutput(Coxmunk));
    brdf_weight_output.reset(new GroundBrdfWeightOutput(Brdf_weight, Hdf_band_names));
}

double kernel_amplitude(boost::shared_ptr<Level1b>& L1b,
                        boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                        int spec_idx)
{
    Unit angle_unit("deg");

    return Coxmunk->kernel_value(spec_idx,
        L1b->solar_zenith(spec_idx).convert(angle_unit).value,
        L1b->sounding_zenith(spec_idx).convert(angle_unit).value,
        L1b->relative_azimuth(spec_idx).convert(angle_unit).value);
}

double kernel_amplitude_giss(boost::shared_ptr<Level1b>& L1b,
                             boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                             int spec_idx)
{
    Unit angle_unit("deg");

    blitz::Array<double, 1> stokes_coeff =
        L1b->stokes_coefficient(spec_idx);

    blitz::Array<double, 1> kernel_value_giss =
        Coxmunk->kernel_value_giss(spec_idx,
            L1b->solar_zenith(spec_idx).convert(angle_unit).value,
            L1b->sounding_zenith(spec_idx).convert(angle_unit).value,
            L1b->relative_azimuth(spec_idx).convert(angle_unit).value);

    double res = 0;

    for (int j = 0; j < kernel_value_giss.rows(); ++j) {
        res += kernel_value_giss(j) * stokes_coeff(j);
    }

    return res;
}

AutoDerivative<double> kernel_amplitude_l_giss(boost::shared_ptr<Level1b>& L1b,
                                               boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                               int spec_idx)
{
    Unit angle_unit("deg");

    blitz::Array<double, 1> stokes_coeff =
        L1b->stokes_coefficient(spec_idx);

    blitz::Array<AutoDerivative<double>, 1> kernel_value_l_giss =
        Coxmunk->kernel_value_l_giss(spec_idx,
            L1b->solar_zenith(spec_idx).convert(angle_unit).value,
            L1b->sounding_zenith(spec_idx).convert(angle_unit).value,
            L1b->relative_azimuth(spec_idx).convert(angle_unit).value);

    AutoDerivative<double> res(kernel_value_l_giss(0) * stokes_coeff(0));

    for (int j = 1; j < kernel_value_l_giss.rows(); ++j) {
        res += kernel_value_l_giss(j) * stokes_coeff(j);
    }

    return res;
}

double reflectance_intercept(boost::shared_ptr<Level1b>& L1b,
                             boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                             boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                             boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                             int spec_idx)
{
  return weight_output->weight_intercept(Brdf_weight, spec_idx) *
         kernel_amplitude_giss(L1b, Coxmunk, spec_idx);
}

double reflectance_slope(boost::shared_ptr<Level1b>& L1b,
                         boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                         boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                         boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                         int spec_idx)
{
  return weight_output->weight_slope(Brdf_weight, spec_idx) *
         kernel_amplitude_giss(L1b, Coxmunk, spec_idx);
}

double reflectance_coeff(boost::shared_ptr<Level1b>& L1b,
                         boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                         boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                         boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                         int spec_idx, int i)
{
  return weight_output->weight_coeff(Brdf_weight, spec_idx, i) *
         kernel_amplitude_giss(L1b, Coxmunk, spec_idx);
}

double reflectance_intercept_uncert(boost::shared_ptr<Level1b>& L1b,
                                    boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                    boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                                    boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                                    int spec_idx)
{
  AutoDerivative<double> kern_amp = kernel_amplitude_l_giss(L1b, Coxmunk, spec_idx);

  return sqrt(pow(weight_output->weight_intercept_uncert(Brdf_weight, spec_idx) * kern_amp.value(), 2) +
              pow(weight_output->weight_intercept       (Brdf_weight, spec_idx) * kern_amp.gradient()(0) *
                  Coxmunk->windspeed_uncert(), 2));
}

double reflectance_slope_uncert(boost::shared_ptr<Level1b>& L1b,
                                boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                                boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                                int spec_idx)
{
  AutoDerivative<double> kern_amp = kernel_amplitude_l_giss(L1b, Coxmunk, spec_idx);

  return sqrt(pow(weight_output->weight_slope_uncert(Brdf_weight, spec_idx) * kern_amp.value(), 2) +
              pow(weight_output->weight_slope       (Brdf_weight, spec_idx) * kern_amp.gradient()(0) *
                  Coxmunk->windspeed_uncert(), 2));
}

double reflectance_coeff_uncert(boost::shared_ptr<Level1b>& L1b,
                                boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                boost::shared_ptr<GroundBrdfWeight>& Brdf_weight,
                                boost::shared_ptr<GroundBrdfWeightOutput>& weight_output,
                                int spec_idx, int i)
{
  AutoDerivative<double> kern_amp = kernel_amplitude_l_giss(L1b, Coxmunk, spec_idx);

  return sqrt(pow(weight_output->weight_coeff_uncert(Brdf_weight, spec_idx, i) * kern_amp.value(), 2) +
              pow(weight_output->weight_coeff       (Brdf_weight, spec_idx, i) * kern_amp.gradient()(0) *
                  Coxmunk->windspeed_uncert(), 2));
}

void GroundCoxmunkScaledOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    coxmunk_output->register_output_apriori(out);
    brdf_weight_output->register_output_apriori(out);

    boost::shared_ptr<GroundCoxmunk> cm_freeze =
      boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk->clone());

    boost::shared_ptr<GroundBrdfWeight> weight_freeze =
      boost::dynamic_pointer_cast<GroundBrdfWeight>(brdf_weight->clone());

    for(int spec_idx = 0; spec_idx < brdf_weight->number_spectrometer(); spec_idx++) {
        std::string band_name = hdf_band_names[spec_idx];

        { boost::function<double ()> f = boost::bind(&reflectance_intercept, l1b, cm_freeze, weight_freeze, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_apriori_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&reflectance_slope, l1b, cm_freeze, weight_freeze, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_slope_apriori_" + band_name, f); }

        for (int i = 2; i < brdf_weight->number_params(); ++i) {
            { boost::function<double ()> f = boost::bind(&reflectance_coeff, l1b, cm_freeze, weight_freeze, brdf_weight_output, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_apriori_" + band_name, f); }
        }
    }
}

void GroundCoxmunkScaledOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    coxmunk_output->register_output(out);
    brdf_weight_output->register_output(out);

    for(int spec_idx = 0; spec_idx < brdf_weight->number_spectrometer(); spec_idx++) {
        std::string band_name = hdf_band_names[spec_idx];

        { boost::function<double ()> f = boost::bind(&reflectance_intercept, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&reflectance_intercept_uncert, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_uncert_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&reflectance_slope, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_slope_" + band_name, f); }

        { boost::function<double ()> f = boost::bind(&reflectance_slope_uncert, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx);
          out->register_data_source("/RetrievalResults/brdf_reflectance_slope_uncert_" + band_name, f); }

        for (int i = 2; i < brdf_weight->number_params(); ++i) {
            { boost::function<double ()> f = boost::bind(&reflectance_coeff, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_" + band_name, f); }

            { boost::function<double ()> f = boost::bind(&reflectance_coeff_uncert, l1b, coxmunk, brdf_weight, brdf_weight_output, spec_idx, i);
              out->register_data_source("/RetrievalResults/brdf_reflectance_" + (i == 2 ? std::string("quadratic") : std::to_string(i)) + "_uncert_" + band_name, f); }
        }
    }

    out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str());
}
