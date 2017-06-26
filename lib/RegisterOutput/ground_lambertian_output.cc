#include "ground_lambertian_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_lamb_output_create(const boost::shared_ptr<Ground>& ground, const std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundLambertianOutput(boost::dynamic_pointer_cast<GroundLambertian>(ground), hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundLambertianOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundLambertian>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_lamb_output_create)
]
REGISTER_LUA_END()
#endif

double albedo(boost::shared_ptr<GroundLambertian>& Lambertian, int spec_idx)
{
  return Lambertian->albedo_coefficients(spec_idx).value()(0);
}

double albedo_slope(boost::shared_ptr<GroundLambertian>& Lambertian, int spec_idx)
{
  return Lambertian->albedo_coefficients(spec_idx).value()(1);
}

double albedo_uncert(boost::shared_ptr<GroundLambertian>& Lambertian, int spec_idx)
{
  Array<double, 2> cov = Lambertian->albedo_covariance(spec_idx);
  if(cov.rows() > 0 and cov.cols() > 0) {
    return (cov(0, 0) < 0 ? 0.0 : sqrt(cov(0, 0)));
  } else {
    return 0.0;
  }
}

double albedo_slope_uncert(boost::shared_ptr<GroundLambertian>& Lambertian, int spec_idx)
{
  Array<double, 2> cov = Lambertian->albedo_covariance(spec_idx);
  if (cov.rows() > 1 and cov.cols() > 1) {
    return (cov(1, 1) < 0 ? 0.0 : sqrt(cov(1, 1)));
  } else {
    return 0.0;
  }
}

void GroundLambertianOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<GroundLambertian> lamb_freeze = 
    boost::dynamic_pointer_cast<GroundLambertian>(lambertian->clone());

  for(int spec_idx = 0; spec_idx < lambertian->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&albedo, lamb_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_apriori_" + band_name + "_fph", f); }

      { boost::function<double ()> f = boost::bind(&albedo_slope, lamb_freeze, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_slope_apriori_" + band_name, f); }
  }
}

void GroundLambertianOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  for(int spec_idx = 0; spec_idx < lambertian->number_spectrometer(); spec_idx++) {
      std::string band_name = hdf_band_names[spec_idx];

      { boost::function<double ()> f = boost::bind(&albedo, lambertian, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_" + band_name + "_fph", f); }

      { boost::function<double ()> f = boost::bind(&albedo_slope, lambertian, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_slope_" + band_name, f); }

      { boost::function<double ()> f = boost::bind(&albedo_uncert, lambertian, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_uncert_" + band_name + "_fph", f); }

      { boost::function<double ()> f = boost::bind(&albedo_slope_uncert, lambertian, spec_idx);
        out->register_data_source("/RetrievalResults/albedo_slope_uncert_" + band_name, f); }
  }

  out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str());
    
}
