#include "dispersion_polynomial_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> disp_pol_create
(const boost::shared_ptr<Dispersion>& Disp, const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new DispersionPolynomialOutput
     (boost::dynamic_pointer_cast<DispersionPolynomial>(Disp), Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(DispersionPolynomialOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<DispersionPolynomial>&, 
			  std::string>())
.scope
[
 luabind::def("create", &disp_pol_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void DispersionPolynomialOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the instrument state
  boost::shared_ptr<DispersionPolynomial> dfreeze = 
    boost::dynamic_pointer_cast<DispersionPolynomial>(d->clone());
  out->register_data_source
    ("/RetrievalResults/dispersion_offset_apriori_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_offset, dfreeze);
  out->register_data_source
    ("/RetrievalResults/dispersion_spacing_apriori_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_spacing, dfreeze);
}

void DispersionPolynomialOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/dispersion_offset_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_offset, d);
  out->register_data_source
    ("/RetrievalResults/dispersion_spacing_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_spacing, d);
  out->register_data_source
    ("/RetrievalResults/dispersion_offset_uncert_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_offset_uncertainty, d);
  out->register_data_source
    ("/RetrievalResults/dispersion_spacing_uncert_" + hdf_band_name, 
     &DispersionPolynomial::dispersion_spacing_uncertainty, d);
}

