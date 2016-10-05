#include "stokes_coefficient_fraction_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> sc_create
(const boost::shared_ptr<StokesCoefficient>& Sc, 
 int Spec_index,
 const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new StokesCoefficientFractionOutput
     (boost::dynamic_pointer_cast<StokesCoefficientFraction>(Sc), Spec_index, Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(StokesCoefficientFractionOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &sc_create)
]
REGISTER_LUA_END()
#endif

class StokesCoefficientFractionOutputHelper {
public:
  StokesCoefficientFractionOutputHelper(const boost::shared_ptr<StokesCoefficientFraction>& Sc, int Spec_index) : f_(Sc), spec_index(Spec_index) {}
  double f() const { return f_->coefficient().value()(spec_index); }
  double f_uncertainty() const {return f_->f_uncertainty(spec_index); }
private:
  boost::shared_ptr<StokesCoefficientFraction> f_;
  int spec_index;
};

// See base class for description

void StokesCoefficientFractionOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the state
  boost::shared_ptr<StokesCoefficientFraction> ffreeze = 
    boost::dynamic_pointer_cast<StokesCoefficientFraction>(f->clone());
  boost::shared_ptr<StokesCoefficientFractionOutputHelper>
    h(new StokesCoefficientFractionOutputHelper(ffreeze, spec_index));
  out->register_data_source
    ("/RetrievalResults/pol_frac_apriori_" + hdf_band_name, 
     &StokesCoefficientFractionOutputHelper::f, h);
}

void StokesCoefficientFractionOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<StokesCoefficientFractionOutputHelper>
    h(new StokesCoefficientFractionOutputHelper(f, spec_index));
  out->register_data_source
    ("/RetrievalResults/pol_frac_" + hdf_band_name, 
     &StokesCoefficientFractionOutputHelper::f, h);
  out->register_data_source
    ("/RetrievalResults/pol_frac_uncert_" + hdf_band_name, 
     &StokesCoefficientFractionOutputHelper::f_uncertainty, h);
}

