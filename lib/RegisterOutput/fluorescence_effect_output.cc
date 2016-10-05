#include "fluorescence_effect_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> fe_create(boost::shared_ptr<SpectrumEffect>& Fe) 
{
  return boost::shared_ptr<RegisterOutputBase>
    (new FluorescenceEffectOutput(boost::dynamic_pointer_cast<FluorescenceEffect>(Fe)));
}
REGISTER_LUA_DERIVED_CLASS(FluorescenceEffectOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &fe_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void FluorescenceEffectOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<FluorescenceEffect> fluor_freeze = 
    boost::dynamic_pointer_cast<FluorescenceEffect>(fluor->clone());
  out->register_data_source
    ("/RetrievalResults/fluorescence_at_reference_apriori", 
     &FluorescenceEffect::fluorescence_at_reference, fluor_freeze);
  out->register_data_source
    ("/RetrievalResults/fluorescence_slope_apriori", 
     &FluorescenceEffect::fluorescence_slope, fluor_freeze);
}

void FluorescenceEffectOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/fluorescence_at_reference",
     &FluorescenceEffect::fluorescence_at_reference, fluor);
  out->register_data_source
    ("/RetrievalResults/fluorescence_at_reference_uncert",
     &FluorescenceEffect::fluorescence_at_reference_uncertainty, fluor);

  out->register_data_source
    ("/RetrievalResults/fluorescence_slope",
     &FluorescenceEffect::fluorescence_slope, fluor);
  out->register_data_source
    ("/RetrievalResults/fluorescence_slope_uncert",
     &FluorescenceEffect::fluorescence_slope_uncertainty, fluor);
}

