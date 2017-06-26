#include "zero_offset_waveform_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> zow_create
(const boost::shared_ptr<InstrumentCorrection>& Ic, 
 const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new ZeroOffsetWaveformOutput
     (boost::dynamic_pointer_cast<ZeroOffsetWaveform>(Ic), Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(ZeroOffsetWaveformOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &zow_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void ZeroOffsetWaveformOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the instrument state
  boost::shared_ptr<ZeroOffsetWaveform> zfreeze = 
    boost::dynamic_pointer_cast<ZeroOffsetWaveform>(z->clone());
  out->register_data_source
    ("/RetrievalResults/zero_level_offset_apriori_" + hdf_band_name, 
     &ZeroOffsetWaveform::offset, zfreeze);
}

void ZeroOffsetWaveformOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/zero_level_offset_" + hdf_band_name, 
     &ZeroOffsetWaveform::offset, z);
  out->register_data_source
    ("/RetrievalResults/zero_level_offset_uncert_" + hdf_band_name, 
     &ZeroOffsetWaveform::offset_uncertainty, z);
}

