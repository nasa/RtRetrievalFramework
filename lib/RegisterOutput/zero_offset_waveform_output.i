// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "zero_offset_waveform_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "zero_offset_waveform.i"
%fp_shared_ptr(FullPhysics::ZeroOffsetWaveformOutput);

namespace FullPhysics {
class ZeroOffsetWaveformOutput : public RegisterOutputBase {
public:
  ZeroOffsetWaveformOutput(const boost::shared_ptr<ZeroOffsetWaveform>& Z,
			     const std::string& Hdf_band_name);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


