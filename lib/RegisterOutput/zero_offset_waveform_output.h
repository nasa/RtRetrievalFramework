#ifndef ZERO_OFFSET_WAVEFORM_OUTPUT_H
#define ZERO_OFFSET_WAVEFORM_OUTPUT_H
#include "register_output_base.h"
#include "zero_offset_waveform.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ZeroOffsetWaveform class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ZeroOffsetWaveform class.
*******************************************************************/
class ZeroOffsetWaveformOutput : public RegisterOutputBase {
public:
  ZeroOffsetWaveformOutput(const boost::shared_ptr<ZeroOffsetWaveform>& Z,
			     const std::string& Hdf_band_name) 
    : z(Z), 
      hdf_band_name(Hdf_band_name) {}
  virtual ~ZeroOffsetWaveformOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<ZeroOffsetWaveform> z;
  std::string hdf_band_name;
};
}
#endif
