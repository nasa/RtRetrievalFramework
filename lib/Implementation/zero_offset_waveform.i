// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "zero_offset_waveform.h"
%}
%base_import(sub_state_vector_array)
%base_import(instrument_correction)
%import "array_with_unit.i"
%import "dispersion.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::ZeroOffsetWaveform);

namespace FullPhysics {
class ZeroOffsetWaveform: public SubStateVectorArray<InstrumentCorrection> {
public:
  ZeroOffsetWaveform(double Coeff, 
		     bool Used_flag,
		     const ArrayWithUnit<double, 1>& Zero_offset_waveform,
		     const std::string& Band_name);
  ZeroOffsetWaveform(double Coeff, 
	     bool Used_flag,
	     const Dispersion& Disp,
	     const HdfFile& Hdf_static_input,
	     int Spec_index,
	     const std::string& Band_name,
	     const std::string& Hdf_group = "Instrument/ZeroLevelOffset");
  virtual std::string state_vector_name_i(int i) const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  %python_attribute(zero_offset, ArrayWithUnit<double, 1>)
};
}
