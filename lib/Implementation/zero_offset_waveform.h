#ifndef ZERO_OFFSET_WAVEFORM_H
#define ZERO_OFFSET_WAVEFORM_H
#include "dispersion.h"
#include "instrument_correction.h"
#include "hdf_file.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This class is a zero level offset. We use a supplied waveform,
  multiple by a single scale factor given by the state vector, and add
  this to the radiance calculated in InstrumentIls.
*******************************************************************/

class ZeroOffsetWaveform : public SubStateVectorArray<InstrumentCorrection> {
public:
//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Zero_offset_waveform - Offset to add for zero level. This
///    is indexed by the instrument pixel.
/// \param Band_name - Name of band
//-----------------------------------------------------------------------
  
  ZeroOffsetWaveform(double Coeff, 
		     bool Used_flag,
		     const ArrayWithUnit<double, 1>& Zero_offset_waveform,
		     const std::string& Band_name)
    : SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
      band_name(Band_name),
      zero_offset_(Zero_offset_waveform)
  { }

  ZeroOffsetWaveform(double Coeff, 
	     bool Used_flag,
	     const Dispersion& Disp,
	     const HdfFile& Hdf_static_input,
	     int Spec_index,
	     const std::string& Band_name,
	     const std::string& Hdf_group = "Instrument/ZeroLevelOffset");
  virtual ~ZeroOffsetWaveform() {}
  virtual std::string state_vector_name_i(int i) const
  { return "Zero offset waveform scale factor " + band_name; }
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Current value of zero level offset, for each pixel number.
//-----------------------------------------------------------------------
  ArrayWithUnit<double, 1> zero_offset() const
  {
    return ArrayWithUnit<double, 1>
      (blitz::Array<double, 1>(coeff(0).value() * zero_offset_.value), 
       zero_offset_.units);
  }

//-----------------------------------------------------------------------
/// Zero level offset. This is just coeff(0), but we wrap this for use
/// by ZeroOffsetWaveformOutput
//-----------------------------------------------------------------------

  double offset() const { return coeff.value()(0); }

//-----------------------------------------------------------------------
/// Zero level offset uncertainty. This is just sqrt(Cov(0,0)), but we
/// wrap this for use  by ZeroOffsetWaveformOutput
//-----------------------------------------------------------------------

  double offset_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(0,0);
    return (t < 0 ? 0 : sqrt(t)); 
  }

private:
  std::string band_name;
  ArrayWithUnit<double, 1> zero_offset_;
};
}
#endif
