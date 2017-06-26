#include "zero_offset_waveform.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ZeroOffsetWaveform, InstrumentCorrection)
.def(luabind::constructor<double, bool, const Dispersion&, const HdfFile&,
			  int, const std::string&, const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Disp - Dispersion.
/// \param Hdf_static_input - File to read data from.
/// \param Spec_index - Spectral index number to for
/// \param Band_name - Name of band
/// \param Hdf_group - The HDF group to read.
///
/// This read the zero offset waveform from an HDF file. We read the
/// field Hdf_group + "/zero_offset_waveform" +
/// (spec_index + 1),
/// e.g. "Instrument/ZeroLevelOffset/zero_offset_waveform_1".
///
/// This table has two columns, the first is the wavenumber and the
/// second is the zero level offset at that wavenumber. We then do a
/// linear interpolation between wavenumbers to get the full waveform.
///
/// Note that we *require* the Units attribute to be set. This gives
/// the units that the table is in, which may not be the normal SI
/// units we use in the Instrument. 
//-----------------------------------------------------------------------
  
ZeroOffsetWaveform::ZeroOffsetWaveform
(double Coeff, 
 bool Used_flag,
 const Dispersion& Disp,
 const HdfFile& Hdf_static_input,
 int Spec_index,
 const std::string& Band_name,
 const std::string& Hdf_group)
: SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
  band_name(Band_name)
{
  std::string fldname = Hdf_group + "/zero_offset_waveform_" 
    + boost::lexical_cast<std::string>(Spec_index + 1);
  ArrayWithUnit<double, 2> table = 
    Hdf_static_input.read_field_with_unit<double, 2>(fldname);
  // Separate columns
  Array<double, 1> wn(table.value(Range::all(), 0));
  Array<double, 1> waveform(table.value(Range::all(), 1));
  // Set up interpolation
  LinearInterpolate<double, double> 
    waveform_interpolate(wn.begin(), wn.end(), waveform.begin(),
			 LinearInterpolate<double, double>::OUT_OF_RANGE_CLIP);
  Array<double, 1> pixel_grid = Disp.pixel_grid().data();
  zero_offset_.units = table.units;
  zero_offset_.value.resize(pixel_grid.shape());
  for(int i = 0; i < pixel_grid.rows(); ++i)
    zero_offset_.value(i) = waveform_interpolate(pixel_grid(i));
}

boost::shared_ptr<InstrumentCorrection> ZeroOffsetWaveform::clone() const
{
  return boost::shared_ptr<InstrumentCorrection>
    (new ZeroOffsetWaveform(coeff.value()(0), used_flag(0), 
			    zero_offset_, band_name));
}

void ZeroOffsetWaveform::apply_correction
(const SpectralDomain& Pixel_grid,
 const std::vector<int>& Pixel_list,
 SpectralRange& Radiance) const
{
  ArrayWithUnit<double, 1> zoff = zero_offset_.convert(Radiance.units());
  if(Radiance.data_ad().number_variable() == 0) {
    for(int i = 0; i < Radiance.data().rows(); ++i) {
      range_check(Pixel_list[i], 0, zoff.value.rows());
      Radiance.data()(i) += coeff.value()(0) * zoff.value(Pixel_list[i]);
    }
  } else {
    for(int i = 0; i < Radiance.data_ad().rows(); ++i) {
      range_check(Pixel_list[i], 0, zoff.value.rows());
      Radiance.data_ad()(i) = Radiance.data_ad()(i) +
	coeff(0) * zoff.value(Pixel_list[i]);
    }
  }
}

void ZeroOffsetWaveform::print(std::ostream& Os) const
{
  Os << "ZeroOffsetWaveform:\n"
     << "  Band: " << band_name << "\n"
     << "  Scale: " << coeff(0).value() << "\n"
     << "  Fit:   " << (used_flag(0) ? "True\n" : "False\n");
}
