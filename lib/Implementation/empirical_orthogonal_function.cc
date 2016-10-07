#include "empirical_orthogonal_function.h"
#include "linear_interpolate.h"
#include <boost/make_shared.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Luabind is restricted to 10 arguments in a function. The
// constructor for EmpiricalOrthogonalFunction takes too many. As an easy
// workaround, just stuff a number of the values into an array.
boost::shared_ptr<InstrumentCorrection> eof_create(
bool Used_flag,
const HdfFile& Hdf_static_input,
const ArrayWithUnit<double, 1>& Uncertainty,
const std::string& Band_name,
const std::string& Hdf_group,
const blitz::Array<double, 1>& Val)
{
  return boost::make_shared<EmpiricalOrthogonalFunction>
    (Val(0), Used_flag, Hdf_static_input, Uncertainty,
     static_cast<int>(Val(1)), static_cast<int>(Val(2)),
     static_cast<int>(Val(3)),
     Band_name, Hdf_group, Val(4));
}
REGISTER_LUA_DERIVED_CLASS(EmpiricalOrthogonalFunction, InstrumentCorrection)
.def(luabind::constructor<double, bool, const Dispersion&, const HdfFile&,
     int, int, int, const std::string&, const std::string&>())
.def(luabind::constructor<double, bool, const HdfFile&,
     int, int, int, const std::string&, const std::string&>())
.scope
[
 luabind::def("create", &eof_create)
]
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
/// \param Sounding_number - The footprint index (e.g., 0 to 7 for
///     OCO). The EOF might be indexed by sounding number, or it might 
///     not. If it isn't
/// \param Order - Order of the eigenvector (e.g., first order
///    correction, second order correction, etc.)
/// \param Band_name - Name of band
/// \param Hdf_group - The HDF group to read.
///
/// This read the empirical orthogonal function (EOF) from an HDF
/// file. We read the field Hdf_group + "/EOF" + order + "_waveform" +
/// (spec_index + 1),
/// e.g. "Instrument/EmpiricalOrthogonalFunction/EOF_1_waveform_1".
///
/// This table has two columns, the first is the wavenumber and the
/// second is the eof at that wavenumber. We then do a
/// linear interpolation between wavenumbers to get the full waveform.
///
/// Note that we *require* the Units attribute to be set. This gives
/// the units that the table is in, which may not be the normal SI
/// units we use in the Instrument. 
//-----------------------------------------------------------------------
  
EmpiricalOrthogonalFunction::EmpiricalOrthogonalFunction
(double Coeff, 
 bool Used_flag,
 const Dispersion& Disp,
 const HdfFile& Hdf_static_input,
 int Spec_index,
 int Sounding_number,
 int Order,
 const std::string& Band_name,
 const std::string& Hdf_group)
: SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
  band_name(Band_name),
  hdf_group(Hdf_group),
  order_(Order),
  sounding_number_(Sounding_number),
  eof_scale_uncertainty_(false),
  scale_to_stddev_(-1)
{
  using namespace H5;
  std::string fldname = hdf_group + "/EOF_"
    + boost::lexical_cast<std::string>(Order) +
    "_waveform_" 
    + boost::lexical_cast<std::string>(Spec_index + 1);
  ArrayWithUnit<double, 2> table;
  // Determine if we have a sounding number index.
  DataSet d = Hdf_static_input.h5_file().openDataSet(fldname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() == 3) {
    ArrayWithUnit<double, 3> full_table =
      Hdf_static_input.read_field_with_unit<double, 3>(fldname);
    table.units = full_table.units;
    table.value.reference(full_table.value(Sounding_number, Range::all(),
					   Range::all()));
    eof_depend_on_sounding_number_ = true;
  } else {
    table = Hdf_static_input.read_field_with_unit<double, 2>(fldname);
    eof_depend_on_sounding_number_ = false;
  }
  // See the we have the Wave Units. If not, we default to
  // "cm^-1". This is for backwards compatibility.
  Unit wave_unit("cm^-1");
  if(Hdf_static_input.has_attribute(fldname + "/WaveUnits"))
    wave_unit = Unit(Hdf_static_input.read_attribute<std::string>(fldname + "/WaveUnits"));
    
  // Separate columns
  Array<double, 1> wn(table.value(Range::all(), 0));
  Array<double, 1> waveform(table.value(Range::all(), 1));
  // Set up interpolation
  LinearInterpolate<double, double> 
    waveform_interpolate(wn.begin(), wn.end(), waveform.begin(),
			 LinearInterpolate<double, double>::OUT_OF_RANGE_CLIP);
  Array<double, 1> pixel_grid = Disp.pixel_grid().convert_wave(wave_unit);
  eof_.units = table.units;
  eof_.value.resize(pixel_grid.shape());
  for(int i = 0; i < pixel_grid.rows(); ++i)
    eof_.value(i) = waveform_interpolate(pixel_grid(i));
}

//-----------------------------------------------------------------------
/// Constructor. This is similar to the other constructor, except we
/// are reading data where the EOF is *not* indexed by
/// wavenumber. Instead, we have the EOF supplied for every
/// sample_index. 
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Hdf_static_input - File to read data from.
/// \param Spec_index - Spectral index number to for
/// \param Sounding_number - The footprint index (e.g., 0 to 7 for
///     OCO). The EOF might be indexed by sounding number, or it might 
///     not. If it isn't
/// \param Order - Order of the eigenvector (e.g., first order
///    correction, second order correction, etc.)
/// \param Band_name - Name of band
/// \param Hdf_group - The HDF group to read.
///
/// This read the empirical orthogonal function (EOF) from an HDF
/// file. We read the field Hdf_group + "/EOF" + order + "_waveform" +
/// (spec_index + 1),
/// e.g. "Instrument/EmpiricalOrthogonalFunction/EOF_1_waveform_1".
///
/// This table has one columns, giving the eof at each pixel. Or, the
/// table has 2 columns where the second one is the sounding_number
/// (so it depends on the footprint).
///
/// Note that we *require* the Units attribute to be set. This gives
/// the units that the table is in, which may not be the normal SI
/// units we use in the Instrument. 
//-----------------------------------------------------------------------
  
EmpiricalOrthogonalFunction::EmpiricalOrthogonalFunction
(double Coeff, 
 bool Used_flag,
 const HdfFile& Hdf_static_input,
 int Spec_index,
 int Sounding_number,
 int Order,
 const std::string& Band_name,
 const std::string& Hdf_group)
: SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
  band_name(Band_name),
  hdf_group(Hdf_group),
  order_(Order),
  sounding_number_(Sounding_number),
  eof_scale_uncertainty_(false),
  scale_to_stddev_(-1)
{
  using namespace H5;
  std::string fldname = hdf_group + "/EOF_"
    + boost::lexical_cast<std::string>(Order) +
    "_waveform_" 
    + boost::lexical_cast<std::string>(Spec_index + 1);
  // Determine if we have a sounding number index.
  DataSet d = Hdf_static_input.h5_file().openDataSet(fldname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() == 2) {
    ArrayWithUnit<double, 2> full_table =
      Hdf_static_input.read_field_with_unit<double, 2>(fldname);
    eof_.units = full_table.units;
    eof_.value.reference(full_table.value(Sounding_number, Range::all()));
    eof_depend_on_sounding_number_ = true;
  } else {
    eof_ = Hdf_static_input.read_field_with_unit<double, 1>(fldname);
    eof_depend_on_sounding_number_ = false;
  }
}

//-----------------------------------------------------------------------
/// Constructor. This is similar to the other constructor, except we
/// are reading data where the EOF is *not* indexed by
/// wavenumber. Instead, we have the EOF supplied for every
/// sample_index.
///
/// This variation scales the EOF waveform by the supplied
/// uncertainty, and then scales the resulting waveform so the
/// standard deviation is Sigma_scale. This more closely matches what
/// L2 actually fits for.
///
/// In this case, the EOF in the HDF file is unitless, and we get the
/// units from the uncertainty.
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Hdf_static_input - File to read data from.
/// \param Uncertainty - The uncertainty to scale EOF by. This
///    contains units.
/// \param Spec_index - Spectral index number to for
/// \param Sounding_number - The footprint index (e.g., 0 to 7 for
///     OCO). The EOF might be indexed by sounding number, or it might 
///     not. If it isn't
/// \param Order - Order of the eigenvector (e.g., first order
///    correction, second order correction, etc.)
/// \param Band_name - Name of band
/// \param Hdf_group - The HDF group to read.
/// \param Scale_to_stddev - Number to scale to standard deviation of EOF
///    is this value.
///
/// This read the empirical orthogonal function (EOF) from an HDF
/// file. We read the field Hdf_group + "/EOF" + order + "_waveform" +
/// (spec_index + 1),
/// e.g. "Instrument/EmpiricalOrthogonalFunction/EOF_1_waveform_1".
/// This is then multiplied by the supplied Uncertainty, and then
/// scaled so the standard deviation is Scale_to_stddev.
///
/// This table has one columns, giving the eof at each pixel. Or, the
/// table has 2 columns where the second one is the sounding_number
/// (so it depends on the footprint).
//-----------------------------------------------------------------------
  
EmpiricalOrthogonalFunction::EmpiricalOrthogonalFunction
(double Coeff, 
 bool Used_flag,
 const HdfFile& Hdf_static_input,
 const ArrayWithUnit<double, 1>& Uncertainty,
 int Spec_index,
 int Sounding_number,
 int Order,
 const std::string& Band_name,
 const std::string& Hdf_group,
 double Scale_to_stddev)
: SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag),
  band_name(Band_name),
  hdf_group(Hdf_group),
  order_(Order),
  sounding_number_(Sounding_number),
  eof_scale_uncertainty_(true),
  scale_to_stddev_(Scale_to_stddev)
{
  using namespace H5;
  std::string fldname = hdf_group + "/EOF_"
    + boost::lexical_cast<std::string>(Order) +
    "_waveform_" 
    + boost::lexical_cast<std::string>(Spec_index + 1);
  // Determine if we have a sounding number index.
  DataSet d = Hdf_static_input.h5_file().openDataSet(fldname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() == 2) {
    Array<double, 2> full_table =
      Hdf_static_input.read_field<double, 2>(fldname);
    eof_.units = Uncertainty.units;
    eof_.value.reference(full_table(Sounding_number, Range::all()));
    eof_depend_on_sounding_number_ = true;
  } else {
    eof_.units = Uncertainty.units;
    eof_.value.reference(Hdf_static_input.read_field<double, 1>(fldname));
    eof_depend_on_sounding_number_ = false;
  }
  eof_.value = eof_.value * Uncertainty.value;
  // Scale eof so stddev is the given value
  double s = scale_to_stddev_ / sqrt(sum(sqr(eof_.value - mean(eof_.value))) /
				     (eof_.value.size() - 1));
  eof_.value = eof_.value * s;
}

boost::shared_ptr<InstrumentCorrection> EmpiricalOrthogonalFunction::clone() 
  const
{
  return boost::shared_ptr<InstrumentCorrection>
    (new EmpiricalOrthogonalFunction(coeff.value()(0), used_flag(0), 
				     eof_, order_,band_name, hdf_group,
				     sounding_number_, 
				     eof_depend_on_sounding_number_));
}

void EmpiricalOrthogonalFunction::apply_correction
(const SpectralDomain& Pixel_grid,
 const std::vector<int>& Pixel_list,
 SpectralRange& Radiance) const
{
  ArrayWithUnit<double, 1> eofv = eof_.convert(Radiance.units());
  if(Radiance.data_ad().number_variable() == 0) {
    for(int i = 0; i < Radiance.data().rows(); ++i) {
      range_check(Pixel_list[i], 0, eofv.value.rows());
      Radiance.data()(i) += coeff.value()(0) * eofv.value(Pixel_list[i]);
    }
  } else {
    for(int i = 0; i < Radiance.data_ad().rows(); ++i) {
      range_check(Pixel_list[i], 0, eofv.value.rows());
      Radiance.data_ad()(i) = Radiance.data_ad()(i) +
	coeff(0) * eofv.value(Pixel_list[i]);
    }
  }
}

void EmpiricalOrthogonalFunction::print(std::ostream& Os) const
{
  Os << "EmpiricalOrthogonalFunction:\n"
     << "  Band:                      " << band_name << "\n"
     << "  HDF group:                 " << hdf_group << "\n"
     << "  Order:                     " << order_ << "\n"
     << "  Sounding number:           " << sounding_number_ << "\n"
     << "  Depend on sounding number: " 
     << (eof_depend_on_sounding_number_ ? "yes" : "no") << "\n"
     << "  Scale by uncertainty:      "
     << (eof_scale_uncertainty_ ? "yes" : "no") << "\n"
     << "  Scale:                     " << coeff(0).value() << "\n"
     << "  Fit:                       " 
     << (used_flag(0) ? "True\n" : "False\n");
  if(eof_scale_uncertainty_)
    Os << "  Scale to stddev:           " << scale_to_stddev_ << "\n";
}
