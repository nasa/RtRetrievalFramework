#include "level_1b_acos.h"
#include "fp_exception.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<NoiseModel> level_1b_acos_noise_model_get(const Level1bAcos& F)
{ return F.noise_model(); }
void level_1b_acos_noise_model_set(Level1bAcos& F, const boost::shared_ptr<NoiseModel>& Nm)
{ F.noise_model(Nm); }
REGISTER_LUA_DERIVED_CLASS(Level1bAcos, Level1b)
.def(luabind::constructor<std::string, 
			  const boost::shared_ptr<HdfSoundingId>&>())
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&,
			  const boost::shared_ptr<HdfSoundingId>&>())
.def("hdf_sounding_id", &Level1bAcos::hdf_sounding_id)
.def("land_fraction", &Level1bAcos::land_fraction)
.def("is_h_gain", &Level1bAcos::is_h_gain)
.def("is_m_gain", &Level1bAcos::is_m_gain)
.property("noise_model", 
	  &level_1b_acos_noise_model_get,
	  &level_1b_acos_noise_model_set)
REGISTER_LUA_END()
#endif

Level1bAcos::Level1bAcos(const std::string& Fname, 
			 const boost::shared_ptr<HdfSoundingId>& Sounding_id) 
: Level1bHdf(Fname, Sounding_id) 
{
  try {
    initialize();
  } catch(Exception& E) {
    E << " in the file " << file_name;
    throw E;
  }
}

Level1bAcos::Level1bAcos(const boost::shared_ptr<HdfFile>& Hfile, 
			 const boost::shared_ptr<HdfSoundingId>& Sounding_id)
: Level1bHdf(Hfile, Sounding_id) 
{
  try {
    initialize();
  } catch(Exception& E) {
    E << " in the file " << file_name;
    throw E;
  }
}

void Level1bAcos::initialize()
{
  firstIndex i1; secondIndex i2; thirdIndex(i3);
//-----------------------------------------------------------------------
/// Read all of the data we need.
//-----------------------------------------------------------------------

  int sindex = hdf_sounding_id_->frame_number();
  int sounding_number = hdf_sounding_id_->sounding_number();
  TinyVector<int, 3> asz = hfile->read_shape<3>("FootprintGeometry/footprint_altitude");
  TinyVector<int, 3> start(sindex, 0, sounding_number);
  TinyVector<int, 3> sz(1, asz[1], 1);
  ArrayWithUnit<double, 3> alt = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_altitude", units::m, start, sz);
  ArrayWithUnit<double, 3> lat = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_latitude", units::deg, start, sz);
  ArrayWithUnit<double, 3> lon = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_longitude", units::deg, start, sz);
  ArrayWithUnit<double, 3> saz = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_solar_azimuth", units::deg, start, sz);
  ArrayWithUnit<double, 3> szn = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_solar_zenith", units::deg, start, sz);
  ArrayWithUnit<double, 3> azm = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_azimuth", units::deg, start, sz);
  ArrayWithUnit<double, 3> zen = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_zenith", units::deg, start, sz);
  ArrayWithUnit<double, 1> relv = 
    hfile->read_field_with_unit<double, 1>("SpacecraftGeometry/relative_velocity", units::m / units::s, TinyVector<int, 1>(sindex), TinyVector<int, 1>(1));
  TinyVector<int, 4> sz2 = hfile->read_shape<4>("SoundingHeader/wavenumber_coefficients");
  ArrayWithUnit<double, 4> wn_coeffs = 
    hfile->read_field_with_unit<double, 4>("SoundingHeader/wavenumber_coefficients", units::inv_cm, TinyVector<int, 4>(sindex, 0, sounding_number, 0), TinyVector<int, 4>(1, sz2[1], 1, sz2[3]));

  TinyVector<int, 4> sz3 = hfile->read_shape<4>("FootprintGeometry/footprint_stokes_coefficients");
  Array<double, 4> st =
    hfile->read_field<double, 4>("FootprintGeometry/footprint_stokes_coefficients",  TinyVector<int, 4>(sindex, 0, sounding_number, 0), TinyVector<int, 4>(1, sz3[1], 1, sz3[3]));
  Array<double, 3> land_frac;
  // This field might be missing for older data, or for simulator data.
  try {
    land_frac.reference
      (hfile->read_field<double, 3>("FootprintGeometry/footprint_land_fraction", start, sz));
  } catch( Exception not_found_error ) {
    land_frac.resize(zen.value.shape());
    land_frac = 100.0;
  }
  char gain_code;
  // This field might be missing for older data, or for simulator data.
  try {
    Array<std::string, 2> ds_gain_codes = 
      hfile->read_field<std::string, 2>("SoundingHeader/gain_swir", 
			      TinyVector<int, 2>(sindex, sounding_number), 
			      TinyVector<int, 2>(1,1));
    gain_code = ds_gain_codes(0, 0)[0];
  } catch( Exception not_found_error ) {
    gain_code = 'H';
  }    
  is_h_gain_ = (gain_code == 'H');

  Array<double, 3> tm;
  try {
    // Try L1B 2.07.00 dataset name
    tm.reference
      (hfile->read_field<double, 3>("FootprintGeometry/footprint_time_tai93",
				   TinyVector<int, 3>(sindex, 0, 0),
				    TinyVector<int, 3>(1, 1, 1)));
  } catch( Exception not_found_error ) {
    tm.reference
      (hfile->read_field<double, 3>("FootprintGeometry/footprint_time",
				    TinyVector<int, 3>(sindex, 0, 0),
				    TinyVector<int, 3>(1, 1, 1)));
  }

  //-----------------------------------------------------------------------
  /// Finally, subset the data for the sounding and mode we are doing
  /// (P, S, or avg).
  //-----------------------------------------------------------------------

  Range ra(Range::all());
  altitude_.units = alt.units;
  latitude_.units = lat.units;
  longitude_.units = lon.units;
  solar_azimuth_.units = saz.units;
  solar_zenith_.units = szn.units;
  sounding_zenith_.units = zen.units;
  sounding_azimuth_.units = azm.units;
  spectral_coefficient_.units = wn_coeffs.units;
  
  altitude_.value.resize(alt.value.cols());
  latitude_.value.resize(altitude_.value.shape());
  longitude_.value.resize(altitude_.value.shape());
  solar_azimuth_.value.resize(altitude_.value.shape());
  solar_zenith_.value.resize(altitude_.value.shape());
  sounding_zenith_.value.resize(altitude_.value.shape());
  sounding_azimuth_.value.resize(altitude_.value.shape());
  spectral_coefficient_.value.resize(wn_coeffs.value.extent(1), wn_coeffs.value.extent(3));
  stokes_coef_.resize(altitude_.value.rows(), 4);
  land_fraction_.resize(altitude_.value.shape());

  altitude_.value = alt.value(0, ra, 0);
  latitude_.value = lat.value(0, ra, 0);
  longitude_.value = lon.value(0, ra, 0);
  solar_azimuth_.value = saz.value(0, ra, 0);
  solar_zenith_.value = szn.value(0, ra, 0);
  sounding_zenith_.value = zen.value(0, ra, 0);
  sounding_azimuth_.value = azm.value(0, ra, 0);
  spectral_coefficient_.value = wn_coeffs.value(0, ra, 0, ra);
  stokes_coef_ = st(0, ra, 0, ra);
  relative_velocity_ = relv(0);
  time_ = Time::time_pgs(tm(0,0,0)); // Use time of first band, first sounding.
  land_fraction_ = land_frac(0, ra, 0);
}

SpectralRange Level1bAcos::radiance_no_uncertainty(int Spec_index) const
{
  firstIndex i1; secondIndex i2;

//-----------------------------------------------------------------------
/// Radiance dataset names hardcoded for now.
//-----------------------------------------------------------------------

  std::string field;
  switch(Spec_index) {
  case 0:
    field = "SoundingSpectra/radiance_o2";
    break;
  case 1:
    field = "SoundingSpectra/radiance_weak_co2";
    break;
  case 2:
    field = "SoundingSpectra/radiance_strong_co2";
    break;
  default:
    throw Exception("Unrecognized Spec_Index");
  }
  TinyVector<int, 3> sz = hfile->read_shape<3>(field);
  ArrayWithUnit<double, 3> rad = hfile->read_field_with_unit<double, 3>
    (field, Unit("W / cm^2 / sr / cm^-1"),
     TinyVector<int, 3>(hdf_sounding_id_->frame_number(), 
			hdf_sounding_id_->sounding_number(), 0),
     TinyVector<int, 3>(1, 1, sz[2]));
  return SpectralRange(rad.value(0, 0, Range::all()), rad.units);
}
