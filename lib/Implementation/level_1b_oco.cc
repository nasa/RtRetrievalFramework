#include "level_1b_oco.h"
#include "fp_exception.h"
#include "old_constant.h"
#include <boost/regex.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<NoiseModel> level_1b_oco_noise_model_get(const Level1bOco& F)
{ return F.noise_model(); }
void level_1b_oco_noise_model_set(Level1bOco& F, const boost::shared_ptr<NoiseModel>& Nm)
{ F.noise_model(Nm); }
REGISTER_LUA_DERIVED_CLASS(Level1bOco, Level1b)
.def(luabind::constructor<const std::string&, 
                          const boost::shared_ptr<HdfSoundingId>&>())
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&, 
                          const boost::shared_ptr<HdfSoundingId>&>())
.property("noise_model", 
          &level_1b_oco_noise_model_get,
          &level_1b_oco_noise_model_set)
.def("has_solar_relative_velocity", &Level1bOco::has_solar_relative_velocity)
.def("land_fraction", &Level1bOco::land_fraction)
.def("solar_distance", &Level1bOco::solar_distance)
.def("solar_velocity", &Level1bOco::solar_velocity)
.def("acquisition_mode", &Level1bOco::acquisition_mode)
.def("has_spike_eof", &Level1bOco::has_spike_eof)
.def("spike_eof", &Level1bOco::spike_eof)
REGISTER_LUA_END()
#endif

Level1bOco::Level1bOco(const std::string& Fname, 
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

Level1bOco::Level1bOco(const boost::shared_ptr<HdfFile>& Hfile, 
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

void Level1bOco::initialize()
{
  firstIndex i1; secondIndex i2; thirdIndex(i3);
  Range ra(Range::all());

//-----------------------------------------------------------------------
/// Read all of the data we need.
//-----------------------------------------------------------------------

  int frame_index = hdf_sounding_id_->frame_number();
  int sounding_index = hdf_sounding_id_->sounding_number();

  TinyVector<int, 3> alt_shape = hfile->read_shape<3>("FootprintGeometry/footprint_altitude");
  TinyVector<int, 3> fg_start, fg_size;
  fg_start = 
    frame_index, sounding_index, 0;
  fg_size = 1, 1, alt_shape(2);

  ArrayWithUnit<double, 3> alt = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_altitude", units::m, fg_start, fg_size);
  ArrayWithUnit<double, 3> lat = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_latitude", units::deg, fg_start, fg_size);
  ArrayWithUnit<double, 3> lon = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_longitude", units::deg, fg_start, fg_size);
  ArrayWithUnit<double, 3> saz = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_solar_azimuth", units::deg, fg_start, fg_size);
  ArrayWithUnit<double, 3> szn = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_solar_zenith", units::deg, fg_start, fg_size);
  ArrayWithUnit<double, 3> azm = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_azimuth", units::deg, fg_start, fg_size);
  ArrayWithUnit<double, 3> zen = 
    hfile->read_field_with_unit<double, 3>("FootprintGeometry/footprint_zenith", units::deg, fg_start, fg_size);

  Array<double, 3> tm =
    hfile->read_field<double, 3>("FootprintGeometry/footprint_time_tai93", fg_start, fg_size);

  TinyVector<int, 1> relv_start, relv_size;
  relv_start = frame_index;
  relv_size = 1;

  if (hfile->has_object("SoundingGeometry/sounding_relative_velocity")) {
    TinyVector<int,2> vstart, vsize;
    vstart = frame_index, sounding_index;
    vsize = 1, 1;
    ArrayWithUnit<double, 2> relv = hfile->read_field_with_unit<double, 2>("SoundingGeometry/sounding_relative_velocity", units::m / units::s, vstart, vsize);
    relative_velocity_ = relv(0, 0);
  } else {
    ArrayWithUnit<double, 1> relv = hfile->read_field_with_unit<double, 1>("FrameGeometry/relative_velocity", units::m / units::s, relv_start, relv_size);
    relative_velocity_ = relv(0);
  }

  ArrayWithUnit<double, 3> wl_coeffs;
  // This field has an old and new name. Try new name first, and then
  // old one if we don't find it
  try {
    wl_coeffs = hfile->read_field_with_unit<double, 3>
      ("/InstrumentHeader/dispersion_coef_samp", units::micron);
  } catch( Exception not_found_error ) {
    wl_coeffs = hfile->read_field_with_unit<double, 3>
      ("Metadata/DispersionCoefSamp", units::micron);
  }

  //-----------------------------------------------------------------------
  /// Finally, subset the data for the sounding and mode we are doing
  //-----------------------------------------------------------------------

  altitude_.units = alt.units;
  latitude_.units = lat.units;
  longitude_.units = lon.units;
  solar_azimuth_.units = saz.units;
  solar_zenith_.units = szn.units;
  sounding_zenith_.units = zen.units;
  sounding_azimuth_.units = azm.units;
  spectral_coefficient_.units = wl_coeffs.units;

  altitude_.value.resize(alt.value.depth());
  latitude_.value.resize(altitude_.value.shape());
  longitude_.value.resize(altitude_.value.shape());
  solar_azimuth_.value.resize(altitude_.value.shape());
  solar_zenith_.value.resize(altitude_.value.shape());
  sounding_zenith_.value.resize(altitude_.value.shape());
  sounding_azimuth_.value.resize(altitude_.value.shape());
  spectral_coefficient_.value.resize(wl_coeffs.value.extent(0), wl_coeffs.value.extent(2));

  // We hyperslabed into the dataset above, so we should have arrays
  // with shape 1 x 1 x 3
  altitude_.value = alt.value(0, 0, ra);
  latitude_.value = lat.value(0, 0, ra);
  longitude_.value = lon.value(0, 0, ra);
  solar_azimuth_.value = saz.value(0, 0, ra);
  solar_zenith_.value = szn.value(0, 0, ra);
  sounding_zenith_.value = zen.value(0, 0, ra);
  sounding_azimuth_.value = azm.value(0, 0, ra);
  time_ = Time::time_pgs(tm(0,0,0)); // Use time of first band, 

  spectral_coefficient_.value = wl_coeffs.value(ra, sounding_index, ra);

  // Read the stokes coefficients from the appropriate dataset
  stokes_coef_.resize(altitude_.value.rows(), 4);
  
  // Try reading one of these datasets
  std::string stokes_dataset = "FootprintGeometry/footprint_stokes_coefficients";
  std::string pol_ang_dataset = "FootprintGeometry/footprint_polarization_angle";

  // Determine if we can read from the stokes dataset
  if (hfile->has_object(stokes_dataset)) {
    TinyVector<int,4> stokes_shape = hfile->read_shape<4>(stokes_dataset);
    TinyVector<int,4> stokes_start, stokes_size;
    stokes_start = frame_index, sounding_index, 0, 0;
    stokes_size = 1, 1, stokes_shape(2), stokes_shape(3);

    Array<double, 4> st =
      hfile->read_field<double, 4>(stokes_dataset, stokes_start, stokes_size);
    stokes_coef_ = st(0, 0, ra, ra);

  } else if(hfile->has_object(pol_ang_dataset)) {
    // OCO1 did not supply stokes coefficients, instead had
    // polarization angle which can be converted to
    // stokes coefficients

    Array<double, 3> pol_ang =
      hfile->read_field<double, 3>(pol_ang_dataset);
    
    // Note that the 1/2 here is due to the fact OCO does not have a polarization
    // filter
    stokes_coef_(ra, 0) = 0.5;
    for(int bidx = 0; bidx < stokes_coef_.extent(firstDim); bidx++) {
      stokes_coef_(bidx,1) = 0.5*cos(2.0*pol_ang(frame_index, sounding_index, bidx)*FullPhysics::conversion(units::deg, units::rad));
      stokes_coef_(bidx,2) = 0.5*sin(2.0*pol_ang(frame_index, sounding_index, bidx)*FullPhysics::conversion(units::deg, units::rad));
    }
    stokes_coef_(ra, 3) = 0.0;

  } else {
    // Polarization angle non existant, for example in uplooking data,
    // Set to default value below: polarization angle == 90.0
    stokes_coef_(ra, 0) = 0.5;
    stokes_coef_(ra, 1) = -0.5;
    stokes_coef_(ra, Range(2, toEnd)) = 0.0;
  }
  if(hfile->has_object("SoundingGeometry/sounding_land_fraction")) {
    TinyVector<int,2> land_start, land_size;
    land_start = frame_index, sounding_index;
    land_size = 1, 1;
    Array<double, 2> lf = hfile->read_field<double, 2>("SoundingGeometry/sounding_land_fraction", land_start, land_size);
    land_fraction_ = lf(0,0);
  } else {
    // For backward compatibility, if we don't have the land fraction field
    // use the land_water_indicator and translate to an fraction
    TinyVector<int,2> land_start, land_size;
    land_start = frame_index, sounding_index;
    land_size = 1, 1;
    Array<int, 2> lind = hfile->read_field<int, 2>("SoundingGeometry/sounding_land_water_indicator", land_start, land_size);
    int lin = lind(0,0);
    switch(lin) {
    case 0:
    case 3:
      land_fraction_ = 100.0;
      break;
    case 1:
    case 2:
      land_fraction_ = 0.0;
      break;
    default:
      throw Exception("Invalid sounding_land_water_indicator value");
    }
  }
  if(hfile->has_object("/SoundingGeometry/sounding_solar_relative_velocity")) {
    has_solar_relative_velocity_ = true;
    TinyVector<int,2> solar_start, solar_size;
    solar_start = frame_index, sounding_index;
    solar_size = 1, 1;
    ArrayWithUnit<double, 2> sdist =
      hfile->read_field_with_unit<double, 2>("/SoundingGeometry/sounding_solar_distance", units::m, solar_start, solar_size);
    ArrayWithUnit<double, 2> svel =
      hfile->read_field_with_unit<double, 2>("/SoundingGeometry/sounding_solar_relative_velocity", units::m / units::s, solar_start, solar_size);
    solar_distance_ = sdist(0,0);
    solar_velocity_ = svel(0,0);
  } else {
    has_solar_relative_velocity_ = false;
    solar_distance_ = 0.0;
    solar_velocity_ = 0.0;
  }
  if(hfile->has_object("/Metadata/AcquisitionMode")) {
    Array<std::string, 1> d = 
      hfile->read_field<std::string, 1>("/Metadata/AcquisitionMode");
    // Sometimes test data add "Sample" to the front. Remove this if
    // found, just so we don't need to have special handing for this.
    acquisition_mode_ = boost::regex_replace(d(0),
					     boost::regex("Sample\\s*"),"");
  } else {
    // Default to "Nadir". Real data will always have the metadata,
    // but the OCO2 simulator data doesn't contain this field.
    acquisition_mode_ = "Nadir";
  }
}

SpectralRange Level1bOco::radiance_no_uncertainty(int Spec_index) const
{
  firstIndex i1; secondIndex i2;

//-----------------------------------------------------------------------
/// Radiance dataset names hardcoded for now.
//-----------------------------------------------------------------------

  std::string field;
  switch(Spec_index) {
  case 0:
    field = "SoundingMeasurements/radiance_o2";
    break;
  case 1:
    field = "SoundingMeasurements/radiance_weak_co2";
    break;
  case 2:
    field = "SoundingMeasurements/radiance_strong_co2";
    break;
  default:
    throw Exception("Unrecognized Spec_Index");
  }

  // Use start and size to keep from having to read entire radiance
  // array each time, OCO datasets get rather large
  TinyVector<int, 3> rad_shape = hfile->read_shape<3>(field);
  TinyVector<int, 3> rad_start, rad_size;
  rad_start = 
    hdf_sounding_id_->frame_number(), hdf_sounding_id_->sounding_number(), 0;
  rad_size = 1, 1, rad_shape(2);

  ArrayWithUnit<double, 3> rad = 
    hfile->read_field_with_unit<double, 3>(field, 
                                           Unit("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}"),
                                           rad_start, rad_size);
  
  return SpectralRange(rad.value(0, 0, Range::all()), rad.units);
}

//-----------------------------------------------------------------------
/// True if we have the SpikeEOF data available.
//-----------------------------------------------------------------------

bool Level1bOco::has_spike_eof(int Spec_index) const
{
//-----------------------------------------------------------------------
/// Dataset names hardcoded for now.
//-----------------------------------------------------------------------

  std::string field;
  switch(Spec_index) {
  case 0:
    field = "SpikeEOF/spike_eof_weighted_residual_o2";
    break;
  case 1:
    field = "SpikeEOF/spike_eof_weighted_residual_weak_co2";
    break;
  case 2:
    field = "SpikeEOF/spike_eof_weighted_residual_strong_co2";
    break;
  default:
    throw Exception("Unrecognized Spec_Index");
  }
  return hfile->has_object(field);
}

//-----------------------------------------------------------------------
/// Return spike eof data, used to perform SAA detection.
//-----------------------------------------------------------------------

blitz::Array<double, 1> Level1bOco::spike_eof(int Spec_index) const
{
  if(!has_spike_eof(Spec_index))
    throw Exception("Spike EOF data is not found in the Level 1b file");

//-----------------------------------------------------------------------
/// Dataset names hardcoded for now.
//-----------------------------------------------------------------------

  std::string field;
  switch(Spec_index) {
  case 0:
    field = "SpikeEOF/spike_eof_weighted_residual_o2";
    break;
  case 1:
    field = "SpikeEOF/spike_eof_weighted_residual_weak_co2";
    break;
  case 2:
    field = "SpikeEOF/spike_eof_weighted_residual_strong_co2";
    break;
  default:
    throw Exception("Unrecognized Spec_Index");
  }

  // Use start and size to keep from having to read entire
  // array each time, OCO datasets get rather large
  TinyVector<int, 3> f_shape = hfile->read_shape<3>(field);
  TinyVector<int, 3> f_start, f_size;
  f_start = 
    hdf_sounding_id_->frame_number(), hdf_sounding_id_->sounding_number(), 0;
  f_size = 1, 1, f_shape(2);

  blitz::Array<double, 3> dat = 
    hfile->read_field<double, 3>(field, f_start, f_size);
  
  return dat(0, 0, Range::all());
}
