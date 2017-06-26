#include "level_1b_uq.h"
#include "fp_exception.h"
#include "old_constant.h"
#include <boost/regex.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<NoiseModel> level_1b_uq_noise_model_get(const Level1bUq& F)
{
    return F.noise_model();
}
void level_1b_uq_noise_model_set(Level1bUq& F, const boost::shared_ptr<NoiseModel>& Nm)
{
    F.noise_model(Nm);
}
REGISTER_LUA_DERIVED_CLASS(Level1bUq, Level1b)
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&,
     const boost::shared_ptr<HdfSoundingId>&>())
.property("noise_model",
          &level_1b_uq_noise_model_get,
          &level_1b_uq_noise_model_set)
.def("has_solar_relative_velocity", &Level1bUq::has_solar_relative_velocity)
.def("land_fraction", &Level1bUq::land_fraction)
.def("solar_distance", &Level1bUq::solar_distance)
.def("solar_velocity", &Level1bUq::solar_velocity)
.def("acquisition_mode", &Level1bUq::acquisition_mode)
.def("has_spike_eof", &Level1bUq::has_spike_eof)
.def("spike_eof", &Level1bUq::spike_eof)
REGISTER_LUA_END()
#endif

Level1bUq::Level1bUq(const boost::shared_ptr<HdfFile>& Hfile,
                     const boost::shared_ptr<HdfSoundingId>& Sounding_id)
{
    try {
        file_name = Hfile->file_name();
        hfile = Hfile;
        hdf_sounding_id_ = Sounding_id;

        initialize();
    } catch(Exception& E) {
        E << " in the file " << file_name;
        throw E;
    }
}

void Level1bUq::initialize()
{
    Range ra(Range::all());

    //-----------------------------------------------------------------------
    /// Read all of the data we need.
    //-----------------------------------------------------------------------

    // All soundings and spectrometers in a UQ file have most of the same values
    // So most datasets are either just arrays or scalars. But need to resize
    // the internal data by spectrometers, so we get the number of spectrometers
    // from a dataset that has this information
    TinyVector<int, 2> disp_shape = hfile->read_shape<2>("InstrumentHeader/dispersion_coef_samp");
    int num_spec = disp_shape(0);

    // Load geometry datasets
    ArrayWithUnit<double, 1> alt = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_altitude", units::m);
    ArrayWithUnit<double, 1> lat = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_latitude", units::deg);
    ArrayWithUnit<double, 1> lon = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_longitude", units::deg);
    ArrayWithUnit<double, 1> saz = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_solar_azimuth", units::deg);
    ArrayWithUnit<double, 1> szn = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_solar_zenith", units::deg);
    ArrayWithUnit<double, 1> azm = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_azimuth", units::deg);
    ArrayWithUnit<double, 1> zen = 
        hfile->read_field_with_unit<double, 1>("FootprintGeometry/footprint_zenith", units::deg);

    altitude_.value.resize(num_spec);
    latitude_.value.resize(altitude_.value.shape());
    longitude_.value.resize(altitude_.value.shape());
    solar_azimuth_.value.resize(altitude_.value.shape());
    solar_zenith_.value.resize(altitude_.value.shape());
    sounding_zenith_.value.resize(altitude_.value.shape());
    sounding_azimuth_.value.resize(altitude_.value.shape());

    altitude_.value = alt.value(0);
    latitude_.value = lat.value(0);
    longitude_.value = lon.value(0);
    solar_azimuth_.value = saz.value(0);
    solar_zenith_.value = szn.value(0);
    sounding_zenith_.value = zen.value(0);
    sounding_azimuth_.value = azm.value(0);

    altitude_.units = alt.units;
    latitude_.units = lat.units;
    longitude_.units = lon.units;
    solar_azimuth_.units = saz.units;
    solar_zenith_.units = szn.units;
    sounding_zenith_.units = zen.units;
    sounding_azimuth_.units = azm.units;

    // Load time dataset
    Array<double, 1> tm =
        hfile->read_field<double, 1>("FootprintGeometry/footprint_time_tai93");
    time_ = Time::time_pgs(tm(0));

    // Load spacecraft relative velocity dataset
    if (hfile->has_object("SoundingGeometry/sounding_relative_velocity")) {
        ArrayWithUnit<double, 1> relv = hfile->read_field_with_unit<double, 1>("SoundingGeometry/sounding_relative_velocity", units::m / units::s);
        relative_velocity_ = relv(0);
    } else {
        relative_velocity_ = 0;
    }

    // Load dispersion coefficients
    ArrayWithUnit<double, 2> wl_coeffs;
    wl_coeffs = hfile->read_field_with_unit<double, 2>("/InstrumentHeader/dispersion_coef_samp", units::micron);
    spectral_coefficient_.value.resize(wl_coeffs.value.extent(0), wl_coeffs.value.extent(1));
    spectral_coefficient_.value = wl_coeffs.value(ra, ra);
    spectral_coefficient_.units = wl_coeffs.units;

    // Read the stokes coefficients
    std::string stokes_dataset = "FootprintGeometry/footprint_stokes_coefficients";
    TinyVector<int, 1> stokes_shape = hfile->read_shape<1>(stokes_dataset);
    stokes_coef_.resize(altitude_.value.rows(), stokes_shape(0));

    Array<double, 1> st = hfile->read_field<double, 1>(stokes_dataset);
    for(int spec_idx = 0; spec_idx < stokes_coef_.rows(); spec_idx++) {
        stokes_coef_(spec_idx, ra) = st(ra);
    }

    // Read land fraction dataset
    Array<double, 1> lf = hfile->read_field<double, 1>("SoundingGeometry/sounding_land_fraction");
    land_fraction_ = lf(0);

    // Read solar relative velocity and solar distance datasets
    if(hfile->has_object("/SoundingGeometry/sounding_solar_relative_velocity")) {
        has_solar_relative_velocity_ = true;
        ArrayWithUnit<double, 1> sdist =
            hfile->read_field_with_unit<double, 1>("/SoundingGeometry/sounding_solar_distance", units::m);
        ArrayWithUnit<double, 1> svel =
            hfile->read_field_with_unit<double, 1>("/SoundingGeometry/sounding_solar_relative_velocity", units::m / units::s);
        solar_distance_ = sdist(0);
        solar_velocity_ = svel(0);
    } else {
        has_solar_relative_velocity_ = false;
        solar_distance_ = 0.0;
        solar_velocity_ = 0.0;
    }

    // Read acquisition mode dataset
    if(hfile->has_object("/Metadata/AcquisitionMode")) {
        Array<std::string, 1> d =
            hfile->read_field<std::string, 1>("/Metadata/AcquisitionMode");
        // Sometimes test data add "Sample" to the front. Remove this if
        // found, just so we don't need to have special handing for this.
        acquisition_mode_ = boost::regex_replace(d(0),
                            boost::regex("Sample\\s*"), "");
    } else {
        // Default to "Nadir". Real data will always have the metadata,
        // but the OCO2 simulator data doesn't contain this field.
        acquisition_mode_ = "Nadir";
    }
}

void Level1bUq::set_radiance(int Spec_index, boost::shared_ptr<SpectralRange>& Rad)
{
  if ((int) radiance.size() < Spec_index + 1) {
        radiance.resize(Spec_index + 1);
    }

    radiance[Spec_index] = Rad;
}

SpectralRange Level1bUq::radiance_no_uncertainty(int Spec_index) const
{
  if((int) radiance.size() < Spec_index + 1) {
        Exception err_msg;
        err_msg << "Radiance for spectrometer: " << Spec_index << " not yet set.";
        throw err_msg;
    }

    return *radiance[Spec_index];
}

