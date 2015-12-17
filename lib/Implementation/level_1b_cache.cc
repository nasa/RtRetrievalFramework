#include "level_1b_cache.h"
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

Level1bCache::Level1bCache(const Level1b& L1_in)
{
  sid = L1_in.sounding_id();
  expind = L1_in.exposure_index();
  for(int i = 0; i < L1_in.number_spectrometer(); ++i) {
    lat.push_back(L1_in.latitude(i));
    lon.push_back(L1_in.longitude(i));
    szen.push_back(L1_in.sounding_zenith(i));
    sazm.push_back(L1_in.sounding_azimuth(i));
    solzen.push_back(L1_in.solar_zenith(i));
    solazm.push_back(L1_in.solar_azimuth(i));
    alt.push_back(L1_in.altitude(i));
    rvel.push_back(L1_in.relative_velocity(i));
    stk_coeff.push_back(L1_in.stokes_coefficient(i).copy());
    spec_coeff.push_back(L1_in.spectral_coefficient(i));
    tm.push_back(L1_in.time(i));
    rad.push_back(L1_in.radiance(i).clone());
  }
}
