#ifndef LEVEL_1B_SCALE_RADIANCE_H
#define LEVEL_1B_SCALE_RADIANCE_H
#include "level_1b.h"
#include <boost/shared_ptr.hpp>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  Scales the measured radiance of another Level1b class for each
  spectral band.
*******************************************************************/

class Level1bScaleRadiance: public Level1b {
public:
    Level1bScaleRadiance(const boost::shared_ptr<Level1b>& L1b, const blitz::Array<double, 1>& Scaling)
    : l1b(L1b), scaling(Scaling) {}

    virtual ~Level1bScaleRadiance() {}

    virtual int number_spectrometer() const 
    { return l1b->number_spectrometer(); }

    virtual DoubleWithUnit relative_velocity(int Spec_index) const
    { return l1b->relative_velocity(Spec_index); }

    virtual Time time(int Spec_index) const
    { return l1b->time(Spec_index); }

    virtual int64_t sounding_id() const
    { return l1b->sounding_id(); }

    virtual int exposure_index() const
    { return l1b->exposure_index(); }

    virtual DoubleWithUnit latitude(int i) const
    { return l1b->latitude(i); }

    virtual DoubleWithUnit longitude(int i) const
    { return l1b->longitude(i); }

    virtual DoubleWithUnit sounding_zenith(int i) const
    { return l1b->sounding_zenith(i); }

    virtual DoubleWithUnit sounding_azimuth(int i) const
    { return l1b->sounding_azimuth(i); }

    virtual blitz::Array<double, 1> stokes_coefficient(int i) const
    { return l1b->stokes_coefficient(i); }

    virtual DoubleWithUnit solar_zenith(int i) const
    { return l1b->solar_zenith(i); }

    virtual DoubleWithUnit solar_azimuth(int i) const
    { return l1b->solar_azimuth(i); }

    virtual DoubleWithUnit altitude(int i) const
    { return l1b->altitude(i); }

    virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const
    { return l1b->spectral_coefficient(Spec_index); }

    virtual SpectralRange radiance(int Spec_index) const;

    virtual void print(std::ostream& Os) const;
private:
    boost::shared_ptr<Level1b> l1b;
    blitz::Array<double, 1> scaling;
};
}
#endif
