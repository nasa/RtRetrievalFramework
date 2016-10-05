#ifndef LEVEL_1B_AVERAGE_H
#define LEVEL_1B_AVERAGE_H
#include "level_1b.h"
#include <boost/shared_ptr.hpp>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This reads averages a set of Level1b classes to get the various
  values. This is used for example on Gosat, where we average the S
  and P data.
*******************************************************************/
class Level1bAverage: public Level1b {
public:
  Level1bAverage(const std::vector<boost::shared_ptr<Level1b> >& Data)
    : l1b(Data) {}
  virtual ~Level1bAverage() {}
  virtual int number_spectrometer() const 
  {return l1b[0]->number_spectrometer();}
  virtual DoubleWithUnit relative_velocity(int Spec_index) const
  {return l1b[0]->relative_velocity(Spec_index);}
  virtual Time time(int Spec_index) const
  {return l1b[0]->time(Spec_index);}
  virtual int64_t sounding_id() const
  {return l1b[0]->sounding_id();}
  virtual int exposure_index() const
  {return l1b[0]->exposure_index();}

  virtual DoubleWithUnit latitude(int i) const;
  virtual DoubleWithUnit longitude(int i) const;
  virtual DoubleWithUnit sounding_zenith(int i) const;
  virtual DoubleWithUnit sounding_azimuth(int i) const;
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const;
  virtual DoubleWithUnit solar_zenith(int i) const;
  virtual DoubleWithUnit solar_azimuth(int i) const;
  virtual DoubleWithUnit altitude(int i) const;
  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const;

  virtual void print(std::ostream& Os) const;
  virtual SpectralRange radiance(int Spec_index) const;
private:
  std::vector<boost::shared_ptr<Level1b> > l1b;
};
}
#endif
