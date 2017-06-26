#ifndef LEVEL_1B_HERITAGE_H
#define LEVEL_1B_HERITAGE_H
#include "level_1b.h"
#include "heritage_file.h"
#include "fp_exception.h"
#include "noise_model.h"

namespace FullPhysics {
/****************************************************************//**
  This reads a Level 1B file that is in the heritage text format.
*******************************************************************/
class Level1bHeritage: public Level1b {
public:
  Level1bHeritage(const std::string& Sounding_info_file, 
		  const std::string& Spectrum_file,
		  const boost::shared_ptr<NoiseModel>& Nm);
  virtual ~Level1bHeritage() {}
  virtual SpectralRange radiance(int Spec_index) const;
  virtual int number_spectrometer() const
  {
    return number_spectrometer_;
  }
  virtual DoubleWithUnit latitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_latitude")[i], units::deg);
  }
  virtual DoubleWithUnit longitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_longitude")[i], units::deg);
  }
  virtual DoubleWithUnit solar_zenith(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_solar_zenith")[i], units::deg);
  }
  virtual DoubleWithUnit sounding_zenith(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_zenith")[i], units::deg);
  }
  virtual DoubleWithUnit sounding_azimuth(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_azimuth")[i], units::deg);
  }
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const
  {
    range_check(i, 0, number_spectrometer());
    blitz::Array<double, 1> res(4);
    std::vector<double> t(hmeta.value<std::vector<double> >("stokes_coefficients"));
    res = t[i * 4], t[i * 4 + 1], t[i * 4 + 2], t[i * 4 + 3];
    return res;
  }
  virtual DoubleWithUnit solar_azimuth(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_solar_azimuth")[i], units::deg);
  }
  virtual DoubleWithUnit altitude(int i) const
  {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(hmeta.value<std::vector<double> >("sounding_altitude")[i], units::m);
  }

  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const
  { throw Exception("Level1BHeritage does not support returning spectral_coefficient data"); }

  virtual DoubleWithUnit relative_velocity(int Spec_index) const
  { return DoubleWithUnit(hmeta.value<double>("relative_velocity"), units::m / units::s); }
  virtual Time time(int Spec_index) const
  { return hmeta.value<Time>("frame_time_stamp"); }
  virtual int64_t sounding_id() const
  { return hmeta.value<int64_t>("sounding_id"); }
  virtual int exposure_index() const { return 1; }
  virtual void print(std::ostream& Os) const
  {
    Os << "Level 1B Heritage file:\n"
       << "  File name: " << hmeta.file_name() << "\n";
  }
public:
  HeritageFile hmeta;
  HeritageFile hspec;
  boost::shared_ptr<NoiseModel> noise_model_;
  int number_spectrometer_;
  std::vector<blitz::Range> r;
};
}
#endif
