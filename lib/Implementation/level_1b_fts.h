#ifndef LEVEL_1B_FTS_H
#define LEVEL_1B_FTS_H
#include "level_1b.h"
#include "array_with_unit.h"
#include "spectral_bound.h"
#include "fts_run_log.h"
#include "fp_exception.h"
#include "logger.h"
#include <boost/lexical_cast.hpp>
#include <vector>
#include <fstream>

namespace FullPhysics {
/****************************************************************//**
  This is the Level 1B data for a FTS.
*******************************************************************/
class Level1bFts: public Level1b {
public:
  Level1bFts(const std::string& Runlog,
             const std::vector<std::string>& Spectra,
             const ArrayWithUnit<double, 2>& Spectral_range);

  Level1bFts(const std::string& Runlog,
             const std::vector<std::string>& Spectra,
             const SpectralBound& Spec_bound);

  Level1bFts(const HdfFile& Hfile,
             const std::vector<std::string>& Band_names,
             const std::string& Radiance_dataset = "/SpectralParameters/modeled_radiance");

  virtual ~Level1bFts() {}
  virtual SpectralRange radiance(int Spec_index) const 
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return radiances_[Spec_index];
  }

  virtual int number_spectrometer() const {
    return (int) run_log_.size();
  }
  virtual DoubleWithUnit latitude(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].latitude, units::deg);
  }
  virtual DoubleWithUnit longitude(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].longitude, units::deg);
  }
  virtual DoubleWithUnit solar_zenith(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].solar_zenith, units::deg);
  }
  virtual DoubleWithUnit sounding_zenith(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].solar_zenith + run_log_[i].zenith_offset, units::deg);
  }
  virtual DoubleWithUnit sounding_azimuth(int i) const {
    range_check(i, 0, number_spectrometer());
    // Not sure where azimuth comes from.
    return DoubleWithUnit(0.0, units::deg);
  }
  virtual blitz::Array<double, 1> stokes_coefficient(int i) const {
    range_check(i, 0, number_spectrometer());
    // Not sure where stokes_coefficients comes from
    //set_stokes_coefs_from_pol_ang i.e. get from instrument
    //FIXME
    blitz::Array<double, 1> dummy(4);
    dummy = 0;
    return dummy;
    }
  virtual DoubleWithUnit solar_azimuth(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].solar_azimuth, units::deg);
  }
  virtual DoubleWithUnit altitude(int i) const {
    range_check(i, 0, number_spectrometer());
    return DoubleWithUnit(run_log_[i].altitude, units::km).convert(units::m);
  }
  virtual DoubleWithUnit relative_velocity(int Spec_index) const { 
    // FTS doesn't move relative to ground
    return DoubleWithUnit(0.0, units::m / units::s); 
  }

  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const
  {
    blitz::Array<double, 1> params(2);
    params = frequency_start(Spec_index), frequency_spacing(Spec_index);

    ArrayWithUnit<double, 1> res;
    res.value.reference(params);
    res.units = units::inv_cm;

    return res;
  }

  virtual Time time(int Spec_index) const { 
    return run_log_[Spec_index].time;
  }
  virtual int64_t sounding_id() const { 
    // Protect code for dying just because spectum_index in the
    // runlog is not always a proper integer
    try {
      return boost::lexical_cast<int64_t>(run_log_[0].spectrum_index);
    } catch(const boost::bad_lexical_cast& e) {
      Logger::warning() << "Could not parse spectrum index value: "
                        << "\""
                        << run_log_[0].spectrum_index 
                        << "\""
                        << " from run log file for determining sounding_id.\n";
      return boost::lexical_cast<int64_t>(0);
    }
  }
  virtual int exposure_index() const {return 1; }
  virtual void print(std::ostream& Os) const {
    Os << "Level1bFts";
  }

  /// Access records read from FTS runlog file
  const FtsRunLogRecord& run_log(int Spec_index) const {
    return run_log_[Spec_index];
  }

  const std::vector<FtsRunLogRecord>& run_log() const {
    return run_log_;
  }

  const double frequency_start(int spec_idx) const { return freq_start(spec_idx); }
  const double frequency_end(int spec_idx) const { return freq_end(spec_idx); };
  const double frequency_spacing(int spec_idx) const { return freq_spacing(spec_idx); }

private:
  void initialize(const std::string& Runlog,
                  const std::vector<std::string>& Spectra,
                  const ArrayWithUnit<double, 2>& Spectral_range);

  std::vector<FtsRunLogRecord> run_log_;
  std::vector<SpectralRange> radiances_;

  blitz::Array<double, 1> freq_start;
  blitz::Array<double, 1> freq_end;
  blitz::Array<double, 1> freq_spacing;
};
}
#endif
