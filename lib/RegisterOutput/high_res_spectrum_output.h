#ifndef HIGH_RES_SPECTRUM_OUTPUT_H
#define HIGH_RES_SPECTRUM_OUTPUT_H
#include "register_output_base.h"
#include "named_spectrum.h"
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
 This class will recieve observer notifications from any class
 that pushes out NamedSpectrum and write them to the output file
 in individual datasets for the spectral domain and spectral range
 parts.

 This is intended to be a way to write high resolution spectrum
 from the forward model without cluttering it with accessor
 functions and without knowing the exact steps that will be 
 performed in it.
*******************************************************************/
class HighResSpectrumOutput : public RegisterOutputBase,
  public Observer<boost::shared_ptr<NamedSpectrum> >,
  public Observer<std::vector<boost::shared_ptr<NamedSpectrum> > > {
public:
  HighResSpectrumOutput() {}
  virtual ~HighResSpectrumOutput() {}

  /// Saves a pointer to all interested output files
  /// The datasets will be actually registered as they are
  /// send via notifications
  virtual void register_output(const boost::shared_ptr<Output>& out) const { 
    output_files.push_back(out);
  }
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const { }

  const blitz::Array<double, 1> saved_spectral_domain(const std::string& spectra_name);
  const blitz::Array<double, 1> saved_spectral_range(const std::string& spectra_name);

  virtual void notify_update(const boost::shared_ptr<NamedSpectrum>& Obs);
  virtual void notify_update(const std::vector<boost::shared_ptr<NamedSpectrum> >& Obs);

  virtual std::string desc() const { return "HighResSpectrumOutput"; }
private:
  void register_named_spectrum(const boost::shared_ptr<FullPhysics::NamedSpectrum>& named_spec,
                               const std::string& domain_prefix, const std::string& range_prefix);

  std::map<std::string, std::vector<boost::shared_ptr<Spectrum> > > saved_spectra;
  mutable std::vector<boost::shared_ptr<Output> > output_files;
};
}
#endif
