// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "instrument.h"
%}

%import "observer.i"
%import "spectrum.i"
%import "spectral_domain.i"
%import "double_with_unit.i"
%base_import(state_vector)

%fp_shared_ptr(FullPhysics::Instrument);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Instrument>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Instrument>)

namespace FullPhysics {
  class Instrument;
}

namespace FullPhysics {
%template(ObservableInstrument) FullPhysics::Observable<Instrument>;
%template(ObserverInstrument) FullPhysics::Observer<Instrument>;
class Instrument : virtual public StateVectorObserver, 
		   public Observable<Instrument> {
public:
  virtual ~Instrument();
  std::string print_to_string() const;
  virtual void add_observer(Observer<Instrument>& Obs);
  virtual void remove_observer(Observer<Instrument>& Obs);
  virtual boost::shared_ptr<Instrument> clone() const = 0;
  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const = 0;
  %python_attribute(number_spectrometer, virtual int);
  virtual SpectralDomain pixel_spectral_domain(int Spec_index) const = 0;
  virtual std::string band_name(int Spec_index) const = 0;
  virtual std::string hdf_band_name(int Spec_index) const;
  virtual DoubleWithUnit ils_half_width(int Spec_index) const;
  virtual void ils_half_width(int Spec_index, DoubleWithUnit& half_width);
};
}
