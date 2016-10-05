// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_instrument.h"
#include "sub_state_vector_array.h"
%}
%base_import(instrument)
%import "instrument_correction.i"
%import "ils.i"
%import "instrument_correction.i"
%import "spectrum.i"
%import "state_vector.i"
%fp_shared_ptr(FullPhysics::IlsInstrument);

namespace FullPhysics {
class IlsInstrument : public Instrument {
public:
  IlsInstrument(const std::vector<boost::shared_ptr<Ils> >& Ils_list,
		const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >&
		Instrument_correction = 
		std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >());
  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const;
  virtual std::string band_name(int Spec_index) const;
  virtual std::string hdf_band_name(int Spec_index) const;
  virtual void notify_update(const StateVector& Sv);
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const Ils& D);
  virtual void notify_update(const InstrumentCorrection& C);
  virtual boost::shared_ptr<Instrument> clone() const;
  boost::shared_ptr<Ils> ils(int Spec_index) const;
  const std::vector<boost::shared_ptr<InstrumentCorrection> >&
  instrument_correction(int Spec_index) const;
};
}
