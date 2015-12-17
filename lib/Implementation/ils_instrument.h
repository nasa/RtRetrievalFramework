#ifndef ILS_INSTRUMENT_H
#define ILS_INSTRUMENT_H
#include "instrument.h"
#include "ils.h"
#include "instrument_correction.h"
#include <vector>

namespace FullPhysics {

/****************************************************************//**
  This is a instrument that uses a Ils object for each spectrometer
  to model the instrument.
*******************************************************************/

class IlsInstrument: public Instrument, public Observer<Ils>, 
		     public Observer<InstrumentCorrection> {
public:
  IlsInstrument(const std::vector<boost::shared_ptr<Ils> >& Ils_list,
		const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >&
		Instrument_correction = 
		std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >());
  virtual ~IlsInstrument() {}
  virtual int number_spectrometer() const {return (int) ils_.size();}
  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const;
  virtual SpectralDomain pixel_spectral_domain(int Spec_index) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index]->pixel_grid();
  }
  virtual std::string band_name(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index]->hdf_band_name();
  }
  virtual std::string hdf_band_name(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index]->hdf_band_name();
  }
  virtual DoubleWithUnit
  ils_half_width(int Spec_index) const 
  {
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index]->ils_half_width();
  }
  virtual void ils_half_width(int Spec_index, DoubleWithUnit& half_width)
  {
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index]->ils_half_width(half_width);
  }
  virtual void print(std::ostream& Os) const;
  virtual void notify_update(const StateVector& Sv) 
  { /* Nothing to do, we have everything done by attached ils. */ }
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const Ils& D)
  { notify_update_do(*this); }
  virtual void notify_update(const InstrumentCorrection& C)
  { notify_update_do(*this); }

  virtual boost::shared_ptr<Instrument> clone() const;

//-----------------------------------------------------------------------
/// Underlying Ils.
//-----------------------------------------------------------------------
  boost::shared_ptr<Ils> ils(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return ils_[Spec_index];
  }

//-----------------------------------------------------------------------
/// Underlying InstrumentCorrection.
//-----------------------------------------------------------------------
  const std::vector<boost::shared_ptr<InstrumentCorrection> >&
  instrument_correction(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return inst_corr[Spec_index];
  }
private:
  std::vector<boost::shared_ptr<Ils> > ils_;
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >
  inst_corr;
};
}
#endif
