#ifndef INSTRUMENT_DOPPLER_H
#define INSTRUMENT_DOPPLER_H
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "spectrum_effect_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
 Corrects Spectrum spectral domain to account for the doppler shift
 due to the relative velocity between earth and the instrument.
 This value is retrievable due to the way the class is structured.
*******************************************************************/
class InstrumentDoppler : public SpectrumEffectImpBase {
public:
  InstrumentDoppler(const DoubleWithUnit& Relative_velocity, 
                    const bool Used_flag = false);
  InstrumentDoppler(const double Relative_velocity_value,
                    const std::string& Relative_velocity_units,
                    const bool Used_flag = false);

  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual std::string state_vector_name_i(int i) const
  { return "Instrument Doppler Correction " + boost::lexical_cast<std::string>(i + 1); }
  virtual void print(std::ostream& Os) const;

  virtual std::string name() const { return "instrument_doppler"; }

private:
  Unit vel_units;
};
}
#endif
