#ifndef FLUOR_EFFECT_IMP_BASE_H
#define FLUOR_EFFECT_IMP_BASE_H
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "spectrum_effect_imp_base.h"
#include "atmosphere_oco.h"
#include "stokes_coefficient.h"

namespace FullPhysics {
/****************************************************************//**
 Implements adding the effect of fluorescence to A-Band spectrum
 by using a retrievable across the band parametrization of the
 effect.
*******************************************************************/
class FluorescenceEffect : public SpectrumEffectImpBase {
public:
  FluorescenceEffect(const blitz::Array<double, 1>& Coeff,
                     const blitz::Array<bool, 1>& Used_flag,
                     const boost::shared_ptr<RtAtmosphere>& Atm,
		     const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                     const DoubleWithUnit& Lza, 
                     const int Spec_index,
                     const DoubleWithUnit& Reference,
                     const Unit& Retrieval_unit);

  virtual void apply_effect(Spectrum& Spec,
		    const ForwardModelSpectralGrid& Forward_model_grid) const;
  virtual ArrayAd<double, 1> contribution() const { return f_contrib_ad; }

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual std::string state_vector_name_i(int i) const
  { return "Fluorescence Surface Coefficient " + boost::lexical_cast<std::string>(i + 1); }
  virtual void print(std::ostream& Os) const;

  virtual std::string name() const { return "fluorescence_effect"; }

  //-----------------------------------------------------------------------
  /// Fluorescence value at reference point 
  //-----------------------------------------------------------------------

  double fluorescence_at_reference() const { return coeff.value()(0); }

  //-----------------------------------------------------------------------
  /// Assumed uncertainty of fluorescence at reference point 
  //-----------------------------------------------------------------------

  double fluorescence_at_reference_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 1)
      return 0;
    double t = sv_cov_sub(0,0);
    return (t < 0 ? 0 : sqrt(t)); 
  }

  //-----------------------------------------------------------------------
  /// Fluorescence slope across band 
  //-----------------------------------------------------------------------

  double fluorescence_slope() const { return coeff.value()(1); }

  //-----------------------------------------------------------------------
  /// Assumed uncertainty of fluorescence slope 
  //-----------------------------------------------------------------------

  double fluorescence_slope_uncertainty() const 
  { 
    if(sv_cov_sub.rows() < 2)
      return 0;
    double t = sv_cov_sub(1,1);
    return (t < 0 ? 0 : sqrt(t)); 
  }

private:
  DoubleWithUnit lza;
  DoubleWithUnit reference;
  Unit retrieval_unit;
  int spec_index;

  // Fluorescence contribution
  mutable ArrayAd<double, 1> f_contrib_ad;
  
  // We need to use AtmosphereOco specific methods
  boost::shared_ptr<AtmosphereOco> atm_oco;

  // The fluorescence contribution needs to be
  // scaled by the stokes value
  boost::shared_ptr<StokesCoefficient> stokes_coef;
};
}
#endif
