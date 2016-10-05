#ifndef RADIANCE_SCALING_FITTED_H
#define RADIANCE_SCALING_FITTED_H
#include "radiance_scaling.h"
#include "double_with_unit.h"
#include "instrument_correction.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {

/****************************************************************//**
  Implements a fitted radiance scaling correction. This correction
  is retrieved by the state vector.
*******************************************************************/

class RadianceScalingSvFit : public RadianceScaling, 
                              public SubStateVectorArray<InstrumentCorrection> {
public:
//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Coeff - Initial value of scale factor
/// \param Used_flag - If true, we update scale factor by values in
///    StateVector. If false, we hold this fixed and just used the
///    initial value.
/// \param Band_ref - The wavenumber/wavelength that the scaling
///    coefficients are relative to.
/// \param Band_name - Name of band
//-----------------------------------------------------------------------
  
  RadianceScalingSvFit(const blitz::Array<double, 1>& Coeff, 
                        const blitz::Array<bool, 1>& Used_flag,
                        const DoubleWithUnit& Band_ref,
                        const std::string& Band_name)
    : RadianceScaling(Coeff, Band_ref, Band_name),
      SubStateVectorArray<InstrumentCorrection>(Coeff, Used_flag)
  { } 

  virtual ~RadianceScalingSvFit() {}

  virtual std::string state_vector_name_i(int i) const;
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;

  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;

  virtual void notify_update(const StateVector& Sv)
  {
    // Call base notify to handle SV stuff
    SubStateVectorObserver::notify_update(Sv);

    // Then set scaling_coeff to reference true coeff array, which is an ArrayAd
    scaling_coeff.reference(coeff);
  };

  //-----------------------------------------------------------------------
  /// Assumed uncertainty of radiance scaling coefficients 
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> radiance_scaling_coeff_uncertainty() const 
  { 
    blitz::Array<double, 1> res(sv_cov_sub.rows());
    if(res.rows() <= 0) {
      res.resize(1);
      res = 0;
    }
    for(int idx = 0; idx < sv_cov_sub.rows(); idx++) {
      double cov_val = sv_cov_sub(idx, idx);
      res(idx) = (cov_val < 0 ? 0 : sqrt(cov_val)); 
    }
    return res;
  }
};
}
#endif
