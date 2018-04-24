// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"
#include "sub_state_vector_array.h"
%}
%base_import(radiative_transfer_fixed_stokes_coefficient)
%import "rt_atmosphere.i"
%fp_shared_ptr(FullPhysics::RadiativeTransferSingleWn);
namespace FullPhysics {
class RadiativeTransferSingleWn : 
  public RadiativeTransferFixedStokesCoefficient {
public:
  %python_attribute_abstract(number_stream, int)
  %python_attribute(atmosphere, boost::shared_ptr<RtAtmosphere>)
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
    (const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index, 
   const ArrayAd<double, 2>& Iv) const = 0;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index, 
   const ArrayAd<double, 2>& Iv) const = 0;
protected:
  RadiativeTransferSingleWn
  (const blitz::Array<double, 2>& Stokes_coef,
   const boost::shared_ptr<RtAtmosphere>& Atm);
  boost::shared_ptr<RtAtmosphere> atm;
};
}
