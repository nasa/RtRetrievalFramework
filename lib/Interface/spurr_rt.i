// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "spurr_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(rt_atmosphere)
%base_import(observer)
%base_import(radiative_transfer_single_wn)
%import "spectral_domain.i"
%fp_shared_ptr(FullPhysics::SpurrRt);

namespace FullPhysics {
class SpurrRt : public RadiativeTransferSingleWn,
		public Observer<RtAtmosphere> {
public:
  SpurrRt(const boost::shared_ptr<RtAtmosphere>& Atm,
	  const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
	  const blitz::Array<double, 1>& Sza, 
	  const blitz::Array<double, 1>& Zen, 
	  const blitz::Array<double, 1>& Azm);
  %python_attribute(number_stokes, virtual int)
  %python_attribute(surface_type, virtual int)
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
