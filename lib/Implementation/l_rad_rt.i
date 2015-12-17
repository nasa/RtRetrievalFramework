// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "l_rad_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%base_import(radiative_transfer_single_wn)
%import "spectral_bound.i"
%import "rt_atmosphere.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::LRadRt);

namespace FullPhysics {
class LRadRt : public RadiativeTransferSingleWn {
public:
    LRadRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt,
           const SpectralBound& Spec_bound,
           const blitz::Array<double, 1>& Sza, 
           const blitz::Array<double, 1>& Zen, 
           const blitz::Array<double, 1>& Azm, 
           bool Pure_nadir,
           bool Use_first_order_scatt_calc = true,
           bool Do_second_order = false,
           double Spectrum_spacing = 0.01,
           const PsMode ps_mode = DETECT);
  
    LRadRt(const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
           const boost::shared_ptr<RtAtmosphere>& Atm,
           const SpectralBound& Spec_bound,
           const blitz::Array<double, 1>& Sza, 
           const blitz::Array<double, 1>& Zen, 
           const blitz::Array<double, 1>& Azm, 
           bool Pure_nadir,
           int Number_stokes,
           bool Do_second_order = false,
           int Number_stream = 4,
           double Spectrum_spacing = 0.01,
           const PsMode ps_mode = DETECT);

  %python_attribute(number_stokes, virtual int)
  %python_attribute(number_stream, virtual int)
  %python_attribute(surface_type, virtual int)
  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  %python_attribute(radiative_transfer, boost::shared_ptr<RadiativeTransfer>)
};
}

