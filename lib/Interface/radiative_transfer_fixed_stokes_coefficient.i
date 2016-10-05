// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "radiative_transfer_fixed_stokes_coefficient.h"
%}

%base_import(radiative_transfer)
%base_import(observer)
%base_import(named_spectrum)
%import "spectral_domain.i"
%import "stokes_coefficient.i"

%fp_shared_ptr(FullPhysics::RadiativeTransferFixedStokesCoefficient);

namespace FullPhysics {
class RadiativeTransferFixedStokesCoefficient: public RadiativeTransfer,
  public Observable<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > > {
public:
  virtual ~RadiativeTransferFixedStokesCoefficient();
  %python_attribute(stokes_coefficient, boost::shared_ptr<StokesCoefficient>);
  %python_attribute(number_spectrometer, virtual int)
  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const;
  virtual void add_observer(Observer<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > > & Obs); 
  virtual void remove_observer(Observer<std::vector<boost::shared_ptr<FullPhysics::NamedSpectrum> > >& Obs); 
protected:
  RadiativeTransferFixedStokesCoefficient(
     const blitz::Array<double, 2>& Stokes_coef);
  RadiativeTransferFixedStokesCoefficient();
  blitz::Array<double, 2> stokes_coef;
};
}
