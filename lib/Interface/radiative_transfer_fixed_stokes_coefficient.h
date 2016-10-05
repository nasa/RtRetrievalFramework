#ifndef RADIATIVE_TRANSFER_FIXED_STOKES_COEFFICIENT_H
#define RADIATIVE_TRANSFER_FIXED_STOKES_COEFFICIENT_H
#include "radiative_transfer.h"
#include "observer.h"
#include "named_spectrum.h"
#include "stokes_coefficient.h"

namespace FullPhysics {
/****************************************************************//**
  For GOSAT and OCO, we have a set of stokes coefficients to go from
  Stokes vector to radiation. This class captures that common behavior.
*******************************************************************/

class RadiativeTransferFixedStokesCoefficient : public RadiativeTransfer, 
  public Observable<std::vector<boost::shared_ptr<NamedSpectrum> > > {
public:
  virtual ~RadiativeTransferFixedStokesCoefficient() {}

//-----------------------------------------------------------------------
/// Stokes coefficients used to go from Stokes vector to scalar
/// reflectance. 
//-----------------------------------------------------------------------

  const boost::shared_ptr<StokesCoefficient>& stokes_coefficient() const 
  {return stokes_coef;}

//-----------------------------------------------------------------------
/// Number of spectrometer we have.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const 
  { return stokes_coef->stokes_coefficient().rows();}
  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const;
  virtual void print(std::ostream& Os, bool Short_form = false) const 
  {
    Os << *stokes_coef;
  }

  /// Required observable functions
  virtual void add_observer(Observer<std::vector<boost::shared_ptr<NamedSpectrum> > > & Obs) 
  { add_observer_do(Obs); }
  virtual void remove_observer(Observer<std::vector<boost::shared_ptr<NamedSpectrum> > >& Obs) 
  { remove_observer_do(Obs); }

protected:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Stokes_coef The stokes coefficients to go from vector stokes
/// parameters to reflectance. 
//-----------------------------------------------------------------------

  RadiativeTransferFixedStokesCoefficient(
       const boost::shared_ptr<StokesCoefficient>& Stokes_coef)
    : stokes_coef(Stokes_coef) {}

//-----------------------------------------------------------------------
/// Default constructor, derived classes should set up stokes_coef. 
//-----------------------------------------------------------------------

  RadiativeTransferFixedStokesCoefficient() {}

  /// Object to go from stokes vector to reflectance. 
  boost::shared_ptr<StokesCoefficient> stokes_coef;
};
}
#endif
