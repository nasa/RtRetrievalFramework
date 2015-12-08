#ifndef RADIATIVE_TRANSFER_SINGLE_WN_H
#define RADIATIVE_TRANSFER_SINGLE_WN_H
#include "radiative_transfer_fixed_stokes_coefficient.h"
#include "rt_atmosphere.h"

namespace FullPhysics {
/****************************************************************//**
  This is a RadiativeTransfer that supplies an interface that can be
  called for a single wavenumber. 

  This class mixes in some of the common functionality needed by
  LRadRt and LidortDriver. We may want to pull some of this out
  into separate classes. But right now we only have a few examples of
  a RadiativeTransfer and it isn't at all clear what a useful
  hierarchy would be. Rather than try to predict what we might need in
  the future, we really just have two categories: "The LSI" and "LRad
  and LIDORT". This class is really the later, and I guess we could
  have called it StuffThatIsInCommonWithLRadAndLidort. 
*******************************************************************/

class RadiativeTransferSingleWn : 
    public RadiativeTransferFixedStokesCoefficient 
{
public:
  virtual ~RadiativeTransferSingleWn() {}
  boost::shared_ptr<RtAtmosphere> atmosphere_ptr() const {return atm;}

//-----------------------------------------------------------------------
/// Number of streams to use in processing. Note that Lidort 3.0 used
/// a less common "full streams" that was twice the more commonly used
/// "half streams".  This function returns the later. This is the same
/// as what is used in Lidort 3.5 and LRad. The "full streams" used in
/// Lidort 3.0 would be twice this.
//-----------------------------------------------------------------------

  virtual int number_stream() const = 0;
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const;

//-----------------------------------------------------------------------
/// Calculate stokes vector for the given wavenumber. You can
/// optionally supply a set of intermediate atmosphere variables
/// (e.g., taug, taur, taua_i) to use instead of with atmosphere_ptr()
/// to calculate this.
///
/// \param Wn Wavenumber to calculate for. This should be in cm^-1
/// \param Spec_index The Spectral index
/// \param Iv Optional intermediate variables to use, rather than
///   calculating. 
/// \return The set of stokes coefficients. This is
///      number_stokes() in size.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv = ArrayAd<double, 2>()) const = 0;

//-----------------------------------------------------------------------
/// Calculate stokes vector and Jacobian for the given wavenumber. You can
/// optionally supply a set of intermediate atmosphere variables
/// (e.g., taug, taur, taua_i) to use instead of with atmosphere_ptr()
/// to calculate this.
///
/// \param Wn Wavenumber to calculate for. This should be in cm^-1
/// \param Spec_index The Spectral index
/// \param Iv Optional intermediate variables to use, rather than
///   calculating.
/// \return The set of stokes coefficients. This is
///      number_stokes() in size.
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv = ArrayAd<double, 2>()) const = 0;
  virtual void print(std::ostream& Os, bool Short_form = false) const;

  const boost::shared_ptr<RtAtmosphere>& atmosphere() const
  { return atm; }
protected:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Stokes_coef The stokes coefficients to go from vector stokes
/// parameters to reflectance. This should be number_spectrometer() x 4.
/// \param Atm The RtAtmosphere to use.
//-----------------------------------------------------------------------

  RadiativeTransferSingleWn
  (const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
   const boost::shared_ptr<RtAtmosphere>& Atm)
    : RadiativeTransferFixedStokesCoefficient(Stokes_coef),
      atm(Atm) {}

//-----------------------------------------------------------------------
/// Default constructor. Derived classes should fill in atm and stokes_coef
//-----------------------------------------------------------------------
  RadiativeTransferSingleWn() {}

  boost::shared_ptr<RtAtmosphere> atm;
};
}
#endif
