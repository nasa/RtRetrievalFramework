#ifndef ILS_FUNCTION_H
#define ILS_FUNCTION_H
#include "state_vector.h"
#include "observer.h"
#include "array_ad.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class models an Instrument Line Shape (ILS) function. This
  returns the response around a given wave number, for a given set of
  wavenumbers. This class is use by IlsConvolution.

  It is not guaranteed that the function is normalized, the calling
  class should normalize this if needed. 
*******************************************************************/

class IlsFunction : virtual public StateVectorObserver,
                    public Observable<IlsFunction> {
public:
  virtual ~IlsFunction() {}
  virtual void add_observer(Observer<IlsFunction>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<IlsFunction>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Return response function.
///
/// Note that is function turns out to be a bit of a bottle neck
/// because it is called so many times. Most of the time the results
/// are the same size from one call to the next, so we pass in the
/// results rather than having this be a return value like we normally
/// do. This avoids recreating the array multiple times. We resize the
/// output, so it is fine if it doesn't happen to be the final result
/// size. But much of the time we avoid and extra allocation and
/// destruction.
///
/// An important optimization is done in IlsConvolution, where instead
/// of calculating dres/dstate we create a short gradient
/// [dwn_center, dscale]. IlsConvolution then applies the chain rule
/// to get the final results in dstate. The flag "jac_optimization"
/// controls this.  
///
/// \param wn_center The wave number of the center of the response
///    function
/// \param wn The wavenumbers to return response function for.
/// \param res Return the response function for each of the wn value.
/// \parmm jacobian_optimization If true, then do the optimization
///    described in this function.  
//-----------------------------------------------------------------------

  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& res,
   bool jacobian_optimization = false) const = 0;

//-----------------------------------------------------------------------
/// Descriptive name of the band.
//-----------------------------------------------------------------------

  virtual std::string band_name() const = 0;

//-----------------------------------------------------------------------
/// In general, the name used in HDF files for a particular band is
/// similar but not identical to the more human readable band_name.
/// For example, with GOSAT we use the HDF field name "weak_co2", but
/// the band name is "WC-Band". This gives the HDF name to use.
///
/// The default implementation just returns the same string as the
/// band name.
//-----------------------------------------------------------------------

  virtual std::string hdf_band_name() const { return band_name();}

  virtual ArrayAd<double, 1> coeff_func() const
  { // Default is no coefficients
    ArrayAd<double, 1> res;
    return res;
  }
  virtual blitz::Array<bool, 1> used_flag_func() const
  { // Default is no coefficients
    blitz::Array<bool, 1> res;
    return res;
  }
};
}
#endif
