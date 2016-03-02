#ifndef AEROSOL_PROPERTY_H
#define AEROSOL_PROPERTY_H
#include "state_vector.h"
#include "pressure.h"
#include "relative_humidity.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This gives the Aerosol properties for an Aerosol.

  Our current AerosolProperty - AerosolPropertyHdf - doesn't make any
  use of our StateVector, we don't have aerosol properties in it.
  But we may want to have this in the future, so we've made this class 
  a StateVectorObserver.

  Other objects may depend on the AerosolProperty, and should be updated
  when the AerosolProperty is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the AerosolProperty is updated.

  When implementing a new class, you almost always will want to derive
  from AerosolPropertyImpBase rather than from this class. See that
  class for a description.
*******************************************************************/
class AerosolProperty : virtual public StateVectorObserver,
			public Observable<AerosolProperty> {
public:
  virtual ~AerosolProperty() {}
  virtual void add_observer(Observer<AerosolProperty>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<AerosolProperty>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Clone a AerosolProperty object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<AerosolProperty>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" AerosolProperty object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AerosolProperty> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a pressure to use. The intent is that
/// the pressure has been cloned from the original pressure (although
/// this class has no way to verify this). This allows sets of objects
/// to be cloned using a common Pressure clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const = 0;

//-----------------------------------------------------------------------
/// Return extinction coefficient for the given wave number, for each
/// layer. 
/// \param wn - Wavenumber
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) 
    const = 0;

//-----------------------------------------------------------------------
/// Return scattering coefficient for the given wave number for each layer.
/// \param wn - Wavenumber
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn)
    const = 0;

//-----------------------------------------------------------------------
/// Return phase function moments for the given wave number for each layer.
///
/// Note that we use the de Rooij convention for the scattering matrix
/// moments. 
///
/// \param wn Wavenumber
/// \param nmom Optional number of moments to return. Default is all
///     moments.
/// \param nscatt Optional number of scattering elements to
///     return. Default is all of them.
/// \return Phase function moment. This is nmom + 1 x nlayer x number
/// scattering elements. 
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const = 0;
};
}
#endif
