// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ils_function.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%fp_shared_ptr(FullPhysics::IlsFunction)
namespace FullPhysics {
  class IlsFunction;
}

%base_import(state_vector)
%import "sub_state_vector_array.i"
%import "pressure.i"

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::IlsFunction>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::IlsFunction>)
%template(ObservableIlsFunction) FullPhysics::Observable<FullPhysics::IlsFunction>;
%template(ObserverIlsFunction) FullPhysics::Observer<FullPhysics::IlsFunction>;

%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>)
%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>;

namespace FullPhysics {


class IlsFunction: virtual public StateVectorObserver,
		   public Observable<IlsFunction> {
public:
  virtual ~IlsFunction();
  virtual void add_observer(Observer<IlsFunction>& Obs);
  virtual void remove_observer(Observer<IlsFunction>& Obs);
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT,
   bool jac_optimization=false) const = 0;
  %python_attribute(band_name, virtual std::string);
  %python_attribute(hdf_band_name, virtual std::string);
};

%template(SubStateVectorArrayIlsFunction) 
FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>;

}

