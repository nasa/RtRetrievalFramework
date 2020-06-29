// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "IlsFunction.h"
#include "sub_state_vector_array.h"
%}

%fp_shared_ptr(FullPhysics::IlsFunction)
namespace FullPhysics {
  class IlsFunction;
}

%base_import(state_vector)
%import "spectral_domain.i"
%import "sub_state_vector_array.i"

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::IlsFunction>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::IlsFunction>)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>)
%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>;

namespace FullPhysics {

%template(ObservableIlsFunction) FullPhysics::Observable<IlsFunction>;
%template(ObserverIlsFunction) FullPhysics::Observer<IlsFunction>;

class IlsFunction: virtual public StateVectorObserver,
		   public Observable<IlsFunction> {
public:
  virtual ~IlsFunction();
  virtual void add_observer(Observer<IlsFunction>& Obs);
  virtual void remove_observer(Observer<IlsFunction>& Obs);
  virtual boost::shared_ptr<IlsFunction> clone() const = 0;
  std::string print_to_string() const;
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT) const = 0;
  %python_attribute(band_name, virtual std::string);
  %python_attribute(hdf_band_name, virtual std::string);
};

%template(SubStateVectorArrayIlsFunction) 
FullPhysics::SubStateVectorArray<FullPhysics::IlsFunction>;

}

