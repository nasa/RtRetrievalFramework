// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"
%{
#include "absorber_vmr.h"
%}
%base_import(state_vector)
%base_import(observer)
%import "pressure.i"
%import "auto_derivative.i"

%fp_shared_ptr(FullPhysics::AbsorberVmr)
namespace FullPhysics {
  class AbsorberVmr;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::AbsorberVmr>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::AbsorberVmr>);

namespace FullPhysics {
%template(ObservableAbsorberVmr) FullPhysics::Observable<AbsorberVmr>;
%template(ObserverAbsorberVmr) FullPhysics::Observer<AbsorberVmr>;

class AbsorberVmr : virtual public StateVectorObserver,
	            public Observable<AbsorberVmr> {
public:
  virtual ~AbsorberVmr();
  virtual void add_observer(Observer<AbsorberVmr>& Obs);
  virtual void remove_observer(Observer<AbsorberVmr>& Obs);
  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
  %python_attribute(gas_name, virtual std::string);
  virtual AutoDerivative<double> 
  volume_mixing_ratio(const AutoDerivative<double>& P) const = 0;
  virtual ArrayAd<double, 1> vmr_grid(const Pressure& P) const;
  %python_attribute(state_used, virtual blitz::Array<bool, 1>);
  std::string print_to_string() const;
};
}

%template(vector_absorber_vmr) std::vector<boost::shared_ptr<FullPhysics::AbsorberVmr> >;
