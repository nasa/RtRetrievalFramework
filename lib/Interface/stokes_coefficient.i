// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "stokes_coefficient.h"
%}

%base_import(state_vector)
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::StokesCoefficient)
namespace FullPhysics {
  class StokesCoefficient;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::StokesCoefficient>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::StokesCoefficient>)

%template(ObservableStokesCoefficient) FullPhysics::Observable<FullPhysics::StokesCoefficient>;
%template(ObserverStokesCoefficient) FullPhysics::Observer<FullPhysics::StokesCoefficient>;

namespace FullPhysics {
class StokesCoefficient : virtual public StateVectorObserver, 
		 public Observable<StokesCoefficient> {
public:
  virtual ~StokesCoefficient();
  virtual void add_observer(Observer<StokesCoefficient>& Obs);
  virtual void remove_observer(Observer<StokesCoefficient>& Obs);
  %python_attribute_abstract(stokes_coefficient, ArrayAd<double, 2>)
  virtual boost::shared_ptr<StokesCoefficient> clone() const = 0;
  std::string print_to_string() const;
};
}

