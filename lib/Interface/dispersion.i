// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "dispersion.h"
#include "sub_state_vector_array.h"
%}

%fp_shared_ptr(FullPhysics::Dispersion)
namespace FullPhysics {
  class Dispersion;
}

%base_import(state_vector)
%import "spectral_domain.i"
%import "sub_state_vector_array.i"

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Dispersion>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Dispersion>)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Dispersion>)
%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::Dispersion>;

namespace FullPhysics {

%template(ObservableDispersion) FullPhysics::Observable<Dispersion>;
%template(ObserverDispersion) FullPhysics::Observer<Dispersion>;

class Dispersion: virtual public StateVectorObserver,
		  public Observable<Dispersion> {
public:
  virtual ~Dispersion();
  virtual void add_observer(Observer<Dispersion>& Obs);
  virtual void remove_observer(Observer<Dispersion>& Obs);
  virtual boost::shared_ptr<Dispersion> clone() const = 0;
  %python_attribute(pixel_grid, virtual SpectralDomain);
};

%template(SubStateVectorArrayDispersion) 
FullPhysics::SubStateVectorArray<FullPhysics::Dispersion>;

}

