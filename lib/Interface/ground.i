// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
#include "sub_state_vector_array.h"
%}

%base_import(observer)
%base_import(state_vector)
%import "sub_state_vector_array.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Ground)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Ground>);
namespace FullPhysics {
  class Ground;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Ground>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Ground>)

namespace FullPhysics {
%template(ObservableGround) FullPhysics::Observable<FullPhysics::Ground>;
%template(ObserverGround) FullPhysics::Observer<FullPhysics::Ground>;

class Ground : public Observable<Ground> {
public:
  virtual ~Ground();
  virtual void add_observer(Observer<Ground>& Obs); 
  virtual void remove_observer(Observer<Ground>& Obs);
  std::string print_to_string() const;
  virtual ArrayAd<double, 1> surface_parameter
    (const double wn, const int spec_index) const = 0;
  virtual boost::shared_ptr<Ground> clone() const = 0;
  virtual void print(std::ostream& Os) const = 0;
};
}

%template(SubStateVectorArrayGround) 
FullPhysics::SubStateVectorArray<FullPhysics::Ground>;
