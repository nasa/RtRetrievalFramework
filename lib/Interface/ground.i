// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

namespace FullPhysics {
  class Ground;
}

%base_import(observer)
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Ground)
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Ground>)

namespace FullPhysics {
%template(ObservableGround) FullPhysics::Observable<Ground>;

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
