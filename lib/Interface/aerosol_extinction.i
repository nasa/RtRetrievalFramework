// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "aerosol_extinction.h"
%}
%base_import(observer)
%base_import(state_vector)
%import "pressure.i"
%import "auto_derivative.i"

%fp_shared_ptr(FullPhysics::AerosolExtinction)

namespace FullPhysics {
  class AerosolExtinction;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::AerosolExtinction>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::AerosolExtinction>);
namespace FullPhysics {
%template(ObservableAerosolExtinction) FullPhysics::Observable<AerosolExtinction>;
%template(ObserverAerosolExtinction) FullPhysics::Observer<AerosolExtinction>;

class AerosolExtinction : virtual public StateVectorObserver, 
		 public Observable<AerosolExtinction> {
public:
  virtual ~AerosolExtinction();
  virtual void add_observer(Observer<AerosolExtinction>& Obs); 
  virtual void remove_observer(Observer<AerosolExtinction>& Obs);
  virtual boost::shared_ptr<AerosolExtinction> clone() const = 0;
  virtual boost::shared_ptr<AerosolExtinction> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
  virtual AutoDerivative<double> extinction_for_layer(int i) const = 0;
  %python_attribute_abstract(aerosol_name, std::string)
  std::string print_to_string() const;
  virtual void print(std::ostream& Os) const;

};
}

%template(vector_aerosol_extinction) std::vector<boost::shared_ptr<FullPhysics::AerosolExtinction> >;

