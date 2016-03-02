// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "sub_state_vector_array.h"
#include "aerosol_property.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(observer)
%base_import(state_vector)
%base_import(generic_object)
%import "pressure.i"
%import "sub_state_vector_array.i"
%import "absorber.i"
%import "relative_humidity.i"

%fp_shared_ptr(FullPhysics::AerosolProperty)
namespace FullPhysics {
  class AerosolProperty;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::AerosolProperty>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::AerosolProperty>);

namespace FullPhysics {
%template(ObservableAerosolProperty) FullPhysics::Observable<AerosolProperty>;
%template(ObserverAerosolProperty) FullPhysics::Observer<AerosolProperty>;

class AerosolProperty : virtual public StateVectorObserver, 
		 public Observable<AerosolProperty> {
public:
  virtual ~AerosolProperty();
  virtual void add_observer(Observer<AerosolProperty>& Obs); 
  virtual void remove_observer(Observer<AerosolProperty>& Obs);
  virtual boost::shared_ptr<AerosolProperty> clone() const = 0;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const = 0;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const = 0;
  std::string print_to_string() const;
};
}
%template(vector_aerosol_property) std::vector<boost::shared_ptr<FullPhysics::AerosolProperty> >;
