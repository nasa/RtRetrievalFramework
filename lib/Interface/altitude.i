// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"

%{
#include "altitude.h"
#include "state_vector.h"
#include "pressure.h"
#include "temperature.h"
%}

%fp_shared_ptr(FullPhysics::Altitude)
namespace FullPhysics {
  class Altitude;
}

%import "observer.i"
%import "auto_derivative_with_unit.i"
%import "pressure.i"
%import "temperature.i"
%import "observer.i"

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Altitude>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Altitude>)

namespace FullPhysics {
// Allow these classes to be derived from in Python.
%feature("director") Altitude;

%template(ObservableAltitude) FullPhysics::Observable<Altitude>;
%template(ObserverAltitude) FullPhysics::Observer<Altitude>;

// Note, a class that is derived from in python needs to declare every virtual function that
// can be called on it, even if all that happens is the base class
// to a director is called. This is because this class is used to
// create the SwigDirector class, and this class needs each of the member functions to
// direct things properly. It is *not* necessary to add these function to the underlying
// C++, only that you declare them here.

class Altitude : public Observable<Altitude> {
public:
  virtual ~Altitude();
  std::string print_to_string() const;
  virtual void print(std::ostream& Os) const;
  virtual void add_observer(Observer<Altitude>& Obs);
  virtual void remove_observer(Observer<Altitude>& Obs);
  virtual AutoDerivativeWithUnit<double> 
  altitude(const AutoDerivativeWithUnit<double>& P) 
    const = 0;
  virtual AutoDerivativeWithUnit<double> 
  gravity(const AutoDerivativeWithUnit<double>& P) 
    const = 0;
  virtual boost::shared_ptr<Altitude> clone() const = 0;
  virtual boost::shared_ptr<Altitude> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<Temperature>& Temp) const = 0;
};
}
%template(vector_altitude) std::vector<boost::shared_ptr<FullPhysics::Altitude> >;
