// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This is a simple test class used to test the SWIG interface.

%include "common.i"
%{
#include "pressure.h"
%}
%import "pressure.i"

%{
namespace FullPhysics {
class PressureHolder {
public:
 PressureHolder(const boost::shared_ptr<Pressure>& P) : p_(P) {}
 const boost::shared_ptr<Pressure>& p() const { return p_; }
 void p(const boost::shared_ptr<Pressure>& P) { p_ = P; }
private:
 boost::shared_ptr<Pressure> p_;
};
}
%}

namespace FullPhysics {
class PressureHolder {
public:
  PressureHolder(const boost::shared_ptr<Pressure>& P);
  %python_attribute_with_set(p, boost::shared_ptr<Pressure>)
};
}		
