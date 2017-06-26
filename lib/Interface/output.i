// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "output.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::Output);
%fp_shared_ptr(FullPhysics::OutputDouble);
%fp_shared_ptr(FullPhysics::OutputBlitz1d);
%fp_shared_ptr(FullPhysics::OutputBlitz2d);

namespace FullPhysics {
%feature("director") OutputDouble;
%feature("director") OutputBlitz1d;
%feature("director") OutputBlitz2d;

class OutputDouble : public GenericObject {
public:
  OutputDouble();
  virtual ~OutputDouble();
  virtual double f() const = 0;
};

class OutputBlitz1d  : public GenericObject {
public:
  OutputBlitz1d();
  virtual ~OutputBlitz1d();
  virtual blitz::Array<double,1> f() const = 0;
};

class OutputBlitz2d  : public GenericObject {
public:
  OutputBlitz2d();
  virtual ~OutputBlitz2d();
  virtual blitz::Array<double,2> f() const = 0;
};

%pythoncode %{
class OutputDoubleWrap(OutputDouble):
    def __init__(self, func):
        OutputDouble.__init__(self)
        self.func = func

    def f(self):
        return self.func()

class OutputBlitz1dWrap(OutputBlitz1d):
    def __init__(self, func):
        OutputBlitz1d.__init__(self)
        self.func = func

    def f(self):
        return self.func()

class OutputBlitz2dWrap(OutputBlitz2d):
    def __init__(self, func):
        OutputBlitz2d.__init__(self)
        self.func = func

    def f(self):
        return self.func()

%}
%nodefaultctor Output;
class Output {
public:
  virtual ~Output();
  std::string print_to_string() const;
  void write();
  void write_best_attempt();
  %extend {
    void _register_data_source(const std::string& Nm, 
			       const boost::shared_ptr<OutputDouble>& F)
    { $self->register_data_source(Nm, &FullPhysics::OutputDouble::f, F); }
    void _register_data_source(const std::string& Nm, 
			       const boost::shared_ptr<OutputBlitz1d>& F)
    { $self->register_data_source(Nm, &FullPhysics::OutputBlitz1d::f, F); }
    void _register_data_source(const std::string& Nm, 
			       const boost::shared_ptr<OutputBlitz2d>& F)
    { $self->register_data_source(Nm, &FullPhysics::OutputBlitz2d::f, F); }
  }
  %pythoncode {
def register_double(self, nm, f):
    self._register_data_source(nm, OutputDoubleWrap(f))

def register_array_1d(self, nm, f):
    self._register_data_source(nm, OutputBlitz1dWrap(f))

def register_array_2d(self, nm, f):
    self._register_data_source(nm, OutputBlitz2dWrap(f))
  }
};


}
