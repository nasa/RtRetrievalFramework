// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "composite_initial_guess.h"
%}
%base_import(generic_object)
%base_import(initial_guess)

%fp_shared_ptr(FullPhysics::CompositeInitialGuess);
%fp_shared_ptr(FullPhysics::InitialGuessBuilder);

namespace FullPhysics {

class InitialGuessBuilder : public GenericObject {
public:
  virtual ~InitialGuessBuilder();
  %python_attribute(number_element, virtual int);
  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) const = 0;
  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const = 0;
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const = 0;
  std::string print_to_string() const;
};

class CompositeInitialGuess : public InitialGuess,
			      public InitialGuessBuilder {
public:
  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const;
  %python_attribute(initial_guess, virtual blitz::Array<double, 1>);
  void add_builder(const boost::shared_ptr<InitialGuessBuilder>& B);
  void remove_builder(const boost::shared_ptr<InitialGuessBuilder>& B);
};
}
