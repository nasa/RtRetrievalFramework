// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "initial_guess_value.h"
%}
%base_import(composite_initial_guess)

%fp_shared_ptr(FullPhysics::InitialGuessValue);

namespace FullPhysics {

class InitialGuessValue : public InitialGuessBuilder {
public:
  virtual ~InitialGuessValue();
  %python_attribute(number_element, int)
  virtual void build_initial_value(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori(blitz::Array<double, 1>& v, int index) const;
  virtual void build_apriori_covariance(blitz::Array<double, 2>& m, 
					int index) const;
  %python_attribute_with_set(apriori, blitz::Array<double, 1>)
  %python_attribute_with_set(apriori_covariance, blitz::Array<double, 2>)
  void apriori_subset(const blitz::Array<bool, 1>& Flag, 
		      const blitz::Array<double, 1>& V);
  void apriori_covariance_subset(const blitz::Array<bool, 1>& Flag, 
				 const blitz::Array<double, 2>& V);
  %python_attribute_with_set(initial_guess, blitz::Array<double, 1>)
};
}
