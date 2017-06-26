// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "composite_perturbation.h"
%}
%base_import(perturbation)
%fp_shared_ptr(FullPhysics::CompositePerturbation);
%fp_shared_ptr(FullPhysics::PerturbationBuilder);

namespace FullPhysics {

class PerturbationBuilder {
public:
  virtual ~PerturbationBuilder();
  %python_attribute_abstract(number_element, int)
  virtual void build_perturbation(blitz::Array<double, 1>& v, int index) const =0;
};

class CompositePerturbation : public Perturbation {
public:
  %python_attribute(perturbation, blitz::Array<double, 1>)
  void add_builder(const boost::shared_ptr<PerturbationBuilder>& B);
  void remove_builder(const boost::shared_ptr<PerturbationBuilder>& B);
};
}
