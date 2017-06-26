// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "l2_fp_configuration.h"
#include "composite_initial_guess.h"
#include "oco_forward_model.h"
%}
%base_import(generic_object)
%import "logger.i"
%import "forward_model.i"
%import "connor_solver.i"
%import "initial_guess.i"
%import "output.i"
%fp_shared_ptr(FullPhysics::L2FpConfiguration)
namespace FullPhysics {
class L2FpConfiguration : public GenericObject {
public:
  virtual ~L2FpConfiguration();
  %python_attribute(logger, boost::shared_ptr<LogImp>)
  %python_attribute(forward_model, boost::shared_ptr<ForwardModel>)
  %python_attribute(solver, boost::shared_ptr<ConnorSolver>)
  %python_attribute(initial_guess, boost::shared_ptr<InitialGuess>)
  virtual void output(boost::shared_ptr<Output>& OUTPUT,
		      boost::shared_ptr<Output>& OUTPUT) const = 0;
  std::string print_to_string() const;
};
}

