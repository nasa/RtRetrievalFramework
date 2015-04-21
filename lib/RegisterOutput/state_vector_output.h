#ifndef STATE_VECTOR_OUTPUT_H
#define STATE_VECTOR_OUTPUT_H
#include "register_output_base.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the StateVector class that should be
  written as output.

  Note that somewhat surprisingly, we don't use the StateVector object
  to write out the state vector results and covariance. This is
  because although the Solver may have a solution, the actual
  StateVector may or may not be set to this final value. Instead, we
  output these values from the Solver in SolverOutput, and StateVector
  just handles the state vector element names.

  See the discussion in RegisterOutputBase why this isn't just part of
  the StateVector class.
*******************************************************************/
class StateVectorOutput : public RegisterOutputBase {
public:
  StateVectorOutput(const boost::shared_ptr<StateVector>& Sv) : sv(Sv) {}
  virtual ~StateVectorOutput() {}
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<StateVector> sv;
};
}
#endif
