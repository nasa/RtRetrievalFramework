#include "connor_solver_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ConnorSolverOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<ConnorSolver>& >())
.def(luabind::constructor<const boost::shared_ptr<ConnorSolver>&, bool>())
REGISTER_LUA_END()
#endif

void ConnorSolverOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalResults/aposteriori_covariance_matrix",
			   &ConnorSolver::aposteriori_covariance,
			   solver);
  out->register_data_source("/RetrievalResults/apriori_covariance_matrix",
			   &ConnorSolver::apriori_covariance,
			   solver);
  out->register_data_source("/RetrievalResults/averaging_kernel_matrix",
  			   &ConnorSolver::averaging_kernel,
  			   solver);
  out->register_data_source("/RetrievalResults/diverging_steps",
  			   &ConnorSolver::number_divergent,
  			   solver);
  out->register_data_source("/RetrievalResults/iterations",
  			   &ConnorSolver::number_iteration,
  			   solver);
  out->register_data_source("/RetrievalResults/outcome_flag",
  			   &ConnorSolver::outcome_flag,
  			   solver);
  out->register_data_source("/RetrievalResults/num_state_vector_elements",
  			   &ConnorSolver::x_solution_size,
  			   solver);


  out->register_data_source("/RetrievedStateVector/state_vector_apriori",
  			   &ConnorSolver::x_apriori,
  			   solver);
  out->register_data_source("/RetrievedStateVector/state_vector_apriori_uncert",
  			   &ConnorSolver::x_apriori_uncertainty,
  			   solver);
  out->register_data_source("/RetrievedStateVector/state_vector_result",
  			   &ConnorSolver::x_solution_zero_unused,
  			   solver);
  out->register_data_source("/RetrievedStateVector/state_vector_aposteriori_uncert",
  			   &ConnorSolver::x_solution_uncertainty,
  			   solver);

  if(write_jacobian)
    out->register_data_source("/RetrievalResults/jacobian",
			     &ConnorSolver::jacobian, solver);
}
