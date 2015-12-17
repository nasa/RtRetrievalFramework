#include "solver_iteration_log.h"

using namespace FullPhysics;
using namespace blitz;

void iter_log_add_as_observer(SolverIterationLog& iter_log, ConnorSolver& solver) {
  solver.add_observer(iter_log);
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SolverIterationLog)
.def(luabind::constructor<const boost::shared_ptr<StateVector>&>())
.def("add_as_observer", &iter_log_add_as_observer)
REGISTER_LUA_END()
#endif

void SolverIterationLog::notify_update(const ConnorSolver& solver)
{ 
  blitz::Array<std::string, 1> sv_names = sv_obj->state_vector_name();
  blitz::Array<double, 1> dx(solver.x_update().copy());
  blitz::Array<double, 1> sv_prev(solver.x_solution().copy());
  sv_prev -= dx; // What it was before the dx update
  
  std::stringstream fit_log;
  fit_log << std::endl 
	  << solver.fit_statistic();
  Logger::info() << fit_log.str();

  std::stringstream sv_log;
  sv_log << std::endl
	 << std::setw(SV_PRINT_WIDTH) 
	 << "State Vector"
	 << std::setw(SV_PRINT_WIDTH) 
	 << "Dx" 
	 << "  Name" << std::endl;
  for(int sv_idx = 0; sv_idx < sv_prev.rows(); sv_idx++)
    sv_log << std::setprecision(SV_PRINT_WIDTH-7) // leave room for -/+, ., exponent and space, etc
	   << std::setw(SV_PRINT_WIDTH) 
	   << sv_prev(sv_idx)
	   << std::setw(SV_PRINT_WIDTH) 
	   << dx(sv_idx) << "  "
	   << sv_names(sv_idx)
	   << std::endl;
  Logger::info() << sv_log.str();
}
