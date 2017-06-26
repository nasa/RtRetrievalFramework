// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_solver.h"
#include "ifstream_cs.h"
%}
%base_import(observer)
%import "cost_function.i"
%import "convergence_check.i"
%fp_shared_ptr(FullPhysics::ConnorSolver);
%fp_shared_ptr(FullPhysics::ConnorSolverState);

namespace FullPhysics {
class ConnorSolver;
class ConnorSolverState;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::ConnorSolver>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::ConnorSolver>);

namespace FullPhysics {
%template(ObservableConnor) FullPhysics::Observable<ConnorSolver>;

class ConnorSolver : public Observable<ConnorSolver> {
public:
  ConnorSolver(const boost::shared_ptr<CostFunction>& Cf,
	       const boost::shared_ptr<ConvergenceCheck>& Conv,
	       double Gamma_initial = 0.0,
	       const std::string& Save_test_data = "");
  virtual void add_observer(Observer<ConnorSolver>& Obs);
  virtual void remove_observer(Observer<ConnorSolver>& Obs);
  void save_test_data(const std::string& Fname);
  void test_do_inversion(const std::string& Fname, 
			 blitz::Array<double, 1>& Dx, 
			 blitz::Array<double, 2>& Kt_se_m1_k);
  virtual bool solve(const blitz::Array<double, 1>& Initial_guess, 
		     const blitz::Array<double, 1>& Apriori, 
		     const blitz::Array<double, 2>& Apriori_cov);
  %python_attribute_with_set_ptr(state, ConnorSolverState);
  %python_attribute(aposteriori_covariance_scaled, blitz::Array<double, 2>)
  %python_attribute(aposteriori_covariance, blitz::Array<double, 2>)
  %python_attribute(x_solution_uncertainty, blitz::Array<double, 1>)
  %python_attribute(averaging_kernel, blitz::Array<double, 2>)
  %python_attribute(gamma_last_step, double)
  %python_attribute(number_iteration, int)
  %python_attribute(number_divergent, int)
  %python_attribute(outcome_flag, int)
  %python_attribute(x_apriori, blitz::Array<double, 1>)
  %python_attribute(x_apriori_uncertainty, blitz::Array<double, 1>)
  %python_attribute(x_solution, blitz::Array<double, 1>)
  %python_attribute(x_solution_zero_unused, blitz::Array<double, 1>)
  %python_attribute(apriori_covariance,blitz::Array<double, 2>)
  %python_attribute(jacobian,blitz::Array<double, 2>)
  %python_attribute(fit_statistic, FitStatistic)
  %python_attribute(residual, blitz::Array<double, 1>)
  %python_attribute(residual_covariance_diagonal, blitz::Array<double, 1>)
  %python_attribute(apriori_covariance_inv_norm, blitz::Array<double, 2>)
  %python_attribute(cost_function, boost::shared_ptr<CostFunction>)
  %python_attribute(convergence_check, boost::shared_ptr<ConvergenceCheck>)
  %extend {
    void save_state(const std::string& Fname)
    { std::ofstream os(Fname.c_str());
      $self->to_stream(os);
    }
    void load_state(const std::string& Fname)
    {
      FullPhysics::IfstreamCs in(Fname);
      $self->from_stream(in);
    }
  }
};

class ConnorSolverState: public GenericObject {
public:
  ConnorSolverState(const blitz::Array<double, 1>& X_i, 
		    const blitz::Array<double, 1>& X_a, 
		    const blitz::Array<double, 2>& Apriori_cov_scaled, 
		    const blitz::Array<double, 2>& Sa_m1_scaled,
		    const blitz::Array<double, 1>& Sigma_ap, 
		    double Gamma, double Gamma_last_step, double Gamma_intial,
		    const blitz::Array<double, 1>& Residual, 
		    const blitz::Array<double, 1>& Se, 
		    const blitz::Array<double, 2>& K, 
		    const blitz::Array<double, 2>& Kt_se_m1_k, 
		    const blitz::Array<double, 1>& Dx, 
		    const FitStatistic& Fstat);
  std::string print_to_string() const;
  %python_attribute(x_i, blitz::Array<double, 1>)
  %python_attribute(x_a, blitz::Array<double, 1>)
  %python_attribute(apriori_cov_scaled, blitz::Array<double, 2>)
  %python_attribute(sa_m1_scaled, blitz::Array<double, 2>)
  %python_attribute(sigma_ap, blitz::Array<double, 1>)
  %python_attribute(gamma, double);
  %python_attribute(gamma_last_step, double);
  %python_attribute(gamma_initial, double);
  %python_attribute(residual, blitz::Array<double, 1>)
  %python_attribute(k, blitz::Array<double, 2>)
  %python_attribute(kt_se_m1_k, blitz::Array<double, 2>)
  %python_attribute(se, blitz::Array<double, 1>)
  %python_attribute(dx, blitz::Array<double, 1>)
  %python_attribute(fstat, FitStatistic)
  %pickle_init(1, self.x_i, self.x_a, self.apriori_cov_scaled,
	       self.sa_m1_scaled, self.sigma_ap, self.gamma,
	       self.gamma_last_step, self.gamma_initial,
	       self.residual, self.se, self.k, self.kt_se_m1_k,
	       self.dx, self.fstat)
};
}

// Do this so we can derive from this and have it able to be used by the C++ code
// Defined here since rename does not like being inside of a namespace
%feature("director") FullPhysics::Observer<FullPhysics::ConnorSolver>;
%rename(ObserverConnorSolver) FullPhysics::Observer<FullPhysics::ConnorSolver>;

namespace FullPhysics {
class FullPhysics::Observer<FullPhysics::ConnorSolver> {
public:
  virtual ~Observer();
  virtual void notify_add(ConnorSolver& Obs);
  virtual void notify_remove(ConnorSolver& Obs);
  virtual void notify_update(const ConnorSolver& Obs);
};
}
