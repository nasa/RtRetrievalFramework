#ifndef ITERATIVE_SOLVER_H
#define ITERATIVE_SOLVER_H
#include <vector>
#include <printable.h>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <cost_func.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all iterative optimizers
///
/// This class is the base class for iterative optimizers.  No
/// derivatives of any order appear in the class implementation
/// because not all optimization problem solvers require 
/// derivatives.
///
/// An iterative solver, during the process of minimizing a 
/// cost function, moves from one point to the next in the 
/// parameter space.  After computing a step, some optimization 
/// algorithms use some criteria to decide whether or not to
/// take the step and to move to the next point.
/// In the context of the class hierarchy rooted at this
/// class, a step that is taken is called an accepted step,
/// and a step that is not taken is called a rejected step.
/// All the classes of the hierarchy that implement the solve()
/// method must correctly record all the points and their associated
/// cost function values after taking accepted steps only.  The 
/// initial starting point and its cost value must also be recorded.
///
/// This class is not associated with any problem (from the 
/// problem class hierarchy).  Problems may appear in different
/// formats (even if equivalent): residual/Jacobian, 
/// cost-function/gradient, cost-function only ...
//-----------------------------------------------------------------------

class IterativeSolver : 
    public Printable<IterativeSolver> {

public:


//-----------------------------------------------------------------------
/// \brief Enum type for the status of the iterative solver
//-----------------------------------------------------------------------

  enum status_t {
    SUCCESS,  ///< solve method called and a solution found
    CONTINUE, ///< solve method called but did not converge to a solution
    ERROR,    ///< solve method called but an error was encountered
    UNTRIED   ///< solve method not called yet
  };


//-----------------------------------------------------------------------
/// \brief Constructor
/// 
/// \param[in] max_cost_function_calls 
///            max number of times to evaluate cost function
///
/// \param[in] dx_tol_abs
///            for convergence check, a tolerance for absolute step size,
///            and may be used in different ways by different derived
///            classes
///
/// \param[in] dx_tol_rel
///            for convergence check, a tolerance for relative step size,
///            and may be used in different ways by different derived
///            classes
///
/// \param[in] vrbs
///            if set to true, then a leaf class of this class hierarchy
///            should display some diagnostic information during each
///            iteration
//-----------------------------------------------------------------------

  IterativeSolver(int max_cost_function_calls, 
                  double dx_tol_abs, double dx_tol_rel, bool vrbs)
    : max_cost_f_calls(max_cost_function_calls),
      dX_tol_abs(dx_tol_abs), dX_tol_rel(dx_tol_rel), 
      stat(UNTRIED), verbose(vrbs)
  {}


  virtual ~IterativeSolver() {}


//-----------------------------------------------------------------------
/// \brief Returns the number of the accepted steps.
///
/// \return Number of the accepted steps
//-----------------------------------------------------------------------

  virtual int num_accepted_steps() const
  { return Accepted_points.size()-1; }


//-----------------------------------------------------------------------
/// \brief Returns a vector (std) of accepted points.
///
/// This method returns a std vector of accepted points in 
/// the parameter space.  The initial starting point is always
/// an accepted point.  Then the second accepted point is the
/// point obtained after taking the first accepted step from 
/// the initial (first) point.  The third accepted point is the
/// point obtained after taking the second accepted step from 
/// the second accepted point and so on.
///
/// In other words, if the initial point and all the accepted
/// points after taking the accepted steps are recorded correctly,
/// then 
///   - accepted_points()[0] is the initial starting point,
///   - accepted_points()[1] is the point obtained after taking
///     the first accepted step from the initial point,
///   - accepted_points()[2] is the point obtained after taking
///     the second accepted step from accepted_points()[1] point,
///   - ...
///   - accepted_points()[num_accepted_steps()] is the last 
///     accepted point obtained after the last accepted step.
///
/// Therefore, if the recording of the accepted points is done
/// correctly, and num_accepted_steps() returns 2, then
///   - accepted_points()[1] - accepted_points()[0] is the first
///     accepted step, and
///   - accepted_points()[2] - accepted_points()[1] is the second
///     accepted step
///
/// \return A vector of accepted points
//-----------------------------------------------------------------------

  virtual std::vector< blitz::Array<double, 1> > accepted_points() const
  { return Accepted_points; }


//-----------------------------------------------------------------------
/// \brief Returns a vector (std) of cost
///        function values at accepted points.
///
/// This method returns a std vector of cost function values
/// computed at the accepted points.  In other words, if the 
/// accepted points and the cost function values at these points
/// are recorded correctly, then
///   - cost_at_accepted_points()[0] is the value of the cost 
///     function at accepted_points()[0]
///   - cost_at_accepted_points()[1] is the value of the cost 
///     function at accepted_points()[1]
///   - ...
///   - and finally cost_at_accepted_points()[num_accepted_steps()]
///     is the value of the cost function at 
///     accepted_points()[num_accepted_steps()]
///
/// \return A vector of cost function values at accepted points
//-----------------------------------------------------------------------

  virtual std::vector<double> cost_at_accepted_points() const
  { return Cost_at_accepted_points; }


//-----------------------------------------------------------------------
/// \brief The method that solves the optimization problem.
///
/// The algorithms that solve the optimization problem are 
/// implemented in this method by the leaf classes of the
/// solver class hierarchy.
//-----------------------------------------------------------------------

  virtual void solve() = 0;


//-----------------------------------------------------------------------
/// \brief Returns a value of IterativeSolver::status_t type.
///
/// This method returns the status of the solver.  The status
/// of the solver is initialized to IterativeSolver::UNTRIED,
/// then it must be set to one of the following values by the
/// implemented version of solve() method:
///   - IterativeSolver::SUCCESS
///   - IterativeSolver::CONTINUE
///   - IterativeSolver::ERROR
///
/// Please, read the comments on IterativeSolver::status_t type
/// and its possible values.
///
/// \return Solver status
//-----------------------------------------------------------------------

  virtual status_t status() const
  { return stat; }


//-----------------------------------------------------------------------
/// \brief Returns the string version of the solver status.
///
/// If the method status() returns
///   - IterativeSolver::UNTRIED,
///   - IterativeSolver::SUCCESS,
///   - IterativeSolver::CONTINUE, or
///   - IterativeSolver::ERROR
///
/// then status_str() will return
///   - "UNTRIED",
///   - "SUCCESS",
///   - "CONTINUE", or 
///   - "ERROR"
///
/// respectively.
///
/// \return Solver status in string form
//-----------------------------------------------------------------------

  virtual const char * const status_str() const;


//-----------------------------------------------------------------------
/// \brief Called to record an accepted point
///
/// This method is called to record an accepted point.
/// It is the responsibility of the implementer of the solve()
/// method to record the accepted points.  The accepted points
/// must be recorded in the same order that they are achieved.
///
/// \param[in] point
///            an accepted point in the parameter space
//-----------------------------------------------------------------------

  void record_accepted_point(const blitz::Array<double, 1>& point)
  { Accepted_points.push_back(point); }


//-----------------------------------------------------------------------
/// \brief Called to record the cost function value at an accepted point
///
/// This method is called to record the cost function value at an
/// accepted point.  It is the responsibility of the implementer 
/// of the solve() method to record the cost function values at the 
/// accepted points.  The cost values must be recorded in the same
/// order that they are evaluated.
///
/// \param[in] cost
///            cost funciotn value at an accepted point
///            in the parameter space
//-----------------------------------------------------------------------

  void record_cost_at_accepted_point(double cost)
  { Cost_at_accepted_points.push_back(cost); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "IterativeSolver"; }


protected:

  int max_cost_f_calls;
  double dX_tol_abs, dX_tol_rel;

  status_t stat;
  bool verbose;

private:

  std::vector< blitz::Array<double, 1> > Accepted_points;
  std::vector< double > Cost_at_accepted_points;

};
}
#endif
