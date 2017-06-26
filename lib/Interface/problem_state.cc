#include <problem_state.h>
#include <fp_exception.h>


using namespace FullPhysics;
using namespace blitz;


bool ProblemState::parameters_different(const Array<double, 1>& x) const
{
  assert_parameter_correct(x);

  bool different = (X.size() <= 0);

  // A test is performed on the parameters. If the absolute
  // relative difference of a single current parameter and
  // its corresponding element from the input parameters vector
  // is larger than abs_rel_diff_tol, then the input parameters
  // are considered to be different (in other words, the new
  // point in the parameter space is considered to be different
  // from the current point in the parameter space.)
  //
  // Most sensitive when less than or equal to zero.
  //
  double abs_rel_diff_tol = 1.0e-12;
  for(int i=0; (i<X.rows()) && !different; i++) {
    if( X(i) == 0.0 )
      different = ( x(i) != 0.0 );
    else
      different = abs((X(i)-x(i))/X(i)) > abs_rel_diff_tol;
  }

  return different;
}

void ProblemState::parameters(const blitz::Array<double, 1>& x)
{
  if( parameters_different(x) ) {
    clear();
    X.reference(x.copy());
  }
}

void ProblemState::assert_parameter_correct(const blitz::Array<double, 1>& x) const
{
  if(x.rows() != expected_parameter_size()) {
    Exception e;
    e << "For this problem state, x does not have the expected size:\n"
      << " Expected size " << expected_parameter_size()
      << " Got size " << x.rows();
    throw e;
  }
}
