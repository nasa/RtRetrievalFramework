#include "convergence_check.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ConvergenceCheck)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void FitStatistic::print(std::ostream& Os) const
{
  Os << "FitStatistic\n"
     << "  number_iteration:  " << number_iteration << "\n"
     << "  number_divergent:  " << number_divergent << "\n"
     << "  outcome:           " << (int) outcome << "\n"
     << "  fit_succeeded:     " << fit_succeeded << "\n"
     << "  d_sigma_sq:        " << d_sigma_sq << "\n"
     << "  d_sigma_sq_scaled: " << d_sigma_sq_scaled << "\n"
     << "  chisq_apriori:     " << chisq_apriori << "\n"
     << "  chisq_measured:    " << chisq_measured << "\n"
     << "  chisq_apriori_fc:  " << chisq_apriori_fc << "\n"
     << "  chisq_measured_fc: " << chisq_measured_fc << "\n"
     << "  gamma2:            " << gamma2() << "\n"
     << "  gamma2_fc:         " << gamma2_fc() << "\n";
}

