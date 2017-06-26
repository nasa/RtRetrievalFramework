#include "stokes_coefficient_fraction.h"
#include "ostream_pad.h"

using namespace FullPhysics;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(StokesCoefficientFraction, StokesCoefficient)
.def(luabind::constructor<const blitz::Array<double, 2>&, 
     const blitz::Array<double, 1>&, const blitz::Array<bool, 1>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
StokesCoefficientFraction::StokesCoefficientFraction
(const blitz::Array<double, 2>& Stokes_coeff_parallel,
 const blitz::Array<double, 1>& Coeffs,
 const blitz::Array<bool, 1>& Flag)
: stokes_coeff_parallel(Stokes_coeff_parallel.copy())
{
  if(stokes_coeff_parallel.rows() != Coeffs.rows() ||
     stokes_coeff_parallel.rows() != Flag.rows())
    throw Exception("Stokes_coeff_parallel, Coeff, and Flag all need to have the same number of rows");
  stokes_coeff.resize(stokes_coeff_parallel.shape(), 0);
  stokes_coeff.value() = stokes_coeff_parallel;
  for(int i = 0; i < stokes_coeff.rows(); ++i) {
    stokes_coeff.value()(i, 1) *= (1 - 2 * Coeffs(i));
    stokes_coeff.value()(i, 2) *= (1 - 2 * Coeffs(i));
  }
  init(Coeffs, Flag);
}

void StokesCoefficientFraction::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "StokesCoefficientFraction:\n";
  Os << "  Initial parallel stokes coefficients:\n";
  opad << stokes_coeff_parallel << "\n";
  opad.strict_sync();
  Os << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
  Os << "  Stokes coefficient:\n";
  opad << stokes_coefficient().value();
  opad.strict_sync();
}

boost::shared_ptr<StokesCoefficient> StokesCoefficientFraction::clone() const
{
  return boost::shared_ptr<StokesCoefficient>
    (new StokesCoefficientFraction(stokes_coeff_parallel,
				   coeff.value(),
				   used_flag));
}

void StokesCoefficientFraction::calc_stokes_coeff() const
{
  stokes_coeff.value() = stokes_coeff_parallel;
  stokes_coeff.jacobian() = 0;
  for(int i = 0; i < stokes_coeff.rows(); ++i) {
    stokes_coeff(i, 1) = stokes_coeff(i, 1) * (1 - 2 * coeff(i));
    stokes_coeff(i, 2) = stokes_coeff(i, 2) * (1 - 2 * coeff(i));
  }
}
