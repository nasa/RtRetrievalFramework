#include "aerosol_shape_gaussian.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolShapeGaussian, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<bool, 1>&, 
			  const blitz::Array<double, 1>&,
			  const std::string&,
			  const bool&>())
REGISTER_LUA_END()
#endif

const double AerosolShapeGaussian::min_aod = 1e-9;
// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolShapeGaussian::clone
(const boost::shared_ptr<Pressure>& Pres) const
{
  return boost::shared_ptr<AerosolExtinction>
    (new AerosolShapeGaussian(Pres, used_flag, coeff.value(), 
			      aerosol_name(), linear_aod));
}

void AerosolShapeGaussian::calc_aerosol_extinction() const
{
  // Input parameters
  AutoDerivative<double> desired_aod;
  if (linear_aod) {
    desired_aod = coefficient()(0);
    // Don't let aod go lower than a minimum value. Not clear if this
    // is actually what we want to do, see ticket #2252 for this
    // suggestion. We might instead want to use a constrained solver.
    if(desired_aod < min_aod)
      desired_aod = min_aod;
  } else
    desired_aod = exp(coefficient()(0));

  int ngaussians = int((coefficient().rows() - 1) / 2);

  ArrayAd<double, 1> pressure_grid = pressure()->pressure_grid().value;
  AutoDerivative<double> surface_press = pressure()->surface_pressure().value;

  aext.resize(pressure()->number_level(),
              pressure_grid.number_variable());

  for(int g_idx = 0; g_idx < ngaussians; g_idx++) {
    AutoDerivative<double> p0 = coefficient()(g_idx*2+1);
    AutoDerivative<double> sigma = coefficient()(g_idx*2+2);

    for(int lev = 0; lev < aext.rows(); lev++) {
      AutoDerivative<double> p = pressure_grid(lev) / surface_press;
      AutoDerivative<double> g_eval = exp( -1 * (p - p0) * (p - p0) / (2 * sigma * sigma) );

      // Because its not that easy to init ArrayAd from python to 0.0
      if (g_idx == 0)
         aext(lev) = g_eval;
      else
         aext(lev) = aext(lev) + g_eval;
    }
  }

  AutoDerivative<double> scaling_N = desired_aod / total_aod();

  for(int lev = 0; lev < aext.rows(); lev++)
    aext(lev) = aext(lev) * scaling_N;
}

void AerosolShapeGaussian::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "AerosolShapeGaussian: (" 
     << (linear_aod ? "Linear" : "Logarithmic")
     << ")\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
