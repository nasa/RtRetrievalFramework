#include "pressure_fixed_level.h"
#include "fp_exception.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;
#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(PressureFixedLevel, Pressure)
.def(luabind::constructor<bool, const boost::shared_ptr<PressureLevelInput>&,
			  double>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Size a new level must be in porportion to the layer above it when
/// resizing pressure grid due to moved surface pressure. Helps prevent
/// small layers causing over representitive VMR values in porportion
/// to a layer size.
//-----------------------------------------------------------------------

const double PressureFixedLevel::new_level_fractional_size = .20;

//-----------------------------------------------------------------------
/// Calculate the new pressure grid. Done each time the surface
/// pressure is updated.
//-----------------------------------------------------------------------

void PressureFixedLevel::calc_pressure_grid() const
{
  const AutoDerivative<double>& surf_p = coefficient()(0);
  const blitz::Array<double, 1>& plevel = press_level->pressure_level();

  blitz::Array<double, 1>::const_iterator lb_iter = std::lower_bound(plevel.begin(), plevel.end(), surf_p);
  int nelem = -1;
  if (lb_iter != plevel.end())
    nelem = lb_iter.position()(0) + 1;

  if(nelem < 2 || nelem > plevel.rows())
    throw Exception("Surface pressure is outside the range of the pressure levels");

  // Group surface pressure as bottom of previous layer if new layer
  // with the lower bound of the surface pressure would be too small
  // But only if doing so would not result in 0 levels
  if (nelem > 2) {
    double psurf_diff = surf_p.value() - plevel(nelem-2);
    double level_frac = (plevel(nelem-1) - plevel(nelem-2)) * 
      new_level_fractional_size;
    if (psurf_diff < level_frac) nelem--;
  }

  pgrid.value.resize(nelem, surf_p.number_variable());
  Range r(0, nelem - 2);
  pgrid.value(r) = plevel(r);
  pgrid.value(nelem - 1) = surf_p;
  pgrid.units = units::Pa;
}

//-----------------------------------------------------------------------
/// Clone a PressureFixedLevel object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<PressureFixedLevel>, although you
/// can of course attach them after receiving the cloned object.
//-----------------------------------------------------------------------

boost::shared_ptr<Pressure> PressureFixedLevel::clone() const
{
  boost::shared_ptr<Pressure> res
    (new PressureFixedLevel(used_flag_value()(0), press_level, 
			    coefficient()(0).value()));
  return res;
}

void PressureFixedLevel::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "PressureFixedLevel:\n"
     << "  Surface pressure:    " << surface_pressure().value.value() << "\n"
     << "  Retrieval Flag:      " << (used_flag_value()(0) ? "True\n": 
				      "False\n")
     << "  Number active level: " << number_active_level() << "\n";
  opad << *press_level;
  opad.strict_sync();
}

