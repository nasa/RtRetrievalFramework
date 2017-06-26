#include "pressure_sigma.h"
#include "fp_exception.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(PressureSigma, Pressure)
.def(luabind::constructor<const blitz::Array<double, 1>&,
			  const blitz::Array<double, 1>&,
			  double, bool>())
.def(luabind::constructor<const blitz::Array<double, 1>&,
			  double, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

PressureSigma::PressureSigma(const blitz::Array<double, 1>& A,
			     const blitz::Array<double, 1>& B,
			     double Surface_pressure, bool Pressure_flag)
: a_(A.copy()), b_(B.copy())
{
  if(A.rows() != B.rows())
    throw Exception("A and B need to be the same size in PressureSigma constructor");

  blitz::Array<double, 1> val(1);
  blitz::Array<bool, 1> flag(1);
  val(0) = Surface_pressure;
  flag(0) = Pressure_flag;
  init(val, flag);
  cov.resize(1, 1);
  cov(0,0) = 1;
}

//-----------------------------------------------------------------------
/// Constructor.
/// Creates A and B parameters from the pressure grid passed in. 
/// A becomes all 0 of the same size as Pressure_grid
/// B becomes Pressure_grid / Pressure_grid[-1]
//-----------------------------------------------------------------------

PressureSigma::PressureSigma(const blitz::Array<double, 1>& Pressure_grid,
                             double Surface_pressure, bool Pressure_flag)
{
  set_levels_from_grid(Pressure_grid);

  blitz::Array<double, 1> val(1);
  blitz::Array<bool, 1> flag(1);
  val(0) = Surface_pressure;
  flag(0) = Pressure_flag;
  init(val, flag);
  cov.resize(1, 1);
  cov(0,0) = 1;
}

//-----------------------------------------------------------------------
/// Creates A and B parameters from the pressure grid passed in. 
/// A becomes all 0 of the same size as Pressure_grid
/// B becomes Pressure_grid / Pressure_grid[-1]
//-----------------------------------------------------------------------

void PressureSigma::set_levels_from_grid(const blitz::Array<double, 1>& Pressure_grid)
{
  a_.resize(Pressure_grid.rows());
  a_ = 0.0;

  b_.resize(Pressure_grid.rows());
  b_ = Pressure_grid;
  b_ = b_ / Pressure_grid(Pressure_grid.rows()-1); 
  
  cache_stale = true;
  Observable<Pressure>::notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// Calculate the new pressure grid. Done each time the surface
/// pressure is updated.
//-----------------------------------------------------------------------

void PressureSigma::calc_pressure_grid() const
{
  pgrid.units = units::Pa;
  pgrid.value.resize(b_.rows(), sv_full.number_variable());
  for(int i = 0; i < b_.rows(); ++i) {
    pgrid.value(i) = b_(i) * coeff(0) + a_(i);

    // Ensure that the pressure grid calculated in increasing
    // Since so many linear interpolations rely on this, this error
    // message will make more sense than what would be throw otherwise:
    // X needs to be sorted
    if(i > 0 and pgrid.value(i-1).value() > pgrid.value(i).value()) {
      stringstream err_msg;
      err_msg << "At level " << i << " pressure is smaller: " 
	      << "(" 
	      << b_(i-1) << " * " << coeff(0).value() << " + " << a_(i-1)
	      << ") = "
	      << pgrid.value(i).value()
	      << " than the value at the previous level: " 
	      << "(" 
	      << b_(i) << " * " << coeff(0).value() << " + " << a_(i)
	      << ") = "
	      << pgrid.value(i-1).value();
      throw Exception(err_msg.str());
    }
  }
}

//-----------------------------------------------------------------------
/// Clone a PressureFixedLevel object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<PressureFixedLevel>, although you
/// can of course attach them after receiving the cloned object.
//-----------------------------------------------------------------------

boost::shared_ptr<Pressure> PressureSigma::clone() const
{
  boost::shared_ptr<Pressure> res
    (new PressureSigma(a_, b_, coefficient()(0).value(), 
		       used_flag_value()(0)));
  return res;
}

void PressureSigma::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "PressureSigma:\n"
     << "  Surface pressure:    " << surface_pressure().value.value() << "\n"
     << "  Retrieval Flag:      " << (used_flag_value()(0) ? "True\n": 
				      "False\n")
     << "  a: \n";
  opad << a() << "\n";
  opad.strict_sync();
  Os << "  b: \n";
  opad << b() << "\n";
  opad.strict_sync();
}

