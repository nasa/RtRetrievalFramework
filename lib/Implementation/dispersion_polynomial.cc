#include "dispersion_polynomial.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(DispersionPolynomial, Dispersion)
.def(luabind::constructor<const blitz::Array<double, 1>&, 
			  const blitz::Array<bool, 1>&,
			  const std::string&,
			  const std::string&,
                          int, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// Units passed seperately since this class subclasses from the 
/// SubStateVectorArray and dispersion values are stored in its array.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial
(const blitz::Array<double, 1>& Coeff,
 const blitz::Array<bool, 1>& Used_flag,
 const Unit& Coeff_unit,
 const std::string& Band_name,
 int Number_pixel,
 bool Is_one_based)
: SubStateVectorArray<Dispersion>(Coeff, Used_flag),
  is_one_based(Is_one_based),
  coeff_unit(Coeff_unit),
  band_name_(Band_name),
  index_array(Number_pixel),
  spectral_index(Number_pixel)
{
  initialize();
}

//-----------------------------------------------------------------------
/// Constructor.
/// Pass units by string because Units class does not compile well
/// with Lua
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial
(const blitz::Array<double, 1>& Coeff,
 const blitz::Array<bool, 1>& Used_flag,
 const std::string& Coeff_unit_name,
 const std::string& Band_name,
 int Number_pixel, bool Is_one_based)
: SubStateVectorArray<Dispersion>(Coeff, Used_flag),
  is_one_based(Is_one_based),
  coeff_unit(Coeff_unit_name),
  band_name_(Band_name),
  index_array(Number_pixel),
  spectral_index(Number_pixel)
{
  initialize();
}

// Initialize class internals
void DispersionPolynomial::initialize() {
  for(int i = 0; i < index_array.rows(); ++i) {
    index_array(i) = i + (is_one_based ? 1 : 0);
    spectral_index(i) = i + 1;
  }
}

// See base class for description.
std::string DispersionPolynomial::state_vector_name_i(int i) const
{
  std::string res = "Instrument Dispersion " + band_name_;
  if(i == 0)
    res += " Offset";
  else if(i == 1)
    res += " Scale";
  else
    res += " Parm " + boost::lexical_cast<std::string>(i + 1);
  return res;
}

// See base class for description.
SpectralDomain DispersionPolynomial::pixel_grid() const
{
  firstIndex i1; secondIndex i2;
  ArrayAd<double, 1> res(index_array.rows(), coeff.number_variable());
  res.value() = 0;
  res.jacobian() = 0;
  for(int i = coeff.rows() - 1; i >= 0; --i) {
    res.value() = res.value() * index_array + coeff.value()(i);
    if(coeff.number_variable() > 0) {
      Array<double, 1> coeff_jac_sub(coeff.jacobian()(i, Range::all()));
      res.jacobian() = res.jacobian() * index_array(i1) + coeff_jac_sub(i2);
    }
  }
  return SpectralDomain(res, spectral_index, coeff_unit);
}

boost::shared_ptr<Dispersion> DispersionPolynomial::clone() const
{
  return boost::shared_ptr<Dispersion>
    (new DispersionPolynomial(coeff.value(), used_flag, coeff_unit,
			      band_name_, index_array.rows(),
			      is_one_based));
}

void DispersionPolynomial::print(std::ostream& Os) const 
{
  Os << "DispersionPolynomial for band " << band_name_ << "\n"
     << "  1 based:  " << (is_one_based ? "True" : "False") << "\n"
     << "  Coeff:    (";
  for(int i = 0; i < coeff.rows() - 1; ++i)
    Os << coeff.value()(i) << ", ";
  Os << coeff.value()(coeff.rows() - 1) << ")\n"
     << "  Retrieve: (";
  for(int i = 0; i < used_flag.rows() - 1; ++i)
    Os << (used_flag(i) ? "true" : "false") << ", ";
  Os << (used_flag(used_flag.rows() - 1) ? "true" : "false") << ")";
}
