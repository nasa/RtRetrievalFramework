#include "rt_atmosphere.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RtAtmosphere)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for RtAtmosphere.
//-----------------------------------------------------------------------

AccumulatedTimer RtAtmosphere::timer("RtAtmosphere calculations");

//-----------------------------------------------------------------------
/// Return timer information
//-----------------------------------------------------------------------

std::string RtAtmosphere::timer_info() const
{
  std::ostringstream os;
  os << timer;
  return os.str();
}


//-----------------------------------------------------------------------
/// The optical depth for each layer, for the given wave number.
///
/// This variation gives the derivatives with respect to the state vector,
/// this just combines optical_depth with the Jacobian of the intermediate
/// variables given by intermediate_variable.
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \return Optical depth for each layer. This is number_layer() in size
//-----------------------------------------------------------------------
  
ArrayAd<double, 1> 
RtAtmosphere::optical_depth_wrt_state_vector(double wn, int spec_index) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  ArrayAd<double, 1> od(optical_depth_wrt_iv(wn, spec_index));
  Array<double, 3> ivjac(intermediate_variable(wn, spec_index).jacobian());
  ArrayAd<double, 1> res(od.rows(), ivjac.depth());
  res.value() = od.value();
  for(int i = 0; i < res.rows(); ++i)
    res.jacobian()(i, ra) = 
      sum(od.jacobian()(i, ra)(i2) * ivjac(i, ra, ra)(i2,i1),i2);
  return res;
}

//-----------------------------------------------------------------------
/// The single scattering albedo for each layer, for the given wave number.
///
/// This variation gives the derivatives with respect to the state
/// vector, this just combines single_scattering_albedo with the
/// Jacobian of the intermediate variables given by
/// intermediate_variable.
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \return Single scattering albedo for each layer. This is
/// number_layer() in size 
//-----------------------------------------------------------------------
  
ArrayAd<double, 1> 
RtAtmosphere::single_scattering_albedo_wrt_state_vector
(double wn, int spec_index) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  ArrayAd<double, 1> ss(single_scattering_albedo_wrt_iv(wn, spec_index));
  Array<double, 3> ivjac(intermediate_variable(wn, spec_index).jacobian());
  ArrayAd<double, 1> res(ss.rows(), ivjac.depth());
  res.value() = ss.value();
  for(int i = 0; i < res.rows(); ++i)
    res.jacobian()(i, ra) = 
      sum(ss.jacobian()(i, ra)(i2) * ivjac(i, ra, ra)(i2,i1),i2);
  return res;
}

//-----------------------------------------------------------------------
/// The scattering moments for for each layer, for the given wave
/// number.
///
/// The scattering moments use the de Rooij convention for the 6
/// scattering matrix element.
///
/// This variation gives the derivatives with respect to the state
/// vector, this just combines single_scattering_albedo with the
/// Jacobian of the intermediate variables given by
/// intermediate_variable.
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \param nummom Number of moments to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \param numscat Number of scattering matrix elements to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \return Scattering moments for each layer. This is 
///         number_moment + 1 x number_layer() x number scattering
///         matrix elements
//-----------------------------------------------------------------------

ArrayAd<double, 3>
RtAtmosphere::scattering_moment_wrt_state_vector
(double wn, int spec_index, int nummom, int numscat) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4; fifthIndex i5;
  ArrayAd<double, 3> 
    pf(scattering_moment_wrt_iv(wn, spec_index, nummom, numscat));
  Array<double, 3> ivjac(intermediate_variable(wn, spec_index).jacobian());
  ArrayAd<double, 3> res(pf.rows(), pf.cols(), pf.depth(), ivjac.depth());
  res.value() = pf.value();
  res.jacobian() = sum(pf.jacobian()(i1, i2, i3, i5) * ivjac(i2,i5, i4), i5);
  return res;
}

