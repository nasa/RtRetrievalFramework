#include "error_analysis.h"
#include "hdf_file.h"
#include "hdf_sounding_id.h"
#include "absorber_absco.h"
#include <new> // std::bad_alloc
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ErrorAnalysis)
.def(luabind::constructor<const boost::shared_ptr<ConnorSolver>&,
     const boost::shared_ptr<RtAtmosphere>&,
     const boost::shared_ptr<ForwardModel>&>())
.def(luabind::constructor<const boost::shared_ptr<MaxAPosteriori>&,
     const boost::shared_ptr<RtAtmosphere>&,
     const boost::shared_ptr<ForwardModel>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
			     const boost::shared_ptr<AtmosphereOco>& Atm,
			     const boost::shared_ptr<ForwardModel>& Fm)
  : solver(Solver), atm(Atm), fm(Fm)
{
}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
			     const boost::shared_ptr<AtmosphereOco>& Atm,
			     const boost::shared_ptr<ForwardModel>& Fm)
  : max_a_posteriori(Max_a_posteriori), atm(Atm), fm(Fm)
{
}

//-----------------------------------------------------------------------
/// This is variation of ErrorAnalysis that takes a general
/// RtAtmosphere. This will *fail* for anything other than a
/// AtmosphereOco, however it is convenient to have this for use with
/// Lua (which has much more limited knowledge of the class structure).
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
			     const boost::shared_ptr<RtAtmosphere>& Atm,
			     const boost::shared_ptr<ForwardModel>& Fm)
: solver(Solver), atm(boost::dynamic_pointer_cast<AtmosphereOco>(Atm)), 
  fm(Fm)
{
}

//-----------------------------------------------------------------------
/// This is variation of ErrorAnalysis that takes a general
/// RtAtmosphere. This will *fail* for anything other than a
/// AtmosphereOco, however it is convenient to have this for use with
/// Lua (which has much more limited knowledge of the class structure).
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
			     const boost::shared_ptr<RtAtmosphere>& Atm,
			     const boost::shared_ptr<ForwardModel>& Fm)
: max_a_posteriori(Max_a_posteriori), 
  atm(boost::dynamic_pointer_cast<AtmosphereOco>(Atm)), 
  fm(Fm)
{
}

//-----------------------------------------------------------------------
/// In order to match error analysis stuff calculated in the old code,
/// we need to zero out any portion of xco2.gradient() that doesn't 
/// directly depend on the portion of the state that explicitly comes
/// from the CO2 parts (e.g., indirect effect of surface pressure, or
/// other parameters). I don't know if this is actually what we want or
/// not, we'll need to sort this out.
///
/// \todo Sort out if this is actually what we want to be doing.
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::dxco2_dstate() const
{
  FeDisableException disable_fp;
  Array<double, 1> res(xco2().gradient().copy());
  res = where(xco2_state_used(), res, 0.0);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate XCO2 measurement error
///
/// This is equation 3-106 in the ATB
///
/// \todo ATB reference is from draft version, make sure it is still
/// the right reference.
//-----------------------------------------------------------------------

double ErrorAnalysis::xco2_measurement_error() const
{
  FeDisableException disable_fp;
  return sum(dxco2_dstate()(i1) * averaging_kernel()(i1, i3) * 
	     aposteriori_covariance()(i3, i2) * dxco2_dstate()(i2));
}

//-----------------------------------------------------------------------
/// Calculate XCO2 smoothing error
///
/// This is equation 3-107 in the ATB
///
/// \todo ATB reference is from draft version, make sure it is still
/// the right reference.
//-----------------------------------------------------------------------

double ErrorAnalysis::xco2_smoothing_error() const
{
  FeDisableException disable_fp;
  Array<bool, 1> xco2_state_usedv = xco2_state_used();
  // t is a multiplier that turns S into S_a,CO2 as described in the ATB.
  Array<double, 1> t(xco2_state_usedv.shape());
  t = where(xco2_state_usedv, 1, 0);
  Array<double, 2> ak(averaging_kernel());
  Array<double, 2> ak_minus_i(ak.shape());
  ak_minus_i = ak;
  for(int i = 0; i < ak_minus_i.rows(); ++i)
    ak_minus_i(i, i) -= 1;
  Array<double, 2> s(apriori_covariance());
  return sum(dxco2_dstate()(i1) * ak_minus_i(i1, i3) * t(i3) * s(i3, i4) * 
	     t(i4) * ak_minus_i(i2, i4) * dxco2_dstate()(i2));
}

//-----------------------------------------------------------------------
/// Calculate XCO2 interference error
///
/// This is equation 3-108 in the ATB
///
/// \todo ATB reference is from draft version, make sure it is still
/// the right reference.
//-----------------------------------------------------------------------

double ErrorAnalysis::xco2_interference_error() const
{
  FeDisableException disable_fp;
  Array<bool, 1> xco2_state_usedv = xco2_state_used();
  // t is a multiplier that turns S into S_ae as described in the ATB.
  Array<double, 1> t(xco2_state_usedv.shape());
  t = where(xco2_state_usedv, 0, 1);
  Array<double, 2> ak(averaging_kernel());
  Array<double, 2> s(apriori_covariance());
  return sum(dxco2_dstate()(i1) * ak(i1, i3) * t(i3) * s(i3, i4) * t(i4) *
	     ak(i2, i4) * dxco2_dstate()(i2));
  return 0;
}

//-----------------------------------------------------------------------
/// This calculates xco2_gain_vector. It is formally the partial 
/// derivative of retrieved XCO2 with respect to input radiance. It has 
/// the same size & shape as SpectralParameters/measured_radiance.  That 
/// is, one entry per sounding, and each sounding has a vector of 
/// length m, which is the total # of channels across the 3 bands. 
///
/// It can be calculated nearly exactly by using the posterior error 
/// covariance matrix S, the jacobian matrix K, the prior error covariance 
/// matrix Sy, and the pressure weighting function h. It is given by
///
/// (h^t, 0^t) S K^t  Sy^-1
///
/// Where ^-1 denotes inverse, ^t denotes matrix transpose.  The thing on 
/// the front in parentheses is a row vector of length n, where n is the
/// number of state vector elements; it has h for its first 20 elements,
/// and 0 after that.  Is is formally the object that transforms the
/// retrieved state vector into XCO2.  Sy^-1 is calculated internally
/// in the code, but it is also easy to construct. It is a diagonal (m,m)
/// matrix with 1/sigma^2 down the diagonal, where sigma is the 1-sigma
/// prior radiance uncertainty for each channel.  K is the (n,m) matrix
/// of jacobians, and S is the (n,n) posterior error covariance matrix.
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::xco2_gain_vector() const
{
    FeDisableException disable_fp;
    Array<double, 1> result(fm->measured_radiance_all().spectral_range().data().shape()); 
    result = 0;
    
    // Temporary, we need to add handling for a iterative_solver.
    if(!solver) {
      Logger::warning() << "Can not calculate xco2_gain_vector since solver is non ConnorSolver" << '\n';
      return result;
    }

    // Make sure the AbsorberVMR for CO2 is a Absorber Hdf so we can use the
    // pressure weighting function method
    boost::shared_ptr<AbsorberAbsco> absorber = boost::dynamic_pointer_cast<AbsorberAbsco>(atm->absorber_ptr());
    if(absorber != 0) {
        // n = # sv items
        // m = # of radiance points

        Array<bool, 1> xco2_state_usedv = xco2_state_used();
        Array<double, 1> sv_press_weighting(xco2_state_usedv.shape());
        sv_press_weighting = 0;
        // 
        // Copy pressure weighting values where CO2 is located in the statevector
        // Shape: n
        Array<double, 1> press_wgt_grid(absorber->pressure_weighting_function_grid().value());
        int p_wgt_idx = 0;
        for(int sv_idx = 0; sv_idx < xco2_state_usedv.rows() && p_wgt_idx < press_wgt_grid.rows(); sv_idx++) {
            if(xco2_state_usedv(sv_idx)) {
                sv_press_weighting(sv_idx) = press_wgt_grid(p_wgt_idx);
                p_wgt_idx++;
            }
        }

        // Shape: n x n
        Array<double, 2> S(aposteriori_covariance());
   
        // Shape: m
        Array<double, 1> Se(solver->residual_covariance_diagonal());

        // Compute Sy_m1 with 1/sqrt(Se) on the diagonal
        // Shape: m x m
        Array<double, 2> Sy_m1;
	// Don't even try for larger sizes, since we don't want to
	// suck up all the memory of a large system. This threshold is
	// arbitrary, and we can change this in the future if needed.
	if(Se.rows() > 20000) {
	  Logger::error() << "Se is too big for calculating xco2_gain_vector. Se.rows() = " << Se.rows() << '\n';
	  return result;
	}
        try {
            Sy_m1.resize(Se.rows(), Se.rows());
        } catch (std::bad_alloc& ba) {
            // If Se is too big then we won't be able to allocated the memory
            // for this operation and we should fail and return an empty result
            Logger::error() << "Se is too big for calculating xco2_gain_vector. Se.rows() = " << Se.rows() << '\n';
            return result;
        }

        Sy_m1 = 0;
        for(int r_idx = 0; r_idx < Se.rows(); r_idx++) {
            Sy_m1(r_idx, r_idx) = 1 / Se(r_idx);
        }

        // Shape: n x m
        Array<double, 2> Kt(solver->jacobian().transpose(secondDim, firstDim));

        Array<double, 1> res1( sum(sv_press_weighting(i2) * S(i2, i1), i2) );
        Array<double, 1> res2( sum(res1(i2) * Kt(i2, i1), i2) );
        result = sum(res2(i2) * Sy_m1(i2, i1), i2);

    } else {
        Logger::warning() << "Can not calculate xco2_gain_vector since atmosphere->absorber() is not of type AbscoAbsorber" << '\n';
    }

    return result;
}

//-----------------------------------------------------------------------
/// This calculates the matrix "H" as described in the ATB, equation
/// 3-95. 
///
/// \todo ATB reference is from draft version, make sure it is still
/// the right reference.
//-----------------------------------------------------------------------

Array<double, 2> ErrorAnalysis::hmat() const
{
  FeDisableException disable_fp;
  int n = (solver ? solver->x_solution().rows() : 
	   max_a_posteriori->parameters().rows());
  Array<double, 2> res(n, n + 1);
  res = 0;
  res(Range::all(), 0) = dxco2_dstate();
  for(int i = 0; i < n; ++i)
    res(i, i + 1) = 1;
  return res;
}

//-----------------------------------------------------------------------
/// This calculates the matrix "H^T C H" as described in the ATB,
/// after equation 3-96. 
///
/// \todo ATB reference is from draft version, make sure it is still
/// the right reference.
//-----------------------------------------------------------------------

Array<double, 2> ErrorAnalysis::ht_c_h() const
{
  FeDisableException disable_fp;
  Array<double, 2> h(hmat());
  Array<double, 2> res(h.cols(), h.cols());
  res = sum(sum(h(i3, i1) * aposteriori_covariance()(i3, i4) * 
		h(i4, i2), i4), i3);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate xco2_correlation_interf
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::xco2_correlation_interf() const
{
  FeDisableException disable_fp;
  Array<double, 2> ht_c_h_v(ht_c_h());
  Array<double, 1> res(ht_c_h_v.rows() - 1);
  Range r(1, ht_c_h_v.rows() - 1);
  double t = ht_c_h_v(0,0);
  res = where(t * ht_c_h_v(r, r)(i1, i1) > 0, 
	      ht_c_h_v(0, r) / sqrt(t * ht_c_h_v(r, r)(i1, i1)), 0);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate the XCO2 averaging kernel
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::xco2_avg_kernel() const
{
  FeDisableException disable_fp;
  Array<double, 1> full_ak(xco2_avg_kernel_full());
  Array<bool, 1> xco2_state_usedv = xco2_state_used();
  Array<double, 1> res(count(xco2_state_usedv));
  int ind = 0;
  for(int i = 0; i < xco2_state_usedv.rows(); ++i)
    if(xco2_state_usedv(i))
      res(ind++) = full_ak(i);
  return res;
}

//-----------------------------------------------------------------------
/// This the XCO2 averaging kernel for the full state vector.
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::xco2_avg_kernel_full() const
{
  FeDisableException disable_fp;
  Array<double, 2> ak(averaging_kernel());
  Array<double, 1> res(ak.cols());
  res = sum(dxco2_dstate()(i2) * ak(i2, i1), i2);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate the normalized XCO2 averaging kernel
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::xco2_avg_kernel_norm() const
{
  FeDisableException disable_fp;
  Array<double, 2> ak(averaging_kernel());
  Array<double, 1> full(ak.cols());
  full = sum(dxco2_dstate()(i2) * ak(i2, i1), i2);
  full = where(abs(dxco2_dstate()) > 1e-20, full / dxco2_dstate(), 0);
  Array<bool, 1> xco2_state_usedv = xco2_state_used();
  Array<double, 1> res(count(xco2_state_usedv));
  int ind = 0;
  for(int i = 0; i < xco2_state_usedv.rows(); ++i)
    if(xco2_state_usedv(i))
      res(ind++) = full(i);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate the interference smoothing uncertainty
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

Array<double, 1> ErrorAnalysis::interference_smoothing_uncertainty() const
{
  FeDisableException disable_fp;
  Array<double, 1> full_ak(xco2_avg_kernel_full());

  Array<double, 1> res(full_ak.shape());
  Array<double, 2> cov(aposteriori_covariance());
  // Allow cov to be slightly negative, e.g., due to round off error.
  res = (full_ak - dxco2_dstate()) * 
    where(cov(i1, i1) > 0, sqrt(cov(i1, i1)), 0);
  return res;
}

//-----------------------------------------------------------------------
/// Calculate an approximation to the size of the continuum signal
/// where there is no significant atmosphere absorption. We
/// approximate this by finding the 10 highest radiance values and
/// averaging them.
//-----------------------------------------------------------------------

double ErrorAnalysis::signal(int band) const
{
  FeDisableException disable_fp;
  const int nrad = 10;
  SpectralRange rad(fm->measured_radiance(band).spectral_range());
  if(rad.data().rows() ==0)
    return 0;
  Array<double, 1> r(rad.data().copy());
  std::sort(r.data(), r.data() + r.rows()); // Min to max value
  r.reverseSelf(firstDim);	     // Now max to min value
  Range r2(0, std::min(nrad - 1, r.rows() - 1));
  return sum(r(r2) / r2.length());
}

//-----------------------------------------------------------------------
/// Helper class for sort done in noise.
//-----------------------------------------------------------------------
// Don't have Doxygen document this class.
/// @cond
class DataRow {
public:
  DataRow(int R, double V) : row(R), value(V) {}
  int row;
  double value;
};

class DataRowCompare {
public:
  bool operator()(const DataRow& D1, const DataRow& D2) const
  { return D1.value > D2.value; }
};

/// @endcond

//-----------------------------------------------------------------------
/// Calculate an approximation to the size noise of the continuum signal
/// where there is no significant atmosphere absorption. We
/// approximate this by finding the 10 highest radiance values and
/// averaging their noise.
//-----------------------------------------------------------------------

double ErrorAnalysis::noise(int band) const
{
  FeDisableException disable_fp;
  const int nrad = 10;
  SpectralRange rad(fm->measured_radiance(band).spectral_range());
  if(rad.data().rows() < nrad || rad.uncertainty().rows() < nrad)
    return 0;
  std::vector<DataRow> dr;
  dr.reserve(rad.data().rows());
  for(int i = 0; i < rad.data().rows(); ++i)
    dr.push_back(DataRow(i, rad.data()(i)));
  std::sort(dr.begin(), dr.end(), DataRowCompare()); // Max to min value
  double sum_noise = 0;
  for(int i = 0; i < nrad; ++i)
    sum_noise += rad.uncertainty()(dr[i].row);

  return sum_noise / nrad;
}

//-----------------------------------------------------------------------
/// XCO2 uncertainty.
//-----------------------------------------------------------------------

double ErrorAnalysis::xco2_uncertainty() const
{
  FeDisableException disable_fp;
  firstIndex i1; secondIndex i2;
  // Special handling if state vector hasn't been set to anything
  // yet. Return a reasonable value.
  Array<double, 2> cov(aposteriori_covariance());
  if(cov.rows() == 0)
    return 0.0;
  return sqrt(sum(dxco2_dstate()(i1) * cov * dxco2_dstate()(i2)));
}

//-----------------------------------------------------------------------
/// Portion of averaging kernel that relates the part of the state vector
/// that is used by the CO2 VMR calculation.
//-----------------------------------------------------------------------

Array<double, 2> ErrorAnalysis::co2_averaging_kernel() const
{
  FeDisableException disable_fp;
  Array<double, 2> ak(averaging_kernel());
  Array<bool, 1> used(xco2_state_used());
  int sz = count(used);
  Array<double, 2> res(sz, sz);
  int ind1 = 0;
  for(int i = 0; i < used.rows(); ++i)
    if(used(i)) {
      int ind2 = 0;
      for(int j = 0; j < used.rows(); ++j)
	if(used(j)) {
	  res(ind1, ind2) = ak(i, j);
	  ++ind2;
	}
      ++ind1;
    }
  return res;
}
