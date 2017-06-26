#include "lsi_rt.h"
#include "ifstream_cs.h"
#include "bin_map.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<RadiativeTransfer>
lsi_rt_create(const boost::shared_ptr<RadiativeTransfer>& Rt_low,
	      const boost::shared_ptr<RadiativeTransfer>& Rt_high,
	      const HdfFile& Config_file,
	      const std::string& Lsi_group)
{
  boost::shared_ptr<RadiativeTransferSingleWn> rt_lows =
    boost::dynamic_pointer_cast<RadiativeTransferSingleWn>(Rt_low);
  boost::shared_ptr<RadiativeTransferSingleWn> rt_highs =
    boost::dynamic_pointer_cast<RadiativeTransferSingleWn>(Rt_high);
  return boost::shared_ptr<RadiativeTransfer>
  (new LsiRt(rt_lows, rt_highs, Config_file, Lsi_group));
}

REGISTER_LUA_DERIVED_CLASS(LsiRt, RadiativeTransfer)
.scope
[
 luabind::def("create", &lsi_rt_create)
]
REGISTER_LUA_END()
#endif

// Threshold for switching from log linear to linear interpolation.
const double log_threshold = 100;

//-----------------------------------------------------------------------
/// Create a object that uses the low stream RT + LSI corrections
/// based on the low and high stream RT.
///
/// Unless you are doing something unusual, you'll want the low and
/// high stream RadiativeTransfer classes to use the same stokes
/// coefficients, StateVector and Atmosphere. However, this class
/// doesn't require that these be the same - if they are different it
/// assumes you know what you are doing. In any case, this class will
/// use the same stokes coefficients, StateVector and Atmosphere as
/// the High_stream_rt.
///
/// This version reads the old ASCII format. This format has comments
/// starting with "#" and going to the end of the line. The data is 
/// <number boundaries band 0> <boundary 1> <boundary 2> ...
/// <number boundaries band 1> ...
///
/// \param Low_stream_rt Low stream RadiativeTransfer object.
/// \param High_stream_rt High stream RadiativeTransfer object.
/// \param Lsi_fname Name of ASCII file containing LSI parameters
/// \param Water_vapor_fraction_threshold Threshold used to determine
///    if a particular spectral point is a water vapor line. We 
///    calculate wv_od / total_od and compare to this threshold.
//-----------------------------------------------------------------------

LsiRt::LsiRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Low_stream_rt,
	     const boost::shared_ptr<RadiativeTransferSingleWn>& High_stream_rt,
	     const std::string& Lsi_fname,
	     double Water_vapor_fraction_threshold)
: RadiativeTransferFixedStokesCoefficient(High_stream_rt->stokes_coefficient()),
  low_stream_rt(Low_stream_rt), high_stream_rt(High_stream_rt),
  wv_threshold(Water_vapor_fraction_threshold)
{
  if(low_stream_rt->number_stokes() != high_stream_rt->number_stokes())
    throw Exception("Low stream and high stream need to have the same number of stokes parameters");

  IfstreamCs in(Lsi_fname);
  if(!in.good())
    throw Exception("Trouble reading LSI file " + Lsi_fname);
  in.exceptions(std::ifstream::badbit);
  for(int i = 0; i < number_spectrometer(); ++i) {
    std::string main_gas_in;
    int nbound;
    in >> main_gas_in >> nbound;
    main_gas.push_back(main_gas_in);
    std::vector<double> t(nbound);
    for(int j = 0; j < nbound; ++j)
      in >> t[j];
    optical_depth_boundary.push_back(t);
  }
}

//-----------------------------------------------------------------------
/// Create a object that uses the low stream RT + LSI corrections
/// based on the low and high stream RT.
///
/// Unless you are doing something unusual, you'll want the low and
/// high stream RadiativeTransfer classes to use the same stokes
/// coefficients, StateVector and Atmosphere. However, this class
/// doesn't require that these be the same - if they are different it
/// assumes you know what you are doing. In any case, this class will
/// use the same stokes coefficients, StateVector and Atmosphere as
/// the High_stream_rt.
///
/// This version reads the HDF format. This takes a HdfFile and
/// optional group. It reads the fields in that group names 
/// "optical_depth_boundary_1", "optical_depth_boundary_2", etc.
///
/// \param Low_stream_rt Low stream RadiativeTransfer object.
/// \param High_stream_rt High stream RadiativeTransfer object.
/// \param Config_file HDF file that contains configuration
///    information
/// \param Lsi_group The group that contains the LSI data to read.
/// \param Water_vapor_fraction_threshold Threshold used to determine
///    if a particular spectral point is a water vapor line. We 
///    calculate wv_od / total_od and compare to this threshold.
//-----------------------------------------------------------------------

LsiRt::LsiRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Low_stream_rt,
	     const boost::shared_ptr<RadiativeTransferSingleWn>& High_stream_rt,
	     const HdfFile& Config_file,
	     const std::string& Lsi_group,
	     double Water_vapor_fraction_threshold)
: RadiativeTransferFixedStokesCoefficient(High_stream_rt->stokes_coefficient()),
  low_stream_rt(Low_stream_rt), high_stream_rt(High_stream_rt),
  wv_threshold(Water_vapor_fraction_threshold)
{
  if(low_stream_rt->number_stokes() != high_stream_rt->number_stokes())
    throw Exception("Low stream and high stream need to have the same number of stokes parameters");

  Array<std::string, 1> mg = Config_file.read_field<std::string, 1>(Lsi_group + "/main_gas");
  for(int i = 0; i < number_spectrometer(); ++i) {
    Array<double, 1> d = Config_file.read_field<double, 1>
      (Lsi_group + "/optical_depth_boundary_" + 
       boost::lexical_cast<std::string>(i + 1));
    optical_depth_boundary.push_back(std::vector<double>(d.begin(), d.end()));
    main_gas.push_back(mg(i));
  }
}

void LsiRt::print(std::ostream& Os, bool Short_form) const 
{
  OstreamPad opad(Os, "    ");
  Os << "LsiRt\n";
  Os << "  ";
  RadiativeTransferFixedStokesCoefficient::print(opad, Short_form);
  opad.strict_sync();
  Os << "\n  Water vapor fraction threshold: " << wv_threshold << "\n";
  for(int i = 0; i < (int) optical_depth_boundary.size(); ++i) {
    Os << "  Main gas[" << i << "]: " << main_gas[i] << "\n";
    Os << "  Optical depth boundary[" << i << "]:\n";
    Os << "   (";
    for(int j = 0; j < (int) optical_depth_boundary[i].size(); ++j) {
      Os << optical_depth_boundary[i][j];
      if(j != (int) optical_depth_boundary[i].size() - 1)
	Os << ", ";
      else
	Os << ")\n";
    }
  }
  Os << "\n  High stream:\n";
  high_stream_rt->print(opad, Short_form);
  opad.strict_sync();
  Os << "  Low stream:\n";
  low_stream_rt->print(opad, true);
  opad.strict_sync();
}

// See base class for description
blitz::Array<double, 2> LsiRt::stokes(const SpectralDomain& Spec_domain,
				      int Spec_index) const

{
  FunctionTimer ft(timer.function_timer(true));
  Logger::info() << "RT for band " << Spec_index + 1 << "\n";
  calc_correction(Spec_domain, Spec_index, false, false);
  for(int i = 0; i < Spec_domain.data().rows(); ++i) {
    Array<AutoDerivative<double>, 1> err_est;
    if(gas_opd.value()(i) > log_threshold)
      err_est.reference(linear_interp[wv_index(i)](gas_opd(i)));
    else
      err_est.reference(log_linear_interp[wv_index(i)](gas_opd(i)));
    // Intensity is scaled correction, others are additive.
    stokes_only(i,0) = stokes_only(i,0) / (1 + err_est(0).value());
    for(int j = 1; j < stokes_only.cols(); ++j)
      stokes_only(i,j) = stokes_only(i,j) - err_est(j).value();
  }
  Logger::info() << low_stream_rt->atmosphere_ptr()->timer_info();
  return stokes_only;
}

// See base class for description
ArrayAd<double, 2> LsiRt::stokes_and_jacobian
(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));
  Logger::info() << "RT + Jac for band " << Spec_index + 1 << "\n";
  calc_correction(Spec_domain, Spec_index, true, false);
  for(int i = 0; i < Spec_domain.data().rows(); ++i) {
    Array<AutoDerivative<double>, 1> err_est;
    if(gas_opd.value()(i) > log_threshold)
      err_est.reference(linear_interp[wv_index(i)](gas_opd(i)));
    else
      err_est.reference(log_linear_interp[wv_index(i)](gas_opd(i)));
    // Intensity is scaled correction, others are additive.
    stokes_and_jac(i,0) = stokes_and_jac(i,0) / (1 + err_est(0));
    for(int j = 1; j < stokes_and_jac.cols(); ++j)
      stokes_and_jac(i,j) = stokes_and_jac(i,j) - err_est(j);
  }
  Logger::info() << low_stream_rt->atmosphere_ptr()->timer_info();
  return stokes_and_jac;
}

//-----------------------------------------------------------------------
/// Normally we calculate both the low streams stokes parameters, the
/// LSI correction, and we apply them. However for testing it can be
/// useful to calculate just the LSI correction. This is much faster
/// to calculate, and allows us to test the LSI without doing the full
/// RT calculation.
//-----------------------------------------------------------------------

ArrayAd<double, 2> LsiRt::correction_only
(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));
  Logger::info() << "LSI corrrection only for band " << Spec_index + 1 << "\n";
  calc_correction(Spec_domain, Spec_index, false, true);
  ArrayAd<double, 2> res;
  bool first = true;
  for(int i = 0; i < Spec_domain.data().rows(); ++i) {
    ArrayAd<double, 1> err_est;
    if(gas_opd.value()(i) > log_threshold)
      err_est.reference
	(ArrayAd<double, 1>(linear_interp[wv_index(i)](gas_opd(i))));
    else
      err_est.reference
	(ArrayAd<double, 1>(log_linear_interp[wv_index(i)](gas_opd(i))));
    if(first) {
      res.resize(Spec_domain.data().rows(), err_est.rows(), err_est.number_variable());
      first = false;
    }
    res(i, Range::all()) = err_est;
  }
  Logger::info() << low_stream_rt->atmosphere_ptr()->timer_info();
  return res;
}

//-----------------------------------------------------------------------
/// Calculate the LSI correction. This fills in the variables gas_opd,
/// wv_index, linear_interp and log_linear_interp.
///
/// We normally also calculate either stokes_only or stokes_and_jac
/// which is the low streams results at every Wn. We do this just
/// because we are already calculating all the atmosphere properties
/// and if we intermingle the low streams code we can avoid a second
/// generation of all the properties. However for testing purpose it
/// can be useful to skip this step and just calculate the corrections
/// since this is much quicker.
///
/// \param Spec_domain The list of wavenumbers/wavelengths to process.
/// \param Spec_index The spectral index we are processing.
/// \param Calc_jacobian If true we calculate we fill in
///    stokes_and_jac, otherwise we do stokes_only.
/// \param Skip_stokes_calc If true we don't fill in either
///    stokes_only or stokes_and_jac. We just calculate the LSI
///    corrections.  This is useful for testing
//-----------------------------------------------------------------------

void LsiRt::calc_correction(const SpectralDomain& Spec_domain,
			    int Spec_index, bool Calc_jacobian, 
			    bool Skip_stokes_calc) const
{
  Range ra(Range::all());
  Array<double, 1> wn(Spec_domain.wavenumber());
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);
  const RtAtmosphere& atm = *low_stream_rt->atmosphere_ptr();
  const std::vector<double>& odb = optical_depth_boundary[Spec_index];

//-----------------------------------------------------------------------
// Go through all the wave numbers and collect the atmosphere
// properties, bin them by the optical_depth_boundary bins, and
// average them. 
//
// We collect information for "water vapor wavenumbers" and
// for regular wavenumbers separately. By convention, there are two of
// each kind of bin, and index "0" is used for water and index "1"
// used for normal wave numbers.
//
// There are a number of variables created here:
//
// atm_sum - 
//   Sum of all the intermediate variables (e.g., taug, taur,
//   taua_i). We use this to create the average value used in the
//   correction calculation. This is indexed by wv_index and then
//   optical depth range.
// log_opd_sum -
//   Sum of the log(gas_opd). This is indexed by wv_index and then
//   optical depth range.
// cnt -
//   Number of entries added to atm_sum and log_opd_sum. May be 0 if
//   we don't have any data for a particular bin.
// wv_index -
//   0 if this is a water vapor wave number, 1 otherwise. This is the
//   same size as wn.
// gas_opd
//   The optical depth for the wavenumber. This is the sum of the
//   water vapor gas column optical depth and the main gas column
//   optical depth.
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Initialize each of the variables mentioned above. We start
// everything out at 0.
//
// Note that "BinMap" is a class that handles binning. If "b" is a bin
// map, then b[x] return the bin that x falls into.
//-----------------------------------------------------------------------

  // Get first intermediate variable so we can figure out the size of
  // the intermediate variables (e.g., number_layer x number
  // variables). The initialize this to all 0s.

  ArrayAd<double, 2> 
    init_atm_sum_value(atm.intermediate_variable(wn(0), Spec_index).copy());
  init_atm_sum_value = 0;
  // Need at least 1 or else FM only mode doesn't work
  int numvar = std::max(init_atm_sum_value.number_variable(), 1);
  std::vector<BinMap<int> > cnt; 
  std::vector<BinMap<AutoDerivative<double> > > log_opd_sum; 
  std::vector<BinMap<ArrayAd<double, 2> > > atm_sum;
  for(int i = 0; i < 2; ++i) {
    cnt.push_back(BinMap<int>(odb.begin(), odb.end(), 0));
    log_opd_sum.push_back
      (BinMap<AutoDerivative<double> >(odb.begin(), odb.end(), 
				       AutoDerivative<double>(0.0)));
    atm_sum.push_back(BinMap<ArrayAd<double, 2> > (odb.begin(), odb.end(), 
      boost::bind(&ArrayAd<double, 2>::copy, &init_atm_sum_value)));
  }

  gas_opd.resize(wn.rows(), numvar);
  wv_index.resize(wn.rows());
  if(!Skip_stokes_calc) {
    if(Calc_jacobian)
      stokes_and_jac.resize(wn.rows(), number_stokes(), numvar);
    else
      stokes_only.resize(wn.rows(), number_stokes());
  }
  
//-----------------------------------------------------------------------
// Now, go through and collect the sums of all the things we are
// averaging. 
//-----------------------------------------------------------------------

  for(int i = 0; i < wn.rows(); ++i) {
    // We want to move this to metadata
    AutoDerivative<double> maingas_opd = 
      atm.column_optical_depth(wn(i), Spec_index, main_gas[Spec_index]);
    AutoDerivative<double> wv_opd = 
      atm.column_optical_depth(wn(i), Spec_index, "H2O");
    gas_opd(i) = maingas_opd + wv_opd;
    // Algorithm can't actually handle gas_opd == 0. We get this
    // sometimes if we go past the end of the ABSCO tables, so just
    // pick a small value if we are zero.
    if(gas_opd(i).value() < 1e-10)
      gas_opd(i) = 1e-10;
    wv_index(i) = (wv_opd / gas_opd(i) >= wv_threshold ? 0 : 1);
    ArrayAd<double, 2> atmv(atm.intermediate_variable(wn(i), Spec_index));
    ++cnt[wv_index(i)][gas_opd.value()(i)];
    log_opd_sum[wv_index(i)][gas_opd.value()(i)] += log(gas_opd(i));
    atm_sum[wv_index(i)][gas_opd.value()(i)].value() += atmv.value();
    atm_sum[wv_index(i)][gas_opd.value()(i)].jacobian() += atmv.jacobian();

    // Since we are already calculating the atmosphere parameters, go
    // ahead and collect low streams data if requested.
    if(!Skip_stokes_calc) {
      if(Calc_jacobian)
	stokes_and_jac(i, ra) = low_stream_rt->stokes_and_jacobian_single_wn
	  (wn(i), Spec_index);
      else
	stokes_only(i, ra) = low_stream_rt->stokes_single_wn
	  (wn(i), Spec_index);
    }

  // Update progress meter in log file, if we are using it.
    if(disp)	
      *disp += 1;
  }

//-----------------------------------------------------------------------
// Now average them, and use to calculate the high and low stream
// values using these averaged values. We then use this to create
// the corrections, and then class to handle interpolating these
// corrections between the gas_opd_avg values.
//
// We actually need two linear interpolations, one that interpolates
// based on gas_opd and one that interpolates based on log(gas_opd).
// We'll use these two separate interpolations when we do the
// correction for each wavenumber. Small gas_opd values get log-linear
// interpolated, and large values get linearly interpolated.
//-----------------------------------------------------------------------

  // We do the RT on the middle wavenumber when calculating the correction.
  double wn_mid = (max(wn) + min(wn)) / 2;
  double wn_mid_closest = wn(0);
  double wn_dist = fabs(wn_mid_closest - wn_mid);
  for(int i = 0; i < wn.rows(); ++i)
    if(fabs(wn(i) - wn_mid) < wn_dist) {
      wn_dist = fabs(wn(i) - wn_mid);
      wn_mid_closest = wn(i);
    }
  linear_interp.clear();
  log_linear_interp.clear();
  for(int i = 0; i < 2; ++i) {
    std::vector<int> cnt_w = cnt[i].value();
    std::vector<AutoDerivative<double> > log_opd_sum_w = log_opd_sum[i].value();
    std::vector<ArrayAd<double, 2> > atm_avg_w = atm_sum[i].value();
    std::vector<AutoDerivative<double> > gas_opd_avg;
    std::vector<Array<AutoDerivative<double>, 1> > err_est;
    for(int j = 0; j < (int) cnt_w.size(); ++j) {
      // We only do a correction for a bin that has data.
      if(cnt_w[j] != 0) {
	gas_opd_avg.push_back(exp(log_opd_sum_w[j] / cnt_w[j]));

	// Get average atmosphere values
	atm_avg_w[j].value() /= cnt_w[j];
	atm_avg_w[j].jacobian() /= cnt_w[j];
	
	// And use to get high and low stream values.
	Array<AutoDerivative<double>, 1> low = 
	  low_stream_rt->stokes_and_jacobian_single_wn
	  (wn_mid_closest, Spec_index, atm_avg_w[j]).to_array();
	Array<AutoDerivative<double>, 1> high = 
	  high_stream_rt->stokes_and_jacobian_single_wn
	  (wn_mid_closest, Spec_index, atm_avg_w[j]).to_array();
	// Use high and low to calculate correction.
	Array<AutoDerivative<double>, 1> c(low.rows());
	// We do a scaled correction for the first stokes coefficient
	// which is the intensity. For the others, we do just a
	// difference. Note sure exactly how much difference this
	// makes, but it is what the original Fortran did. Note that
	// you do *not* want to do a scaled correction for the Q and U
	// stokes parameters because these can be negative so
	// 1/1+err_est can produce infinity if err_est is -1 (or a
	// large number of it is close).
	c(0) = (fabs(high(0).value()) > 1e-15 ? 
		(low(0) - high(0)) / high(0) : 0);
	for(int i =1; i < c.rows(); ++i)
	  c(i) = low(i) - high(i);
	err_est.push_back(c);
      }
    }
    // Now use corrections values to create Linear and Log-linear
    // interpolaters.
    linear_interp.push_back
      (linear_interp_type(gas_opd_avg.begin(), gas_opd_avg.end(), 
			  err_est.begin(), 
			  linear_interp_type::OUT_OF_RANGE_CLIP));
    log_linear_interp.push_back
      (log_linear_interp_type(gas_opd_avg.begin(), gas_opd_avg.end(), 
			      err_est.begin(),
			      linear_interp_type::OUT_OF_RANGE_CLIP));
  }
}

