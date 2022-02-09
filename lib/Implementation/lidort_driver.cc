#include "lidort_driver.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "ground.h"
#include "ostream_pad.h"
#include "fe_disable_exception.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

//=======================================================================
// LidortBrdfInterface
//=======================================================================

//-----------------------------------------------------------------------
/// Initialize Lidort BRDF interface
//-----------------------------------------------------------------------

LidortBrdfDriver::LidortBrdfDriver(int nstream, int nmoment) : nmoment_(nmoment)
{
  brdf_interface_.reset( new Brdf_Lin_Sup_Masters() );

  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();

  brdf_inputs.bs_do_solar_sources(true);

  // Only use 1 beam meaning only one set of sza, azm
  brdf_inputs.bs_nbeams(1);
  brdf_inputs.bs_n_user_streams(1);
  brdf_inputs.bs_n_user_relazms(1);

  // This MUST be consistent with streams used for
  // LIDORT RT calculation
  brdf_inputs.bs_nstreams(nstream);

  // Recommended value from LIDORT manual
  // Number of quadtrature streams for BRDF calculation
  brdf_inputs.bs_nstreams_brdf(50);

  brdf_params.reference( brdf_inputs.bs_brdf_parameters() );
  brdf_factors.reference( brdf_inputs.bs_brdf_factors() );
}

void LidortBrdfDriver::setup_geometry(double sza, double azm, double zen) const
{

  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();

  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> brdf_sza( brdf_inputs.bs_beam_szas() );
  brdf_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> brdf_azm( brdf_inputs.bs_user_relazms() );
  brdf_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> brdf_zen( brdf_inputs.bs_user_angles_input() );
  brdf_zen(0) = zen;
}

void LidortBrdfDriver::calculate_brdf() const
{
  // Lidort may cause floating point exceptions when doing setup. This
  // is because it may copy garbage value, which are never used. By
  // chance the garbage values may cause a overflow. We
  // suspend floating point exceptions when doing setup
  FeDisableException disable_fp;

  // Process BRDF inputs
  bool do_debug_restoration = false;
  brdf_interface_->run(do_debug_restoration, nmoment_);
}

int LidortBrdfDriver::n_brdf_kernels() const
{
  return brdf_interface_->brdf_sup_in().bs_n_brdf_kernels();
}

void LidortBrdfDriver::n_brdf_kernels(const int n_kernels)
{
  brdf_interface_->brdf_sup_in().bs_n_brdf_kernels(n_kernels);
}

int LidortBrdfDriver::n_kernel_factor_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_kernel_factor_wfs();
}

void LidortBrdfDriver::n_kernel_factor_wfs(const int n_factors) {
  brdf_interface_->brdf_linsup_in().bs_n_kernel_factor_wfs(n_factors);
}

int LidortBrdfDriver::n_kernel_params_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_kernel_params_wfs();
}

void LidortBrdfDriver::n_kernel_params_wfs(const int n_params) {
  brdf_interface_->brdf_linsup_in().bs_n_kernel_params_wfs(n_params);
}

int LidortBrdfDriver::n_surface_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_surface_wfs();
}

void LidortBrdfDriver::n_surface_wfs(const int n_wfs) {
  brdf_interface_->brdf_linsup_in().bs_n_surface_wfs(n_wfs);
}

bool LidortBrdfDriver::do_kparams_derivs(const int kernel_index) const
{
  return brdf_interface_->brdf_linsup_in().bs_do_kparams_derivs()(kernel_index);
}

void LidortBrdfDriver::do_kparams_derivs(const int kernel_index, const bool do_kparams)
{
  Brdf_Linsup_Inputs& brdf_lin_inputs = brdf_interface_->brdf_linsup_in();
  Array<bool, 1> do_kparams_derivs( brdf_lin_inputs.bs_do_kparams_derivs() );
  do_kparams_derivs(kernel_index) = do_kparams;
  brdf_lin_inputs.bs_do_kparams_derivs(do_kparams_derivs);
}

bool LidortBrdfDriver::do_shadow_effect() const {
  return brdf_interface_->brdf_sup_in().bs_do_shadow_effect();
}

void LidortBrdfDriver::do_shadow_effect(const bool do_shadow) const {
  brdf_interface_->brdf_sup_in().bs_do_shadow_effect(do_shadow);
}

void LidortBrdfDriver::initialize_kernel_parameters(const int kernel_index,
                                                    const int which_brdf,
                                                    const bool lambertian_flag,
                                                    const int n_brdf_parameters,
                                                    const bool do_factor_wfs,
                                                    const blitz::Array<bool, 1>& do_params_wfs)
{
  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();
  Brdf_Linsup_Inputs& brdf_lin_inputs = brdf_interface_->brdf_linsup_in();

  Array<int, 1> bs_which_brdf( brdf_inputs.bs_which_brdf() );
  bs_which_brdf(kernel_index) = which_brdf;

  Array<bool, 1> bs_lambertian_flag = brdf_inputs.bs_lambertian_kernel_flag();
  bs_lambertian_flag(kernel_index) = lambertian_flag;
  brdf_inputs.bs_lambertian_kernel_flag(bs_lambertian_flag);

  Array<int, 1> bs_n_brdf_parameters( brdf_inputs.bs_n_brdf_parameters() );
  bs_n_brdf_parameters(kernel_index) = n_brdf_parameters;

  Array<bool, 1> bs_do_factor_wfs( brdf_lin_inputs.bs_do_kernel_factor_wfs() );
  bs_do_factor_wfs(kernel_index) = do_factor_wfs;
  brdf_lin_inputs.bs_do_kernel_factor_wfs(bs_do_factor_wfs);

  Array<bool, 2> bs_do_params_wfs( brdf_lin_inputs.bs_do_kernel_params_wfs() );
  bs_do_params_wfs(kernel_index, Range(0, do_params_wfs.rows()-1)) = do_params_wfs;
  brdf_lin_inputs.bs_do_kernel_params_wfs(bs_do_params_wfs);
}

//=======================================================================
// LidortRtDriver
//=======================================================================

LidortRtDriver::LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type, const blitz::Array<double, 1>& zen, bool pure_nadir)
  : nstream_(nstream), nmoment_(nmoment), do_multi_scatt_only_(do_multi_scatt_only), surface_type_(surface_type), pure_nadir_(pure_nadir)
{
  brdf_driver_.reset( new LidortBrdfDriver(nstream, nmoment) );
  lidort_interface_.reset( new Lidort_Lps_Masters() );

  // Check inputs against sizes allowed by LIDORT
  Lidort_Pars lid_pars = Lidort_Pars::instance();
  range_check(nstream, 1, lid_pars.maxstreams+1);
  range_check(nmoment, 3, lid_pars.maxmoments_input+1);

  // Initialize BRDF data structure
  brdf_driver()->initialize_brdf_inputs(surface_type_);

  initialize_rt();

  // Set up scatting mode based on viewing zenith angle
  setup_sphericity(max(zen));
}

int LidortRtDriver::number_moment() const
{
  if (!lidort_interface_) {
    throw Exception("Lidort Interface not initialized");
  }
  return lidort_interface_->lidort_modin().mcont().ts_nmoments_input();
}

int LidortRtDriver::number_stream() const
{
  if (!lidort_interface_) {
    throw Exception("Lidort Interface not initialized");
  }
  return lidort_interface_->lidort_fixin().cont().ts_nstreams();
}

void LidortRtDriver::initialize_rt()
{
  // Set up these references for convienence
  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();
  Lidort_Modified_Control& mcontrol_inputs = lidort_interface_->lidort_modin().mcont();
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();

  Lidort_Fixed_Control& fcontrol_inputs = lidort_interface_->lidort_fixin().cont();
  Lidort_Fixed_Sunrays& fbeam_inputs = lidort_interface_->lidort_fixin().sunrays();
  Lidort_Modified_Sunrays& mbeam_inputs = lidort_interface_->lidort_modin().msunrays();
  Lidort_Modified_Chapman& mchapman_inputs = lidort_interface_->lidort_modin().mchapman();

  Lidort_Modified_Uservalues& muser_inputs = lidort_interface_->lidort_modin().muserval();
  Lidort_Fixed_Uservalues& fuser_inputs = lidort_interface_->lidort_fixin().userval();

  Lidort_Fixed_Lincontrol& lincontrol = lidort_interface_->lidort_linfixin().cont();

  // Set values used for all calculations
  // Number of quadtature streams in the cosine half space
  fcontrol_inputs.ts_nstreams(nstream_);

  // Number of Legendre expansion coefficients for the phase function
  mcontrol_inputs.ts_nmoments_input(nmoment_);

  // Needed for atmospheric scattering of sunlight
  mboolean_inputs.ts_do_solar_sources(true);

  fboolean_inputs.ts_do_thermal_emission(false);
  fboolean_inputs.ts_do_surface_emission(false);

  // Always use BRDF supplement, don't use specialized lambertian_albedo mode
  brdf_interface()->brdf_sup_in().bs_do_brdf_surface(true);
  fboolean_inputs.ts_do_brdf_surface(true);

  // Flags for viewing mode
  fboolean_inputs.ts_do_upwelling(true);
  fboolean_inputs.ts_do_dnwelling(false);

  // Fourier azimuth series is examined twice for convergence
  mboolean_inputs.ts_do_double_convtest(true);

  // Do internal calculation of slanth path optical depths
  mboolean_inputs.ts_do_chapman_function(true);

  // In most instances this flag for delta-m scaling should be set
  mboolean_inputs.ts_do_deltam_scaling(true);

  // There will be output at a number of off-quadrature
  // zenith angles specified by the user, this is the normal case
  mboolean_inputs.ts_do_user_streams(true);
  brdf_interface()->brdf_sup_in().bs_do_user_streams(true);

  // Accuracy criterion for convergence of Fourier series in
  // relative azimuth. Set to recommended value.
  fcontrol_inputs.ts_lidort_accuracy(1e-8);

  // New value in LIDORT 3.8.3 that if left at zero will cause floating point errors in asymtx
  fcontrol_inputs.ts_asymtx_tolerance(1e-20);

  // Beam source flux, same value used for all solar angles
  // Normally set to 1 for "sun-normalized" output
  fbeam_inputs.ts_flux_factor(1.0);

  // Number of solar beams
  mbeam_inputs.ts_nbeams(1);

  // Earth's radius in km
  mchapman_inputs.ts_earth_radius( OldConstant::wgs84_a.
                                  convert(units::km).value);

  // Number of user-defined relative azimuth angles
  muser_inputs.ts_n_user_relazms(1);

  // Number of user-defined viewing zenith angles
  muser_inputs.ts_n_user_streams(1);

  // Number of vertical output levels
  fuser_inputs.ts_n_user_levels(1);

  // Flag for calculating profile Jacobians in layer n
  // Calculate jacobians for all layers
  Array<bool, 1> layer_jac_flag(lincontrol.ts_layer_vary_flag());
  layer_jac_flag = true;
  lincontrol.ts_layer_vary_flag(layer_jac_flag);

  // Needs to match the number set in the BRDF structure
  lincontrol.ts_n_surface_wfs(brdf_interface()->brdf_linsup_in().bs_n_surface_wfs());
}

/// Set up/reset sphericity mode which may be affected by
/// the current zenith viewing angle
void LidortRtDriver::setup_sphericity(double zen) const
{

  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  Lidort_Fixed_Control& fcontrol_inputs = lidort_interface_->lidort_fixin().cont();

  if(do_multi_scatt_only_) {
    // Do not do full reflectance calculation if only multiple scattering needed
    // These flags should already be false, but just in case..

    // Only return diffuse multiple scattering component
    fboolean_inputs.ts_do_fullrad_mode(false);

    fboolean_inputs.ts_do_plane_parallel(false);

    // Also SS corr and plane-parallel flags are false
    mboolean_inputs.ts_do_focorr(false);
    mboolean_inputs.ts_do_focorr_nadir(false);
    mboolean_inputs.ts_do_focorr_outgoing(false);

  } else { // if not do multi_scatt_only
    // Do full SS + MS calculation and use LOS correction

    // Do a full reflectance calculation
    fboolean_inputs.ts_do_fullrad_mode(true);

    fboolean_inputs.ts_do_plane_parallel(false);

    // Pseudo-spherical + Line of Sight correction
    mboolean_inputs.ts_do_focorr(true);
    mboolean_inputs.ts_do_focorr_nadir(false);
    mboolean_inputs.ts_do_focorr_outgoing(true);

    // Number of fine layers subdividing coarse layering
    // Only used during LOS correction
    fcontrol_inputs.ts_nfinelayers(4);
  }

  // Flag for controlling azimuth dependence in the output
  // LIDORT will complain if user zenith is 0 and this is not
  // set, when not using ss correction mode
  if( pure_nadir_ ) {
    if (mboolean_inputs.ts_do_focorr_outgoing()) {
      // Use Pseudo-spherical correction instead for these small viewing zenith angles
      mboolean_inputs.ts_do_focorr_outgoing(false);
      mboolean_inputs.ts_do_focorr_nadir(true);
    }
  }
}

/// Set plane parallel sphericity
void LidortRtDriver::set_plane_parallel() const
{
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();

  mboolean_inputs.ts_do_focorr(false);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(false);
  fboolean_inputs.ts_do_plane_parallel(true);
  mboolean_inputs.ts_do_no_azimuth(true);
}

/// Set pseudo spherical sphericity
void LidortRtDriver::set_pseudo_spherical() const
{
  // Lidort may cause floating point exceptions when doing setup. This
  // is because it may copy garbage value, which are never used. By
  // chance the garbage values may cause a overflow. We
  // suspend floating point exceptions when doing setup
  FeDisableException disable_fp;

  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(true);
  fboolean_inputs.ts_do_plane_parallel(false);
  mboolean_inputs.ts_do_no_azimuth(false);
}

/// Set plane parallel plus single scattering correction
void LidortRtDriver::set_plane_parallel_plus_ss_correction() const
{
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(false);
  mboolean_inputs.ts_do_focorr_nadir(true);
  fboolean_inputs.ts_do_plane_parallel(true);
  mboolean_inputs.ts_do_no_azimuth(false);
}

/// Set line of sight mode
void LidortRtDriver::set_line_of_sight() const
{
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  Lidort_Fixed_Boolean& fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();

  mboolean_inputs.ts_do_focorr(true);
  mboolean_inputs.ts_do_focorr_outgoing(true);
  mboolean_inputs.ts_do_focorr_nadir(false);
  fboolean_inputs.ts_do_plane_parallel(false);
  mboolean_inputs.ts_do_no_azimuth(false);
}

void LidortRtDriver::setup_height_grid(const blitz::Array<double, 1>& in_height_grid) const
{

  Lidort_Fixed_Chapman& fchapman_inputs = lidort_interface_->lidort_fixin().chapman();
  Lidort_Modified_Uservalues& muser_inputs = lidort_interface_->lidort_modin().muserval();

  Array<double, 1> lidort_height_grid( fchapman_inputs.ts_height_grid() );
  int nlayer = in_height_grid.extent(firstDim) - 1;
  Range lay_range = Range(0,nlayer);
  lidort_height_grid(lay_range) = in_height_grid;

  // Set GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
  muser_inputs.ts_geometry_specheight( lidort_height_grid(nlayer) );

  // Tell LIDORT number of layers
  lidort_interface_->lidort_fixin().cont().ts_nlayers(nlayer);
}

void LidortRtDriver::setup_geometry(double sza, double azm, double zen) const
{

  Lidort_Modified_Sunrays& mbeam_inputs = lidort_interface_->lidort_modin().msunrays();
  Lidort_Modified_Uservalues& muser_inputs = lidort_interface_->lidort_modin().muserval();

  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> ld_sza( mbeam_inputs.ts_beam_szas() );
  ld_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> ld_azm( muser_inputs.ts_user_relazms() );
  ld_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> ld_zen( muser_inputs.ts_user_angles_input() );
  ld_zen(0) = zen;
}

void LidortRtDriver::setup_optical_inputs(const blitz::Array<double, 1>& od,
                                          const blitz::Array<double, 1>& ssa,
                                          const blitz::Array<double, 2>& pf) const
{
  // Ranges for copying inputs to method
  Range rlay(0, od.extent(firstDim) - 1);
  Range rmom(0, pf.extent(firstDim) - 1);

  // Convienence references
  Lidort_Fixed_Optical& foptical_inputs = lidort_interface_->lidort_fixin().optical();
  Lidort_Modified_Optical& moptical_inputs = lidort_interface_->lidort_modin().moptical();

  // Vertical optical depth thicness values for all layers and threads
  Array<double, 1> deltau( foptical_inputs.ts_deltau_vert_input() );
  deltau(rlay) = od;

  // Single scattering albedos for all layers and threads
  Array<double, 1> omega( moptical_inputs.ts_omega_total_input() );
  omega(rlay) = where(ssa > 0.999, 0.999999, ssa);

  // For all layers n and threads t, Legrenre moments of
  // the phase function expansion multiplied by (2L+1);
  // initial value (L=0) should always be 1
  // phasmoms_total_input(n, L, t)
  // n = moments, L = layers, t = threads
  Array<double, 2> phasmoms( foptical_inputs.ts_phasmoms_total_input() );
  phasmoms(rmom, rlay) = where(abs(pf) > 1e-11, pf, 1e-11);
}

void LidortRtDriver::clear_linear_inputs() const
{
  Lidort_Modified_Lincontrol& lincontrol = lidort_interface_->lidort_linmodin().mcont();

  // Flag for output of profile Jacobians
  lincontrol.ts_do_profile_linearization(false);

  // Flag for output of surface Jacobians
  lincontrol.ts_do_surface_linearization(false);
}

void LidortRtDriver::setup_linear_inputs(const ArrayAd<double, 1>& od,
                                         const ArrayAd<double, 1>& ssa,
                                         const ArrayAd<double, 2>& pf,
                                         bool do_surface_linearization) const
{

  if(od.number_variable() > Lidort_Pars::instance().max_atmoswfs) {
    Exception err;
    err << "LIDORT has been compiled to allow a maximum of " << Lidort_Pars::instance().max_atmoswfs
        << " atmosphere derivatives to be calculated. We are trying to calculate "
        << od.number_variable() << " atmosphere derivatives";
    throw err;
  }

  Lidort_Fixed_Lincontrol& lincontrol = lidort_interface_->lidort_linfixin().cont();
  Lidort_Modified_Lincontrol& mlincontrol = lidort_interface_->lidort_linmodin().mcont();
  Lidort_Fixed_Linoptical& linoptical = lidort_interface_->lidort_linfixin().optical();

  // Number of profile weighting functions in layer n
  int natm_jac = od.number_variable();
  Array<int, 1> layer_jac_number( lincontrol.ts_layer_vary_number() );
  layer_jac_number = natm_jac;

  // Ranges for copying inputs to method
  Range ra(Range::all());
  Range rlay(0, od.rows() - 1);         // number of layers
  Range rjac(0, natm_jac - 1);

  Range rmom(0, pf.rows() - 1); // number phase function moments
  Range all(Range::all());

  // Flag for output of profile Jacobians
  // Unlikely we won't need these
  mlincontrol.ts_do_profile_linearization(true);

  // Flag for output of surface Jacobians
  // Certainly wouldn't need this if not retrieving ground
  mlincontrol.ts_do_surface_linearization(do_surface_linearization);

    // Check that we fit within the LIDORT configuration
  if(linoptical.ts_l_deltau_vert_input().cols() < od.rows()) {
    Exception e;
    e << "The number of layers you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << od.rows() << "\nand the maximum allowed is "
      << linoptical.ts_l_deltau_vert_input().cols() << "\n"
      << "\n"
      << "You might try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxlayer=value set to a larger value.\n";
    throw e;
  }
  if(linoptical.ts_l_deltau_vert_input().rows() < natm_jac) {
    Exception e;
    e << "The number of jacobians you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << natm_jac << "\nand the maximum allowed is "
      << linoptical.ts_l_deltau_vert_input().rows() << "\n"
      << "\n"
      << "This number of jacobians is a function of the number of aerosols\n"
      << "in your state vector, so you can reduce the number of aerosols\n"
      << "\n"
      << "You might also try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxatmoswfs=value set to a larger value.\n";
    throw e;
  }

  // Setup optical linear inputs
  // LIDORT expects the following:
  // l_deltau = xi/tau * dtau/dxi
  // l_omega = xi/omega * domega/dxi
  // ... etc
  // The driver is handed dtau/dxi, domega/dxi ... etc
  // Where it is normally set up for xi being taug, taur, taua[0], taua[1]...
  //
  // LIDORT returns to us:
  // xi * dI/dxi
  // Therefore you will notice that by not multiplying dtau/dxi by xi and only
  // dividing by xi, we are cancelling out the xi in the result and hence
  // the driver really return dI/dxi
  Array<double, 2> l_deltau( linoptical.ts_l_deltau_vert_input()(rjac,rlay) );
  Array<double, 2> l_omega( linoptical.ts_l_omega_total_input()(rjac,rlay) );
  Array<double, 3> l_phasmoms( linoptical.ts_l_phasmoms_total_input()(rjac,rmom,rlay) );

  // Transpose these to match dimensions used internally
  l_deltau.transposeSelf(secondDim, firstDim);
  l_omega.transposeSelf(secondDim, firstDim);

  // Note that LIDORT does *not* take l_deltau etc. Rather it takes
  // what it calls "normalized derivative inputs". So for an optical
  // quantity like taug wrt to tau, this would be taug / tau *
  // dtau/dtaug.
  //
  // It then returns a normalized jacobian, which for taug would be
  // taug * d I /dtaug.
  //
  // Note that we actually leave out one of these factors, because it
  // effectively cancels out. We pass in 1 / tau * dtau / dtaug and
  // get back d I / dtaug.

  firstIndex i1; secondIndex i2;
  l_deltau = where(od.value()(i1) != 0, od.jacobian() / od.value()(i1), 0.0);

  Array<double, 1> ssa_limit(ssa.rows());
  ssa_limit = where(ssa.value() > 0.999, 0.999999, ssa.value());
  l_omega  = where(ssa_limit(i1) != 0, ssa.jacobian() / ssa_limit(i1), 0.0);

  if(pf.is_constant())
    l_phasmoms(rjac, rmom, rlay) = 0.0;
  else {
    blitz::Array<double, 2> pf_in( where(abs(pf.value()) > 1e-11, pf.value(), 1e-11) );

    // We need this loop since l_phasmoms and pf variables have jacobian data in different dimensions
    for (int jidx = 0; jidx < pf.number_variable(); jidx++)
      l_phasmoms(jidx, rmom, rlay) = pf.jacobian()(rmom, rlay, jidx) / pf_in(rmom, rlay)(i1,i2);
  }
}

/// Copy outputs from BRDF supplement into LIDORT Sup inputs types
void LidortRtDriver::copy_brdf_sup_outputs() const {

  // Copy BRDF outputs to LIDORT's BRDF inputs
  Brdf_Sup_Outputs& brdf_outputs = brdf_interface()->brdf_sup_out();
  Brdf_Linsup_Outputs& brdf_lin_outputs = brdf_interface()->brdf_linsup_out();

  Lidort_Sup_Brdf& lid_brdf = lidort_interface_->lidort_sup().brdf();
  Lidort_Linsup_Brdf& lid_lin_brdf = lidort_interface_->lidort_linsup().brdf();

  lid_brdf.copy_from_sup( brdf_outputs );
  lid_lin_brdf.copy_from_sup( brdf_lin_outputs );
}

void LidortRtDriver::calculate_rt() const
{
  // Must copy current BRDF supplement outputs into datastructures used by LIDORT
  copy_brdf_sup_outputs();

  // Call LIDORT for calculations
  lidort_interface_->run(false);
}

double LidortRtDriver::get_intensity() const
{
  // So we know index of intensity
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  // Total Intensity I(t,v,d,T) at output level t, output geometry v,
  // direction d
  return lidort_interface_->lidort_out().main().ts_intensity()(0,0,lid_pars.upidx-1);
}

void LidortRtDriver::copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf) const
{
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  Lidort_Linatmos& lpoutputs = lidort_interface_->lidort_linout().atmos();
  Lidort_Linsurf& lsoutputs = lidort_interface_->lidort_linout().surf();

  Range ra(Range::all());

  // Surface Jacobians KR(r,t,v,d) with respect to surface variable r
  // at output level t, geometry v, direction d
  jac_surf.reference( lsoutputs.ts_surfacewf()(ra, 0, 0, lid_pars.upidx-1).copy() );

  // Get profile jacobians
  // Jacobians K(q,n,t,v,d) with respect to profile atmospheric variable
  // q in layer n, at output level t, geometry v, direction d
  jac_atm.reference( lpoutputs.ts_profilewf()(ra, ra, 0, 0, lid_pars.upidx-1).copy() );
}
