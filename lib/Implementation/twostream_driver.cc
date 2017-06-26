#include "twostream_driver.h"
#include "old_constant.h"
#include "wgs84_constant.h"
#include "ground.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;


//-----------------------------------------------------------------------
/// Initialize Twostream BRDF interface
//-----------------------------------------------------------------------

TwostreamBrdfDriver::TwostreamBrdfDriver(int surface_type)
{
  int nbeams = 1;
  int n_user_streams = 1;
  int n_user_relazms = 1;

  int nspars;
  switch(surface_type) {
  case LAMBERTIAN:
    nspars = 1;
    break;
  case BREONVEG:
  case BREONSOIL:
    nspars = 5;
    break;
  case COXMUNK:
    nspars = 3;
    break;
  default:
    Exception e("Unhandled BRDF type: ");
    e << surface_type;
    throw e;
  }

  twostream_brdf_.reset( new Twostream_Ls_Brdf_Supplement(nbeams, n_user_streams, n_user_relazms, nspars) );

  // Number of quadtrature streams for BRDF calculation
  // Use value consistent with LIDORT
  twostream_brdf_->nstreams_brdf(50);

  brdf_params.reference( twostream_brdf_->brdf_parameters() );
  brdf_factors.reference( twostream_brdf_->brdf_factors() );
}

void TwostreamBrdfDriver::setup_geometry(double sza, double azm, double zen) const
{
  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> beam_szas( twostream_brdf_->beam_szas() );
  beam_szas(0) = sza;

  // Viewing zenith angle (degrees)
  Array<double, 1> user_angles( twostream_brdf_->user_angles() );
  user_angles(0) = zen;

  // Azimuth angle (degrees)
  Array<double, 1> user_relazms( twostream_brdf_->user_relazms() );
  user_relazms(0) = azm;
}

void TwostreamBrdfDriver::calculate_brdf() const
{
  twostream_brdf_->run();
}

int TwostreamBrdfDriver::n_brdf_kernels() const
{
  return twostream_brdf_->n_brdf_kernels();
}

void TwostreamBrdfDriver::n_brdf_kernels(const int n_kernels)
{
  twostream_brdf_->n_brdf_kernels(n_kernels);
}

int TwostreamBrdfDriver::n_kernel_factor_wfs() const 
{
  return twostream_brdf_->n_kernel_factor_wfs();
}

void TwostreamBrdfDriver::n_kernel_factor_wfs(const int n_factors) 
{
  // Do not really set this, twostream calculates this internally
  // It adds to the passed in value, so book keeping values
  // will occur if this is non-zero
  //twostream_brdf_->n_kernel_factor_wfs(n_factors);
}

int TwostreamBrdfDriver::n_kernel_params_wfs() const 
{
  return twostream_brdf_->n_kernel_params_wfs();
}

void TwostreamBrdfDriver::n_kernel_params_wfs(const int n_params) 
{
  // Do not really set this, twostream calculates this internally
  // It adds to the passed in value, so book keeping values
  // will occur if this is non-zero
  //twostream_brdf_->n_kernel_params_wfs(n_params);
}

int TwostreamBrdfDriver::n_surface_wfs() const 
{
  return twostream_brdf_->n_surface_wfs();
}

void TwostreamBrdfDriver::n_surface_wfs(const int n_wfs) 
{
  // Do not really set this, twostream calculates this internally
  // It adds to the passed in value, so book keeping values
  // will occur if this is non-zero
  //twostream_brdf_->n_surface_wfs(n_wfs);
}

bool TwostreamBrdfDriver::do_kparams_derivs(const int kernel_index) const
{
  return twostream_brdf_->do_kparams_derivs()(kernel_index);
}

void TwostreamBrdfDriver::do_kparams_derivs(const int kernel_index, const bool do_kparams)
{
  // Do not really set this, twostream calculates this internally
  // It adds to the passed in value, so book keeping values
  // will occur if this is non-zero
  //Array<bool, 1> ts_do_kparams_derivs(twostream_brdf_->do_kparams_derivs());
  //ts_do_kparams_derivs(kernel_index) = do_kparams;
  //twostream_brdf_->do_kparams_derivs(ts_do_kparams_derivs);
}

bool TwostreamBrdfDriver::do_shadow_effect() const {
  return twostream_brdf_->do_shadow_effect();
}

void TwostreamBrdfDriver::do_shadow_effect(const bool do_shadow) const {
  twostream_brdf_->do_shadow_effect(do_shadow);
}

void TwostreamBrdfDriver::initialize_kernel_parameters(const int kernel_index,
						       const int which_brdf,
						       const bool lambertian_flag,
						       const int n_brdf_parameters,
						       const bool do_factor_wfs,
						       const blitz::Array<bool, 1>& do_params_wfs)
{
  Array<int, 1> ts_which_brdf( twostream_brdf_->which_brdf() );
  ts_which_brdf(kernel_index) = which_brdf;

  Array<bool, 1> ts_lambertian_flag = twostream_brdf_->lambertian_kernel_flag();
  ts_lambertian_flag(kernel_index) = lambertian_flag;
  twostream_brdf_->lambertian_kernel_flag(ts_lambertian_flag);

  Array<int, 1> ts_n_brdf_parameters( twostream_brdf_->n_brdf_parameters() );
  ts_n_brdf_parameters(kernel_index) = n_brdf_parameters;

  Array<bool, 1> ts_do_factor_wfs( twostream_brdf_->do_kernel_factor_wfs() );
  ts_do_factor_wfs(kernel_index) = do_factor_wfs;
  twostream_brdf_->do_kernel_factor_wfs(ts_do_factor_wfs);

  Array<bool, 2> ts_do_params_wfs( twostream_brdf_->do_kernel_params_wfs() );
  ts_do_params_wfs(kernel_index, Range(0, do_params_wfs.rows()-1)) = do_params_wfs;
  twostream_brdf_->do_kernel_params_wfs(ts_do_params_wfs);
}

//=======================================================================
/// TwostreamRtDriver
/// Sizes of layers, and number of jacobians must be set up in construtor
/// as seen in signature.
/// Use do_fullquadratures = false only for comparison against LIDORT 
//=======================================================================

TwostreamRtDriver::TwostreamRtDriver(int nlayers, int npars, int surface_type, bool do_fullquadrature, bool pure_nadir)
  : surface_type_(surface_type), do_fullquadrature_(do_fullquadrature), pure_nadir_(pure_nadir)
{
  brdf_driver_.reset( new TwostreamBrdfDriver(surface_type_) );

  int nthreads       = 1;
  int nbeams         = 1;
  int n_user_angles  = 1;
  int n_user_relazms = 1;

  int nspars         = brdf_interface()->nspars();
  int n_geometries   = nbeams * n_user_angles * n_user_relazms;
  int ntotal         = 2 * nlayers;

  int thread = 1; // 1 because this number is a Fortran index
  double earth_radius = OldConstant::wgs84_a.convert(units::km).value;

  twostream_interface_.reset( new Twostream_L_Master( thread, nthreads, nlayers, ntotal, n_geometries, n_user_angles, n_user_relazms, nbeams, earth_radius, npars, nspars) );

  // Lambertian albedo values are stored seperate from BRDF data structures
  // Do this after TwoStream has been instantiated
  brdf_driver_->set_lambertian_albedo( twostream_interface_->lambertian_albedo() );

  // Initialize BRDF data structure
  brdf_driver_->initialize_brdf_inputs(surface_type_);

  initialize_rt();
}

void TwostreamRtDriver::initialize_rt()
{
  // We are doing solar sources only
  twostream_interface_->do_solar_sources(true);

  // Two choices of stream value................ CHOOSE One !!!!!
  if (do_fullquadrature_)
    twostream_interface_->stream_value( sqrt(1.0e0 / 3.0e0 ) );
  else
    twostream_interface_->stream_value( 0.5e0 );

  // Set pure nadir flag
  twostream_interface_->pure_nadir(pure_nadir_);

  // Set stream value for BRDF interface to the same value used by 2stream interface
  brdf_interface()->stream_value( twostream_interface_->stream_value() );

  // Leave false if doing lambertian
  if(surface_type_ == LAMBERTIAN)
    twostream_interface_->do_brdf_surface(false);
  else
    twostream_interface_->do_brdf_surface(true);
   
  // Flags for viewing mode
  twostream_interface_->do_upwelling(true);
  twostream_interface_->do_dnwelling(false);

  // In most instances this flag for delta-m scaling should be set
  twostream_interface_->do_d2s_scaling(true);

  // Beam source flux, same value used for all solar angles
  // Normally set to 1 for "sun-normalized" output
  twostream_interface_->flux_factor(1.0);

  // Flag for calculating profile Jacobians in layer n
  // Calculate jacobians for all layers
  Array<bool, 1> layer_jac_flag(twostream_interface_->layer_vary_flag());
  layer_jac_flag = true;
  twostream_interface_->layer_vary_flag(layer_jac_flag);

  // Link twostream interface brdf inputs to brdf interface
  twostream_interface_->brdf_f_0().reference(brdf_interface()->brdf_f_0());
  twostream_interface_->brdf_f().reference(brdf_interface()->brdf_f());
  twostream_interface_->ubrdf_f().reference(brdf_interface()->ubrdf_f());

  twostream_interface_->ls_brdf_f_0().reference(brdf_interface()->ls_brdf_f_0());
  twostream_interface_->ls_brdf_f().reference(brdf_interface()->ls_brdf_f());
  twostream_interface_->ls_ubrdf_f().reference(brdf_interface()->ls_ubrdf_f());
}

void TwostreamRtDriver::setup_height_grid(const blitz::Array<double, 1>& in_height_grid) const
{

  int nlayers = in_height_grid.rows() - 1;

  // Make sure incoming height grid is compatible with initialized size
  if (twostream_interface_->nlayers() != nlayers)
    throw Exception("TwostreamDriver does not support the number of layers changing during the retrieval.");

  // Set only the layers being used
  Array<double, 1> ts_height_grid( twostream_interface_->height_grid() );
  ts_height_grid(Range(0, nlayers)) = in_height_grid;
}

void TwostreamRtDriver::setup_geometry(double sza, double azm, double zen) const
{
  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> ts_sza( twostream_interface_->beam_szas() );
  ts_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> ts_azm( twostream_interface_->user_relazms() );
  ts_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> ts_zen( twostream_interface_->user_angles() );
  ts_zen(0) = zen;
}

void TwostreamRtDriver::setup_optical_inputs(const blitz::Array<double, 1>& od, 
					     const blitz::Array<double, 1>& ssa,
					     const blitz::Array<double, 2>& pf) const
{
  // Ranges for copying inputs to method
  Range rlay(0, od.extent(firstDim) - 1);

  // Vertical optical depth thicness values for all layers and threads
  Array<double, 2> deltau( twostream_interface_->deltau_input() );
  deltau(rlay, 0) = od;

  // Single scattering albedos for all layers and threads 
  Array<double, 2> omega( twostream_interface_->omega_input() );
  omega(rlay, 0) = where(ssa > 0.999, 0.999999, ssa);

  // Assymetry factor
  // Equal to one third of the first phase function moment
  Array<double, 2> asymm_input( twostream_interface_->asymm_input() );
  asymm_input(rlay, 0) = pf(1, rlay) / 3.0;

  // Delta-m scaling factor for 2-stream
  // Equal to one fifth of the second  phase function moment
  if (twostream_interface_->do_d2s_scaling()) {
    Array<double, 2> d2s_scaling( twostream_interface_->d2s_scaling() );
    d2s_scaling(rlay, 0) = pf(2, rlay) / 5.0;
  }
}

void TwostreamRtDriver::clear_linear_inputs() const
{
  // Disable wfs for both types supported
  twostream_interface_->do_profile_wfs(false);
  twostream_interface_->do_surface_wfs(false);
}

void TwostreamRtDriver::setup_linear_inputs(const ArrayAd<double, 1>& od, 
                                            const ArrayAd<double, 1>& ssa,
                                            const ArrayAd<double, 2>& pf,
                                            bool do_surface_linearization) const
{
  // Number of profile weighting functions in layer n
  int natm_jac = od.number_variable();
  Array<int, 1> layer_jac_number( twostream_interface_->layer_vary_number() );
  layer_jac_number = natm_jac;

  // Ranges for copying inputs to method
  Range rlay(0, od.rows() - 1);         // number of layers
  Range rjac(0, natm_jac - 1);
  Range all(Range::all());

  // Flag for output of profile Jacobians
  // Unlikely we won't need these
  twostream_interface_->do_profile_wfs(true);

  // Needs to match the number set in the BRDF structure
  twostream_interface_->n_surface_wfs(brdf_driver_->n_surface_wfs());

  // Flag for output of surface Jacobians
  // Certainly wouldn't need this if not retrieving ground
  twostream_interface_->do_surface_wfs(do_surface_linearization);

  // Check that we fit within the LIDORT configuration
  if(twostream_interface_->l_deltau_input().rows() < od.rows()) {
    Exception e;
    e << "The number of layers you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << od.rows() << "\nand the maximum allowed is "
      << twostream_interface_->l_deltau_input().rows() << "\n"
      << "\n"
      << "You might try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxlayer=value set to a larger value.\n";
    throw e;
  }
  if(twostream_interface_->l_deltau_input().cols() < natm_jac) {
    Exception e;
    e << "The number of jacobians you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << natm_jac << "\nand the maximum allowed is "
      << twostream_interface_->l_deltau_input().cols() << "\n"
      << "\n"
      << "This number of jacobians is a function of the number of aerosols\n"
      << "in your state vector, so you can reduce the number of aerosols\n"
      << "\n"
      << "You might also try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxatmoswfs=value set to a larger value.\n";
    throw e;
  }

  // Setup optical linear inputs
  Array<double, 2> l_deltau( twostream_interface_->l_deltau_input()(rlay,rjac,0) );
  Array<double, 2> l_omega( twostream_interface_->l_omega_input()(rlay,rjac,0) );
  Array<double, 2> l_asymm_input( twostream_interface_->l_asymm_input()(rlay,rjac,0) );
  Array<double, 2> l_d2s_scaling( twostream_interface_->l_d2s_scaling()(rlay,rjac,0) );

  l_deltau = od.jacobian();
  l_omega  = ssa.jacobian();
  
  l_asymm_input = pf.jacobian()(1, all, all) / 3.0;
  l_d2s_scaling = pf.jacobian()(2, all, all) / 5.0;
}

void TwostreamRtDriver::calculate_rt() const
{
  twostream_interface_->run();

  //  Exception handling
  if (twostream_interface_->status_inputcheck() == 1) {
    stringstream err_msg;
    err_msg << "TwostreamInterace input check failed:" << std::endl;
    for(int k = 0; k < twostream_interface_->c_nmessages(); k++) {
      err_msg << " - Message # " << k << ": " << twostream_interface_->c_messages()[k] << std::endl
              << " - Action  # " << k << ": " << twostream_interface_->c_actions()[k] << std::endl;
    }
    throw Exception(err_msg.str());
  }
  if (twostream_interface_->status_execution() == 1) {
    stringstream err_msg;
    err_msg << "TwostreamInterace execution failed:" << std::endl
            << twostream_interface_->e_message() << std::endl
            << twostream_interface_->e_trace_1() << std::endl
            << twostream_interface_->e_trace_2() << std::endl;
    throw Exception(err_msg.str());
  }

}

double TwostreamRtDriver::get_intensity() const
{
  return twostream_interface_->intensity_toa()(0, 0);
}

void TwostreamRtDriver::copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf) const
{
  // Copy out jacobian values
  Range ra(Range::all());
  jac_atm.reference( twostream_interface_->profilewf_toa()(0, ra, ra, 0).copy() );
  jac_atm.transposeSelf(secondDim, firstDim); // swap to same ordering as lidort
    
  jac_surf.reference( twostream_interface_->surfacewf_toa()(0, ra, 0).copy() );
}

