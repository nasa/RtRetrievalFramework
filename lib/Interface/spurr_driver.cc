#include "spurr_driver.h"
#include "spurr_brdf_types.h"
#include "ground.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Initializes the BRDF kernels for the given Ground surface type 
/// integer. A surface type might consist of multiple kernels.
//----------------------------------------------------------------------

void SpurrBrdfDriver::initialize_brdf_inputs(int surface_type) 
{
  switch (surface_type) {
  case LAMBERTIAN:
    initialize_brdf_kernel(LAMBERTIAN);
    break;
  case COXMUNK:
    initialize_brdf_kernel(COXMUNK);
    initialize_brdf_kernel(LAMBERTIAN);
    break;
  case BREONVEG:
  case BREONSOIL:
    initialize_brdf_kernel(RAHMAN);
    initialize_brdf_kernel(surface_type);
    break;
  default:
    Exception e("Unhandled BRDF type index: ");
    e << surface_type;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Initializes a specific BRDF kernel based on the kernel type integer
/// Each call adds a new kernel setup to the BRDF interface, meaning
/// multiple calls are set up by successive calls to this method.
//----------------------------------------------------------------------

void SpurrBrdfDriver::initialize_brdf_kernel(int which_brdf) {
  // Increment kernel index
  int kernel_index = n_brdf_kernels();
  n_brdf_kernels(n_brdf_kernels()+1);

  bool lambertian_flag = false;
  int n_brdf_parameters;
  bool do_factor_wfs;
  Array<bool, 1> do_params_wfs;

  switch (which_brdf) {
  case LAMBERTIAN:
    lambertian_flag = true;
    n_brdf_parameters = 1;

    do_factor_wfs = true;

    do_params_wfs.resize(n_brdf_parameters);
    do_params_wfs(0) = false;
    break;
  case COXMUNK:
    n_brdf_parameters = 3;

    do_factor_wfs = true;

    do_params_wfs.resize(n_brdf_parameters);
    do_params_wfs(0) = true;
    do_params_wfs(1) = true;
    do_params_wfs(2) = false;
    break;
  case RAHMAN:
    n_brdf_parameters = 3;

    do_factor_wfs = true;

    do_params_wfs.resize(n_brdf_parameters);
    do_params_wfs = true;
    break;
  case BREONVEG:
  case BREONSOIL:
    n_brdf_parameters = 1;

    do_factor_wfs = true;

    do_params_wfs.resize(n_brdf_parameters);
    do_params_wfs(0) = false;
    break;
  default:
    Exception e("Unhandled BRDF kernel type: ");
    e << which_brdf;
    throw e;
  }

  initialize_kernel_parameters(kernel_index, which_brdf, lambertian_flag, n_brdf_parameters, do_factor_wfs, do_params_wfs);

  // Count number of wfs and add to total(s)
  int n_factor_wfs = do_factor_wfs ? 1 : 0;
  int n_param_wfs = sum(where(do_params_wfs(Range::all()) == true, 1, 0));

  n_kernel_factor_wfs(n_kernel_factor_wfs() + n_factor_wfs);
  n_kernel_params_wfs(n_kernel_params_wfs() + n_param_wfs);
  n_surface_wfs(n_surface_wfs() + n_factor_wfs + n_param_wfs); 

  // Determine if any parameter wfs are computed for current kernel
  do_kparams_derivs(kernel_index, any(do_params_wfs(Range::all())));
}

//-----------------------------------------------------------------------
/// Sets up the BRDF inputs to be used by the BRDF calculation code
/// This routine is intended to be called for each spectral point
//-----------------------------------------------------------------------

ArrayAd<double, 1> SpurrBrdfDriver::setup_brdf_inputs(int surface_type, const ArrayAd<double, 1>& surface_parameters) const
{
  // Copy input surface parameters as the returned value may
  // be modified to properly account for changes to parameters
  // made for use inside of RT
  ArrayAd<double, 1> rt_surf_params(surface_parameters);

  // Map kernel parameter indexes to indexes in surface_parameter input array
  // Due to heritiage code these two will not always match up, but we need
  // to know which indexes of the paremeters array to properly account for jacobians
  Array<int, 1> parameter_indexes;

  switch (surface_type) {
  case LAMBERTIAN:
    parameter_indexes.resize(1);
    parameter_indexes(0) = 0;
    setup_lambertian_inputs(0, rt_surf_params, parameter_indexes);

    break;
  case COXMUNK:
    parameter_indexes.resize(4);
    parameter_indexes(0) = 0; // scale factor
    parameter_indexes(1) = 1; // windspeed
    parameter_indexes(2) = 2; // refractive index
    parameter_indexes(3) = 4; // shadowing
    setup_coxmunk_inputs(0, rt_surf_params, parameter_indexes);

    // lamberitan component to coxmunk
    parameter_indexes.resize(1);
    parameter_indexes(0) = 3; // albedo
    setup_lambertian_inputs(1, rt_surf_params, parameter_indexes);
    break;
  case BREONVEG:
  case BREONSOIL:
    parameter_indexes.resize(4);
    parameter_indexes(0) = 0; // rahman kernel factor
    parameter_indexes(1) = 1; // hotspot parameter
    parameter_indexes(2) = 2; // asymmetry
    parameter_indexes(3) = 3; // anisotropy_parameter
    setup_rahman_inputs(0, rt_surf_params, parameter_indexes);

    parameter_indexes.resize(1);
    parameter_indexes(0) = 4; // breon kernel factor
    setup_breon_inputs(1, rt_surf_params, parameter_indexes);
    break;
  default:
    Exception e("Unhandled BRDF type index: ");
    e << surface_type;
    throw e;
  }

  // Calculate BRDF base on all input setup 
  calculate_brdf();

  return rt_surf_params;
}

void SpurrBrdfDriver::setup_lambertian_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const Array<int, 1>& parameter_indexes) const
{
  // brdf_params value only used if do_brdf_surface = True
  // albedo set to brdf weighting function so that weighting function
  // calculated automatically. The brdf parameters for Lambertian
  // do not have a wf defined for them
  int albedo_idx = parameter_indexes(0);
  brdf_factors(kernel_index) = surface_parameters(albedo_idx).value();
  
  // According to LIDORT user's guide section "2.7.1 BRDFs as a sum of kernel functions" this should be 1.0
  brdf_params(kernel_index, 0) = 1.0;
}

void SpurrBrdfDriver::setup_coxmunk_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const Array<int, 1>& parameter_indexes) const
{
  // Modify surface_parameters in place so that jacobians reflect modifications
  int scale_idx  = parameter_indexes(0);
  int ws_idx     = parameter_indexes(1);
  int refr_idx   = parameter_indexes(2);
  int shadow_idx = parameter_indexes(3);

  brdf_factors(kernel_index) = surface_parameters(scale_idx).value();

  // windspeed
  surface_parameters(ws_idx) = 0.003 + 0.00512 * surface_parameters(ws_idx);
  brdf_params(kernel_index, 0) = surface_parameters(ws_idx).value();

  // refactive index of water (squared)
  surface_parameters(refr_idx) = surface_parameters(refr_idx) * surface_parameters(refr_idx); // refractive index of air
  brdf_params(kernel_index, 1) = surface_parameters(refr_idx).value();

  // shadowing flag, for kernel routine
  // kernel routine seems to use this value despite shadowing flag below
  brdf_params(kernel_index, 2) = surface_parameters(shadow_idx).value();

  // Shadowing flag that seems independent of value passed to kernel function,
  // or possibly value above is overwritten by this flag
  do_shadow_effect(surface_parameters(shadow_idx).value() > 0.0);
}

void SpurrBrdfDriver::setup_rahman_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const Array<int, 1>& parameter_indexes) const
{
  // Modify surface_parameters in place so that jacobians reflect modifications
  int kf_idx   = parameter_indexes(0);
  int ampl_idx = parameter_indexes(1);
  int asym_idx = parameter_indexes(2);
  int geom_idx = parameter_indexes(3);

  brdf_factors(kernel_index) = surface_parameters(kf_idx).value();

  // Overall amplitude
  brdf_params(kernel_index, 0) = surface_parameters(ampl_idx).value();
  // Asymmetry parameter
  brdf_params(kernel_index, 1) = surface_parameters(asym_idx).value();
  // Geometric factor
  brdf_params(kernel_index, 2) = surface_parameters(geom_idx).value();
}

void SpurrBrdfDriver::setup_breon_inputs(int kernel_index, ArrayAd<double, 1>& surface_parameters, const Array<int, 1>& parameter_indexes) const
{
  int kf_idx = parameter_indexes(0);
  
  brdf_factors(kernel_index) = surface_parameters(kf_idx).value();
  brdf_params(kernel_index, 0) = 2.25; // Refractive index squared, same as hardcoded value inside lrad, should be same as value in ground_brdf.cc
}

/********************************************************************/

//-----------------------------------------------------------------------
/// Calculates intensity value with the given inputs
//-----------------------------------------------------------------------

double SpurrRtDriver::reflectance_calculate(const Array<double, 1>& height_grid,
					 double sza, double azm, double zen,
					 int surface_type,
					 const Array<double, 1>& surface_parameters,
					 const Array<double, 1>& od, 
					 const Array<double, 1>& ssa,
					 const Array<double, 2>& pf)
{
  // Initialize scene 
  setup_height_grid(height_grid);
  brdf_driver_-> setup_geometry(sza, azm, zen);
  setup_geometry(sza, azm, zen);

  // Set up BRDF inputs, here we throw away the jacobian
  // value of the surface parameters
  ArrayAd<double, 1> surf_param_ad(surface_parameters.rows(), 0);
  surf_param_ad.value() = surface_parameters;
  ArrayAd<double, 1> lidort_surf = brdf_driver_->setup_brdf_inputs(surface_type, surf_param_ad);

  // Set up LIDORT inputs and run
  setup_optical_inputs(od, ssa, pf);
  clear_linear_inputs();
  calculate_rt();

  // Return answer
  return get_intensity();
}

//-----------------------------------------------------------------------
/// Calculates intensity, profile and surface weighting factors (jacobians)
/// with the given inputs
//-----------------------------------------------------------------------

void SpurrRtDriver::reflectance_and_jacobian_calculate(const Array<double, 1>& height_grid,
						    double sza, double azm, double zen,
						    int surface_type,
						    ArrayAd<double, 1>& surface_parameters,
						    const ArrayAd<double, 1>& od, 
						    const ArrayAd<double, 1>& ssa,
						    const ArrayAd<double, 2>& pf,
						    double& reflectance,
						    Array<double, 2>& jac_atm, 
						    Array<double, 1>& jac_surf)
{
  // Initialize scene 
  setup_height_grid(height_grid);
  brdf_driver_->setup_geometry(sza, azm, zen);
  setup_geometry(sza, azm, zen);

  // Set up BRDF inputs and run
  surface_parameters = brdf_driver_->setup_brdf_inputs(surface_type, surface_parameters);

  setup_optical_inputs(od.value(), ssa.value(), pf.value());
  bool do_surface_pd = true;
  setup_linear_inputs(od, ssa, pf, do_surface_pd);
  calculate_rt();

  // Copy values from LIDORT
  reflectance = get_intensity();
  copy_jacobians(jac_atm, jac_surf);
}
