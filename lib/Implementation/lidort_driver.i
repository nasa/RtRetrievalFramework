// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%include "spurr_driver.i"
%{
#include "lidort_driver.h"
%}
%import "lidort_interface_masters.i"
%import "lidort_interface_types.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::LidortBrdfDriver);
%fp_shared_ptr(FullPhysics::LidortRtDriver);

namespace FullPhysics {
class LidortBrdfDriver : public SpurrBrdfDriver {
public:
  LidortBrdfDriver(int nstream, int nmoment);
  virtual ~LidortBrdfDriver();

  %python_attribute(brdf_interface, boost::shared_ptr<Brdf_Lin_Sup_Masters>)

  virtual void setup_geometry(double sza, double azm, double zen) const;

  %python_attribute(n_brdf_kernels, int)
  %python_attribute(n_kernel_factor_wfs, int)
  %python_attribute(n_kernel_params_wfs, int)
  %python_attribute(n_surface_wfs, int)
  %python_attribute(do_shadow_effect, bool)
  virtual bool do_kparams_derivs(const int kernel_index) const;
};

class LidortRtDriver : public SpurrRtDriver {
public:
  LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type, const blitz::Array<double, 1>& zen, bool pure_nadir);

  %python_attribute(number_moment, int)
  %python_attribute(number_stream, int)
  void setup_sphericity(double zen) const;
  void set_plane_parallel() const;
  void set_pseudo_spherical() const;
  void set_plane_parallel_plus_ss_correction() const;
  void set_line_of_sight() const;

  %python_attribute(do_multi_scatt_only, bool)
  %python_attribute(pure_nadir, bool)

  /// Access to BRDF driver
  %python_attribute(lidort_brdf_driver, boost::shared_ptr<LidortBrdfDriver>)
  %python_attribute(brdf_interface, boost::shared_ptr<Brdf_Lin_Sup_Masters>)
  %python_attribute(lidort_interface, boost::shared_ptr<Lidort_Lps_Masters>)
  void setup_height_grid(const blitz::Array<double, 1>& height_grid) const;
  void setup_geometry(double sza, double azm, double zen) const;
  void setup_optical_inputs(const blitz::Array<double, 1>& od, 
			    const blitz::Array<double, 1>& ssa,
			    const blitz::Array<double, 2>& pf) const;
  void clear_linear_inputs() const;
  void setup_linear_inputs(const ArrayAd<double, 1>& od,
			   const ArrayAd<double, 1>& ssa,
			   const ArrayAd<double, 2>& pf,
			   bool do_surface_linearization) const;
  void calculate_rt() const;
  double get_intensity() const;
  void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf) const;
};

}
