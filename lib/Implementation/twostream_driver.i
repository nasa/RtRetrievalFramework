// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "twostream_driver.h"
%}
%import "twostream_interface.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::TwostreamBrdfDriver);
%fp_shared_ptr(FullPhysics::TwostreamRtDriver);

namespace FullPhysics {
class TwostreamBrdfDriver {
public:
  TwostreamBrdfDriver(int surface_type);
  virtual ~TwostreamBrdfDriver();

  virtual void setup_geometry(double sza, double azm, double zen) const;

  %python_attribute(n_brdf_kernels,int)
  %python_attribute(n_kernel_factor_wfs,int)
  %python_attribute(n_kernel_params_wfs,int)
  %python_attribute(n_surface_wfs,int)
  %python_attribute(do_shadow_effect,bool)
  %python_attribute(brdf_interface, boost::shared_ptr<Twostream_Ls_Brdf_Supplement>)
  virtual bool do_kparams_derivs(const int kernel_index) const;
};

class TwostreamRtDriver {
public:
  TwostreamRtDriver(int nlayers, int surface_type, bool do_fullquadrature = true);

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

  %python_attribute(twostream_brdf_driver, boost::shared_ptr<TwostreamBrdfDriver>)
  %python_attribute(brdf_interface, boost::shared_ptr<Twostream_Ls_Brdf_Supplement>)
  %python_attribute(twostream_interface, boost::shared_ptr<Twostream_Lps_Master>)

  %python_attribute(do_full_quadrature, bool)
};
}
