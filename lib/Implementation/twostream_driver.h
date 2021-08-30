#ifndef TWOSTREAM_DRIVER_H
#define TWOSTREAM_DRIVER_H

#include "spurr_driver.h"
#include "twostream_interface.h"

namespace FullPhysics {

/****************************************************************//**
  TwoStream specific BRDF driver implementation
 *******************************************************************/

class TwostreamBrdfDriver : public SpurrBrdfDriver {
public:
  TwostreamBrdfDriver(int surface_type);
  virtual ~TwostreamBrdfDriver() {}

  virtual void setup_geometry(double sza, double azm, double zen) const;

  virtual int n_brdf_kernels() const;

  virtual int n_kernel_factor_wfs() const;
  virtual int n_kernel_params_wfs() const;
  virtual int n_surface_wfs() const;
  virtual bool do_kparams_derivs(const int kernel_index) const;
  virtual bool do_shadow_effect() const;

  boost::shared_ptr<Twostream_Ls_Brdf_Supplement> brdf_interface() const { return twostream_brdf_; }

protected:
  virtual void calculate_brdf() const;

  virtual void n_brdf_kernels(const int n_kernels);
  virtual void n_kernel_factor_wfs(const int n_factors);
  virtual void n_kernel_params_wfs(const int n_params);
  virtual void n_surface_wfs(const int n_wfs);
  virtual void do_kparams_derivs(const int kernel_index, const bool do_kparams);
  virtual void do_shadow_effect(const bool do_shadow) const;

  virtual void initialize_kernel_parameters(const int kernel_index,
                                            const int which_brdf,
                                            const bool lambertian_flag,
                                            const int n_brdf_parameters,
                                            const bool do_factor_wfs,
                                            const blitz::Array<bool, 1>& do_params_wfs);

  boost::shared_ptr<Twostream_Ls_Brdf_Supplement> twostream_brdf_;
};

/****************************************************************//**
  TwoStream specific Radiative transfer driver implementation
 *******************************************************************/

class TwostreamRtDriver : public SpurrRtDriver {
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

  boost::shared_ptr<TwostreamBrdfDriver> twostream_brdf_driver() const { 
    return boost::shared_ptr<TwostreamBrdfDriver>(boost::dynamic_pointer_cast<TwostreamBrdfDriver>(brdf_driver_)); 
  }

  boost::shared_ptr<Twostream_Ls_Brdf_Supplement> brdf_interface() const { return twostream_brdf_driver()->brdf_interface(); }
 
  boost::shared_ptr<Twostream_Lps_Master> twostream_interface() const { return twostream_interface_; }
  
  bool do_full_quadrature() const { return do_fullquadrature_; };

protected:
  void initialize_rt();

  int surface_type_;
  bool do_fullquadrature_;

  boost::shared_ptr<Twostream_Lps_Master> twostream_interface_;
};

}
#endif
