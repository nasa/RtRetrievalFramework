#ifndef LIDORT_DRIVER_H
#define LIDORT_DRIVER_H

#include "spurr_driver.h"
#include "lidort_interface_masters.h"

namespace FullPhysics {

/****************************************************************//**
  LIDORT specific BRDF driver implementation
 *******************************************************************/

class LidortBrdfDriver : public SpurrBrdfDriver {
public:
  LidortBrdfDriver(int nstream, int nmoment);
  virtual ~LidortBrdfDriver() {}

  /// Interface to BRDF interface to allow changing configuration to values
  // other than default
  const boost::shared_ptr<Brdf_Lin_Sup_Masters> brdf_interface() const { return brdf_interface_; };

  virtual void setup_geometry(double sza, double azm, double zen) const;

  virtual int n_brdf_kernels() const;

  virtual int n_kernel_factor_wfs() const;
  virtual int n_kernel_params_wfs() const;
  virtual int n_surface_wfs() const;
  virtual bool do_kparams_derivs(const int kernel_index) const;
  virtual bool do_shadow_effect() const;

protected:
  int nmoment_;

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

  boost::shared_ptr<Brdf_Lin_Sup_Masters> brdf_interface_;
};

/****************************************************************//**
  LIDORT specific Radiative transfer interface implementation
 *******************************************************************/

class LidortRtDriver : public SpurrRtDriver {
public:
  LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type, const blitz::Array<double, 1>& zen, bool pure_nadir);

  int number_moment() const;
  int number_stream() const;

  void setup_sphericity(double zen) const;
  void set_plane_parallel() const;
  void set_pseudo_spherical() const;
  void set_plane_parallel_plus_ss_correction() const;
  void set_line_of_sight() const;

  bool do_multi_scatt_only() const { return do_multi_scatt_only_; }

  bool pure_nadir() const { return pure_nadir_; }

  /// Access to BRDF driver
  const boost::shared_ptr<LidortBrdfDriver> lidort_brdf_driver() const
  { return boost::shared_ptr<LidortBrdfDriver>(boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver_)); }

  const boost::shared_ptr<Brdf_Lin_Sup_Masters> brdf_interface() const
  { return boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver_)->brdf_interface(); }

  /// Interface to LIDORT RT software inputs to allow changing LIDORT configuration to values other than default
  const boost::shared_ptr<Lidort_Lps_Masters> lidort_interface() const { return lidort_interface_; }

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

protected:
  void initialize_rt();
  void copy_brdf_sup_outputs() const;

  int nstream_, nmoment_;
  bool do_multi_scatt_only_;
  int surface_type_;
  bool pure_nadir_;
  boost::shared_ptr<Lidort_Lps_Masters> lidort_interface_;
};

}
#endif
