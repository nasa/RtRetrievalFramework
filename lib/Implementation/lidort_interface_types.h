#ifndef LIDORT_INTERFACE_TYPES_H
#define LIDORT_INTERFACE_TYPES_H

#include <iostream>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "fp_exception.h"
#include <vector>

/* This file was auto-generated */

#define FORTRAN_TRUE_INT 1

#define BYTE_SIZE_ERROR_CHECK(var_name, c_size, f_size) \
  if(c_size != f_size) { \
    std::stringstream err_msg; \
    err_msg << "Size of C variable: " << c_size \
            << " for " << var_name \
            << " does not match size of Fortran variable: " << f_size; \
    throw Exception(err_msg.str()); \
  }

namespace FullPhysics {

extern "C" {
  void set_lidort_pars(struct Lidort_Pars *lidort_pars_struct_c);
}

/* This struct cannot inherit any properties from another struct/class which is a requirement of Fortran interoperability */

struct Lidort_Pars {

  const char lidort_version_number[3];
  const int lidort_inunit;
  const int lidort_scenunit;
  const int lidort_funit;
  const int lidort_resunit;
  const int lidort_errunit;
  const int lidort_dbgunit;
  const int max_messages;
  const int maxthreads;
  const int maxstreams;
  const int maxlayers;
  const int maxfinelayers;
  const int maxmoments_input;
  const int max_thermal_coeffs;
  const int maxbeams;
  const int max_user_streams;
  const int max_user_relazms;
  const int max_user_levels;
  const int max_partlayers;
  const int max_directions;
  const int max_brdf_kernels;
  const int max_brdf_parameters;
  const int maxstreams_brdf;
  const int max_msrs_muquad;
  const int max_msrs_phiquad;
  const int max_atmoswfs;
  const int max_surfacewfs;
  const int max_sleavewfs;
  const int max_geometries;
  const int max_allstrms;
  const int max_allstrms_p1;
  const int maxmoments;
  const int maxfourier;
  const int maxsthalf_brdf;
  const int maxstreams_2;
  const int maxstreams_p1;
  const int maxtotal;
  const int maxbandtotal;
  const double one;
  const double zero;
  const double onep5;
  const double two;
  const double three;
  const double four;
  const double quarter;
  const double half;
  const double minus_one;
  const double minus_two;
  const double pie;
  const double deg_to_rad;
  const double pi2;
  const double pi4;
  const double pio2;
  const double pio4;
  const double eps3;
  const double eps4;
  const double eps5;
  const double smallnum;
  const double bigexp;
  const double hopital_tolerance;
  const double omega_smallnum;
  const double max_tau_spath;
  const double max_tau_upath;
  const double max_tau_qpath;
  const int lidort_serious;
  const int lidort_warning;
  const int lidort_info;
  const int lidort_debug;
  const int lidort_success;
  const int upidx;
  const int dnidx;
  const int lambertian_idx;
  const int rossthin_idx;
  const int rossthick_idx;
  const int lisparse_idx;
  const int lidense_idx;
  const int hapke_idx;
  const int roujean_idx;
  const int rahman_idx;
  const int coxmunk_idx;
  const int breonveg_idx;
  const int breonsoil_idx;
  const int maxbrdf_idx;
  
  static Lidort_Pars& instance() {
    static Lidort_Pars obj;
    return obj;
  }

  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Pars &obj) {
    output_stream << "Lidort_Pars:" << std::endl
      << "lidort_version_number: " << "\"" << obj.lidort_version_number << "\"" << std::endl
      << "        lidort_inunit: " << obj.lidort_inunit  << std::endl
      << "      lidort_scenunit: " << obj.lidort_scenunit  << std::endl
      << "         lidort_funit: " << obj.lidort_funit  << std::endl
      << "       lidort_resunit: " << obj.lidort_resunit  << std::endl
      << "       lidort_errunit: " << obj.lidort_errunit  << std::endl
      << "       lidort_dbgunit: " << obj.lidort_dbgunit  << std::endl
      << "         max_messages: " << obj.max_messages  << std::endl
      << "           maxthreads: " << obj.maxthreads  << std::endl
      << "           maxstreams: " << obj.maxstreams  << std::endl
      << "            maxlayers: " << obj.maxlayers  << std::endl
      << "        maxfinelayers: " << obj.maxfinelayers  << std::endl
      << "     maxmoments_input: " << obj.maxmoments_input  << std::endl
      << "   max_thermal_coeffs: " << obj.max_thermal_coeffs  << std::endl
      << "             maxbeams: " << obj.maxbeams  << std::endl
      << "     max_user_streams: " << obj.max_user_streams  << std::endl
      << "     max_user_relazms: " << obj.max_user_relazms  << std::endl
      << "      max_user_levels: " << obj.max_user_levels  << std::endl
      << "       max_partlayers: " << obj.max_partlayers  << std::endl
      << "       max_directions: " << obj.max_directions  << std::endl
      << "     max_brdf_kernels: " << obj.max_brdf_kernels  << std::endl
      << "  max_brdf_parameters: " << obj.max_brdf_parameters  << std::endl
      << "      maxstreams_brdf: " << obj.maxstreams_brdf  << std::endl
      << "      max_msrs_muquad: " << obj.max_msrs_muquad  << std::endl
      << "     max_msrs_phiquad: " << obj.max_msrs_phiquad  << std::endl
      << "         max_atmoswfs: " << obj.max_atmoswfs  << std::endl
      << "       max_surfacewfs: " << obj.max_surfacewfs  << std::endl
      << "        max_sleavewfs: " << obj.max_sleavewfs  << std::endl
      << "       max_geometries: " << obj.max_geometries  << std::endl
      << "         max_allstrms: " << obj.max_allstrms  << std::endl
      << "      max_allstrms_p1: " << obj.max_allstrms_p1  << std::endl
      << "           maxmoments: " << obj.maxmoments  << std::endl
      << "           maxfourier: " << obj.maxfourier  << std::endl
      << "       maxsthalf_brdf: " << obj.maxsthalf_brdf  << std::endl
      << "         maxstreams_2: " << obj.maxstreams_2  << std::endl
      << "        maxstreams_p1: " << obj.maxstreams_p1  << std::endl
      << "             maxtotal: " << obj.maxtotal  << std::endl
      << "         maxbandtotal: " << obj.maxbandtotal  << std::endl
      << "                  one: " << obj.one  << std::endl
      << "                 zero: " << obj.zero  << std::endl
      << "                onep5: " << obj.onep5  << std::endl
      << "                  two: " << obj.two  << std::endl
      << "                three: " << obj.three  << std::endl
      << "                 four: " << obj.four  << std::endl
      << "              quarter: " << obj.quarter  << std::endl
      << "                 half: " << obj.half  << std::endl
      << "            minus_one: " << obj.minus_one  << std::endl
      << "            minus_two: " << obj.minus_two  << std::endl
      << "                  pie: " << obj.pie  << std::endl
      << "           deg_to_rad: " << obj.deg_to_rad  << std::endl
      << "                  pi2: " << obj.pi2  << std::endl
      << "                  pi4: " << obj.pi4  << std::endl
      << "                 pio2: " << obj.pio2  << std::endl
      << "                 pio4: " << obj.pio4  << std::endl
      << "                 eps3: " << obj.eps3  << std::endl
      << "                 eps4: " << obj.eps4  << std::endl
      << "                 eps5: " << obj.eps5  << std::endl
      << "             smallnum: " << obj.smallnum  << std::endl
      << "               bigexp: " << obj.bigexp  << std::endl
      << "    hopital_tolerance: " << obj.hopital_tolerance  << std::endl
      << "       omega_smallnum: " << obj.omega_smallnum  << std::endl
      << "        max_tau_spath: " << obj.max_tau_spath  << std::endl
      << "        max_tau_upath: " << obj.max_tau_upath  << std::endl
      << "        max_tau_qpath: " << obj.max_tau_qpath  << std::endl
      << "       lidort_serious: " << obj.lidort_serious  << std::endl
      << "       lidort_warning: " << obj.lidort_warning  << std::endl
      << "          lidort_info: " << obj.lidort_info  << std::endl
      << "         lidort_debug: " << obj.lidort_debug  << std::endl
      << "       lidort_success: " << obj.lidort_success  << std::endl
      << "                upidx: " << obj.upidx  << std::endl
      << "                dnidx: " << obj.dnidx  << std::endl
      << "       lambertian_idx: " << obj.lambertian_idx  << std::endl
      << "         rossthin_idx: " << obj.rossthin_idx  << std::endl
      << "        rossthick_idx: " << obj.rossthick_idx  << std::endl
      << "         lisparse_idx: " << obj.lisparse_idx  << std::endl
      << "          lidense_idx: " << obj.lidense_idx  << std::endl
      << "            hapke_idx: " << obj.hapke_idx  << std::endl
      << "          roujean_idx: " << obj.roujean_idx  << std::endl
      << "           rahman_idx: " << obj.rahman_idx  << std::endl
      << "          coxmunk_idx: " << obj.coxmunk_idx  << std::endl
      << "         breonveg_idx: " << obj.breonveg_idx  << std::endl
      << "        breonsoil_idx: " << obj.breonsoil_idx  << std::endl
      << "          maxbrdf_idx: " << obj.maxbrdf_idx  << std::endl;
    return output_stream;
  }

private:
  Lidort_Pars() : lidort_version_number(), lidort_inunit(0), lidort_scenunit(0), lidort_funit(0), lidort_resunit(0), lidort_errunit(0), lidort_dbgunit(0), max_messages(0), maxthreads(0), maxstreams(0), maxlayers(0), maxfinelayers(0), maxmoments_input(0), max_thermal_coeffs(0), maxbeams(0), max_user_streams(0), max_user_relazms(0), max_user_levels(0), max_partlayers(0), max_directions(0), max_brdf_kernels(0), max_brdf_parameters(0), maxstreams_brdf(0), max_msrs_muquad(0), max_msrs_phiquad(0), max_atmoswfs(0), max_surfacewfs(0), max_sleavewfs(0), max_geometries(0), max_allstrms(0), max_allstrms_p1(0), maxmoments(0), maxfourier(0), maxsthalf_brdf(0), maxstreams_2(0), maxstreams_p1(0), maxtotal(0), maxbandtotal(0), one(0), zero(0), onep5(0), two(0), three(0), four(0), quarter(0), half(0), minus_one(0), minus_two(0), pie(0), deg_to_rad(0), pi2(0), pi4(0), pio2(0), pio4(0), eps3(0), eps4(0), eps5(0), smallnum(0), bigexp(0), hopital_tolerance(0), omega_smallnum(0), max_tau_spath(0), max_tau_upath(0), max_tau_qpath(0), lidort_serious(0), lidort_warning(0), lidort_info(0), lidort_debug(0), lidort_success(0), upidx(0), dnidx(0), lambertian_idx(0), rossthin_idx(0), rossthick_idx(0), lisparse_idx(0), lidense_idx(0), hapke_idx(0), roujean_idx(0), rahman_idx(0), coxmunk_idx(0), breonveg_idx(0), breonsoil_idx(0), maxbrdf_idx(0) { 
    set_lidort_pars(this);
  }
};


class Lidort_Structure : public Printable<Lidort_Structure> {
public:
  Lidort_Structure() : fortran_type_c(0), owns_pointer(true) {}
  Lidort_Structure(void* allocated_f_type_c) : fortran_type_c(allocated_f_type_c), owns_pointer(false) {}
  void* fortran_type_ptr() { return fortran_type_c; }

  virtual void print(std::ostream &output_stream) const = 0;

protected:
  void *fortran_type_c;
  bool owns_pointer;
};

// Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
extern "C" {
  void brdf_linsup_inputs_c_alloc_init(struct brdf_linsup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_inputs_c_init_only(struct brdf_linsup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_linsup_inputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_linsup_inputs {
  int* bs_do_kernel_factor_wfs_;
  int bs_do_kernel_factor_wfs__f_shapes[1];
  int bs_do_kernel_factor_wfs__f_byte_size;

  int* bs_do_kernel_params_wfs_;
  int bs_do_kernel_params_wfs__f_shapes[2];
  int bs_do_kernel_params_wfs__f_byte_size;

  int* bs_do_kparams_derivs_;
  int bs_do_kparams_derivs__f_shapes[1];
  int bs_do_kparams_derivs__f_byte_size;

  int* bs_n_surface_wfs_;
  int bs_n_surface_wfs__f_byte_size;

  int* bs_n_kernel_factor_wfs_;
  int bs_n_kernel_factor_wfs__f_byte_size;

  int* bs_n_kernel_params_wfs_;
  int bs_n_kernel_params_wfs__f_byte_size;

  
};

// Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
class Brdf_Linsup_Inputs : public Lidort_Structure {
public:
  // Allocating constructor
  Brdf_Linsup_Inputs() : Lidort_Structure() {
    brdf_linsup_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Linsup_Inputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    brdf_linsup_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Linsup_Inputs() {
    if (owns_pointer)
      brdf_linsup_inputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<bool, 1> bs_do_kernel_factor_wfs() const {
    blitz::Array<bool,1> as_bool(bs_do_kernel_factor_wfs_.shape());
    as_bool = blitz::where(bs_do_kernel_factor_wfs_ != 0, true, false);
    return as_bool;
  }
  void bs_do_kernel_factor_wfs(const blitz::Array<bool, 1>& bs_do_kernel_factor_wfs_in) {
    blitz::Array<int,1> as_int(bs_do_kernel_factor_wfs_.shape());
    as_int = blitz::where(bs_do_kernel_factor_wfs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kernel_factor_wfs_ = as_int;
  }
  
  const blitz::Array<bool, 2> bs_do_kernel_params_wfs() const {
    blitz::Array<bool,2> as_bool(bs_do_kernel_params_wfs_.shape());
    as_bool = blitz::where(bs_do_kernel_params_wfs_ != 0, true, false);
    return as_bool;
  }
  void bs_do_kernel_params_wfs(const blitz::Array<bool, 2>& bs_do_kernel_params_wfs_in) {
    blitz::Array<int,2> as_int(bs_do_kernel_params_wfs_.shape());
    as_int = blitz::where(bs_do_kernel_params_wfs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kernel_params_wfs_ = as_int;
  }
  
  const blitz::Array<bool, 1> bs_do_kparams_derivs() const {
    blitz::Array<bool,1> as_bool(bs_do_kparams_derivs_.shape());
    as_bool = blitz::where(bs_do_kparams_derivs_ != 0, true, false);
    return as_bool;
  }
  void bs_do_kparams_derivs(const blitz::Array<bool, 1>& bs_do_kparams_derivs_in) {
    blitz::Array<int,1> as_int(bs_do_kparams_derivs_.shape());
    as_int = blitz::where(bs_do_kparams_derivs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kparams_derivs_ = as_int;
  }
  
  const int& bs_n_surface_wfs() const {
    return *transfer_struct_c.bs_n_surface_wfs_;
  }
  void bs_n_surface_wfs(const int& bs_n_surface_wfs_in) {
    *transfer_struct_c.bs_n_surface_wfs_ = bs_n_surface_wfs_in;
  }
  
  const int& bs_n_kernel_factor_wfs() const {
    return *transfer_struct_c.bs_n_kernel_factor_wfs_;
  }
  void bs_n_kernel_factor_wfs(const int& bs_n_kernel_factor_wfs_in) {
    *transfer_struct_c.bs_n_kernel_factor_wfs_ = bs_n_kernel_factor_wfs_in;
  }
  
  const int& bs_n_kernel_params_wfs() const {
    return *transfer_struct_c.bs_n_kernel_params_wfs_;
  }
  void bs_n_kernel_params_wfs(const int& bs_n_kernel_params_wfs_in) {
    *transfer_struct_c.bs_n_kernel_params_wfs_ = bs_n_kernel_params_wfs_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Linsup_Inputs:" << std::endl
      << "bs_do_kernel_factor_wfs: " << std::endl << bs_do_kernel_factor_wfs()  << std::endl
      << "bs_do_kernel_params_wfs: " << std::endl << bs_do_kernel_params_wfs()  << std::endl
      << "   bs_do_kparams_derivs: " << std::endl << bs_do_kparams_derivs()  << std::endl
      << "       bs_n_surface_wfs: " << bs_n_surface_wfs()  << std::endl
      << " bs_n_kernel_factor_wfs: " << bs_n_kernel_factor_wfs()  << std::endl
      << " bs_n_kernel_params_wfs: " << bs_n_kernel_params_wfs()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_do_kernel_factor_wfs_",sizeof(*transfer_struct_c.bs_do_kernel_factor_wfs_),transfer_struct_c.bs_do_kernel_factor_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_kernel_params_wfs_",sizeof(*transfer_struct_c.bs_do_kernel_params_wfs_),transfer_struct_c.bs_do_kernel_params_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_kparams_derivs_",sizeof(*transfer_struct_c.bs_do_kparams_derivs_),transfer_struct_c.bs_do_kparams_derivs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_surface_wfs_",sizeof(*transfer_struct_c.bs_n_surface_wfs_),transfer_struct_c.bs_n_surface_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_kernel_factor_wfs_",sizeof(*transfer_struct_c.bs_n_kernel_factor_wfs_),transfer_struct_c.bs_n_kernel_factor_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_kernel_params_wfs_",sizeof(*transfer_struct_c.bs_n_kernel_params_wfs_),transfer_struct_c.bs_n_kernel_params_wfs__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_do_kernel_factor_wfs_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_do_kernel_factor_wfs_,
      blitz::shape(transfer_struct_c.bs_do_kernel_factor_wfs__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_do_kernel_params_wfs_.reference(blitz::Array<int, 2>(transfer_struct_c.bs_do_kernel_params_wfs_,
      blitz::shape(transfer_struct_c.bs_do_kernel_params_wfs__f_shapes[0],
                   transfer_struct_c.bs_do_kernel_params_wfs__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_do_kparams_derivs_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_do_kparams_derivs_,
      blitz::shape(transfer_struct_c.bs_do_kparams_derivs__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_linsup_inputs transfer_struct_c;

  blitz::Array<int, 1> bs_do_kernel_factor_wfs_;
  blitz::Array<int, 2> bs_do_kernel_params_wfs_;
  blitz::Array<int, 1> bs_do_kparams_derivs_;
  
};

// Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
extern "C" {
  void brdf_linsup_outputs_c_alloc_init(struct brdf_linsup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_outputs_c_init_only(struct brdf_linsup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_linsup_outputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_linsup_outputs {
  double* bs_ls_exactdb_brdfunc_;
  int bs_ls_exactdb_brdfunc__f_shapes[4];
  int bs_ls_exactdb_brdfunc__f_byte_size;

  double* bs_ls_brdf_f_0_;
  int bs_ls_brdf_f_0__f_shapes[4];
  int bs_ls_brdf_f_0__f_byte_size;

  double* bs_ls_brdf_f_;
  int bs_ls_brdf_f__f_shapes[4];
  int bs_ls_brdf_f__f_byte_size;

  double* bs_ls_user_brdf_f_0_;
  int bs_ls_user_brdf_f_0__f_shapes[4];
  int bs_ls_user_brdf_f_0__f_byte_size;

  double* bs_ls_user_brdf_f_;
  int bs_ls_user_brdf_f__f_shapes[4];
  int bs_ls_user_brdf_f__f_byte_size;

  double* bs_ls_emissivity_;
  int bs_ls_emissivity__f_shapes[3];
  int bs_ls_emissivity__f_byte_size;

  double* bs_ls_user_emissivity_;
  int bs_ls_user_emissivity__f_shapes[3];
  int bs_ls_user_emissivity__f_byte_size;

  
};

// Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
class Brdf_Linsup_Outputs : public Lidort_Structure {
public:
  // Allocating constructor
  Brdf_Linsup_Outputs() : Lidort_Structure() {
    brdf_linsup_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Linsup_Outputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    brdf_linsup_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Linsup_Outputs() {
    if (owns_pointer)
      brdf_linsup_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& bs_ls_exactdb_brdfunc() const {
    return bs_ls_exactdb_brdfunc_;
  }
  void bs_ls_exactdb_brdfunc(const blitz::Array<double, 4>& bs_ls_exactdb_brdfunc_in) {
    bs_ls_exactdb_brdfunc_ = bs_ls_exactdb_brdfunc_in;
  }
  
  const blitz::Array<double, 4>& bs_ls_brdf_f_0() const {
    return bs_ls_brdf_f_0_;
  }
  void bs_ls_brdf_f_0(const blitz::Array<double, 4>& bs_ls_brdf_f_0_in) {
    bs_ls_brdf_f_0_ = bs_ls_brdf_f_0_in;
  }
  
  const blitz::Array<double, 4>& bs_ls_brdf_f() const {
    return bs_ls_brdf_f_;
  }
  void bs_ls_brdf_f(const blitz::Array<double, 4>& bs_ls_brdf_f_in) {
    bs_ls_brdf_f_ = bs_ls_brdf_f_in;
  }
  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f_0() const {
    return bs_ls_user_brdf_f_0_;
  }
  void bs_ls_user_brdf_f_0(const blitz::Array<double, 4>& bs_ls_user_brdf_f_0_in) {
    bs_ls_user_brdf_f_0_ = bs_ls_user_brdf_f_0_in;
  }
  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f() const {
    return bs_ls_user_brdf_f_;
  }
  void bs_ls_user_brdf_f(const blitz::Array<double, 4>& bs_ls_user_brdf_f_in) {
    bs_ls_user_brdf_f_ = bs_ls_user_brdf_f_in;
  }
  
  const blitz::Array<double, 3>& bs_ls_emissivity() const {
    return bs_ls_emissivity_;
  }
  void bs_ls_emissivity(const blitz::Array<double, 3>& bs_ls_emissivity_in) {
    bs_ls_emissivity_ = bs_ls_emissivity_in;
  }
  
  const blitz::Array<double, 3>& bs_ls_user_emissivity() const {
    return bs_ls_user_emissivity_;
  }
  void bs_ls_user_emissivity(const blitz::Array<double, 3>& bs_ls_user_emissivity_in) {
    bs_ls_user_emissivity_ = bs_ls_user_emissivity_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Linsup_Outputs:" << std::endl
      << "bs_ls_exactdb_brdfunc: " << std::endl << bs_ls_exactdb_brdfunc()  << std::endl
      << "       bs_ls_brdf_f_0: " << std::endl << bs_ls_brdf_f_0()  << std::endl
      << "         bs_ls_brdf_f: " << std::endl << bs_ls_brdf_f()  << std::endl
      << "  bs_ls_user_brdf_f_0: " << std::endl << bs_ls_user_brdf_f_0()  << std::endl
      << "    bs_ls_user_brdf_f: " << std::endl << bs_ls_user_brdf_f()  << std::endl
      << "     bs_ls_emissivity: " << std::endl << bs_ls_emissivity()  << std::endl
      << "bs_ls_user_emissivity: " << std::endl << bs_ls_user_emissivity()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_ls_exactdb_brdfunc_",sizeof(*transfer_struct_c.bs_ls_exactdb_brdfunc_),transfer_struct_c.bs_ls_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_brdf_f_0_",sizeof(*transfer_struct_c.bs_ls_brdf_f_0_),transfer_struct_c.bs_ls_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_brdf_f_",sizeof(*transfer_struct_c.bs_ls_brdf_f_),transfer_struct_c.bs_ls_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_brdf_f_0_",sizeof(*transfer_struct_c.bs_ls_user_brdf_f_0_),transfer_struct_c.bs_ls_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_brdf_f_",sizeof(*transfer_struct_c.bs_ls_user_brdf_f_),transfer_struct_c.bs_ls_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_emissivity_",sizeof(*transfer_struct_c.bs_ls_emissivity_),transfer_struct_c.bs_ls_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_emissivity_",sizeof(*transfer_struct_c.bs_ls_user_emissivity_),transfer_struct_c.bs_ls_user_emissivity__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_ls_exactdb_brdfunc_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.bs_ls_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.bs_ls_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.bs_ls_exactdb_brdfunc__f_shapes[2],
                   transfer_struct_c.bs_ls_exactdb_brdfunc__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_ls_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[2],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_brdf_f_,
      blitz::shape(transfer_struct_c.bs_ls_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[2],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_user_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[2],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_user_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_user_brdf_f_,
      blitz::shape(transfer_struct_c.bs_ls_user_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[2],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_emissivity_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_ls_emissivity_,
      blitz::shape(transfer_struct_c.bs_ls_emissivity__f_shapes[0],
                   transfer_struct_c.bs_ls_emissivity__f_shapes[1],
                   transfer_struct_c.bs_ls_emissivity__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_ls_user_emissivity_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_ls_user_emissivity_,
      blitz::shape(transfer_struct_c.bs_ls_user_emissivity__f_shapes[0],
                   transfer_struct_c.bs_ls_user_emissivity__f_shapes[1],
                   transfer_struct_c.bs_ls_user_emissivity__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_linsup_outputs transfer_struct_c;

  blitz::Array<double, 4> bs_ls_exactdb_brdfunc_;
  blitz::Array<double, 4> bs_ls_brdf_f_0_;
  blitz::Array<double, 4> bs_ls_brdf_f_;
  blitz::Array<double, 4> bs_ls_user_brdf_f_0_;
  blitz::Array<double, 4> bs_ls_user_brdf_f_;
  blitz::Array<double, 3> bs_ls_emissivity_;
  blitz::Array<double, 3> bs_ls_user_emissivity_;
  
};

// Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
extern "C" {
  void brdf_sup_inputs_c_alloc_init(struct brdf_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_inputs_c_init_only(struct brdf_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_sup_inputs_c_destroy(void **fortran_type_c);
  void brdf_sup_inputs_bs_brdf_names_get(void **fortran_type_c, const int* bs_brdf_names_in_shape_1, const int* bs_brdf_names_in_len, const char* bs_brdf_names_in);
  
}

struct brdf_sup_inputs {
  int* bs_do_user_streams_;
  int bs_do_user_streams__f_byte_size;

  int* bs_do_brdf_surface_;
  int bs_do_brdf_surface__f_byte_size;

  int* bs_do_surface_emission_;
  int bs_do_surface_emission__f_byte_size;

  int* bs_nstreams_;
  int bs_nstreams__f_byte_size;

  int* bs_nbeams_;
  int bs_nbeams__f_byte_size;

  double* bs_beam_szas_;
  int bs_beam_szas__f_shapes[1];
  int bs_beam_szas__f_byte_size;

  int* bs_n_user_relazms_;
  int bs_n_user_relazms__f_byte_size;

  double* bs_user_relazms_;
  int bs_user_relazms__f_shapes[1];
  int bs_user_relazms__f_byte_size;

  int* bs_n_user_streams_;
  int bs_n_user_streams__f_byte_size;

  double* bs_user_angles_input_;
  int bs_user_angles_input__f_shapes[1];
  int bs_user_angles_input__f_byte_size;

  int* bs_n_brdf_kernels_;
  int bs_n_brdf_kernels__f_byte_size;

  
  int bs_brdf_names__f_shapes[1];
  int bs_brdf_names__f_len;

  int* bs_which_brdf_;
  int bs_which_brdf__f_shapes[1];
  int bs_which_brdf__f_byte_size;

  int* bs_n_brdf_parameters_;
  int bs_n_brdf_parameters__f_shapes[1];
  int bs_n_brdf_parameters__f_byte_size;

  double* bs_brdf_parameters_;
  int bs_brdf_parameters__f_shapes[2];
  int bs_brdf_parameters__f_byte_size;

  int* bs_lambertian_kernel_flag_;
  int bs_lambertian_kernel_flag__f_shapes[1];
  int bs_lambertian_kernel_flag__f_byte_size;

  double* bs_brdf_factors_;
  int bs_brdf_factors__f_shapes[1];
  int bs_brdf_factors__f_byte_size;

  int* bs_nstreams_brdf_;
  int bs_nstreams_brdf__f_byte_size;

  int* bs_do_shadow_effect_;
  int bs_do_shadow_effect__f_byte_size;

  int* bs_do_exactonly_;
  int bs_do_exactonly__f_byte_size;

  int* bs_do_glitter_msrcorr_;
  int bs_do_glitter_msrcorr__f_byte_size;

  int* bs_do_glitter_msrcorr_exactonly_;
  int bs_do_glitter_msrcorr_exactonly__f_byte_size;

  int* bs_glitter_msrcorr_order_;
  int bs_glitter_msrcorr_order__f_byte_size;

  int* bs_glitter_msrcorr_nmuquad_;
  int bs_glitter_msrcorr_nmuquad__f_byte_size;

  int* bs_glitter_msrcorr_nphiquad_;
  int bs_glitter_msrcorr_nphiquad__f_byte_size;

  
};

// Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
class Brdf_Sup_Inputs : public Lidort_Structure {
public:
  // Allocating constructor
  Brdf_Sup_Inputs() : Lidort_Structure() {
    brdf_sup_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Sup_Inputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    brdf_sup_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Sup_Inputs() {
    if (owns_pointer)
      brdf_sup_inputs_c_destroy(&fortran_type_c);
  }

  const bool bs_do_user_streams() const {
    return *transfer_struct_c.bs_do_user_streams_ != 0;
  }
  void bs_do_user_streams(const bool& bs_do_user_streams_in) {
    *transfer_struct_c.bs_do_user_streams_ = bs_do_user_streams_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool bs_do_brdf_surface() const {
    return *transfer_struct_c.bs_do_brdf_surface_ != 0;
  }
  void bs_do_brdf_surface(const bool& bs_do_brdf_surface_in) {
    *transfer_struct_c.bs_do_brdf_surface_ = bs_do_brdf_surface_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool bs_do_surface_emission() const {
    return *transfer_struct_c.bs_do_surface_emission_ != 0;
  }
  void bs_do_surface_emission(const bool& bs_do_surface_emission_in) {
    *transfer_struct_c.bs_do_surface_emission_ = bs_do_surface_emission_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const int& bs_nstreams() const {
    return *transfer_struct_c.bs_nstreams_;
  }
  void bs_nstreams(const int& bs_nstreams_in) {
    *transfer_struct_c.bs_nstreams_ = bs_nstreams_in;
  }
  
  const int& bs_nbeams() const {
    return *transfer_struct_c.bs_nbeams_;
  }
  void bs_nbeams(const int& bs_nbeams_in) {
    *transfer_struct_c.bs_nbeams_ = bs_nbeams_in;
  }
  
  const blitz::Array<double, 1>& bs_beam_szas() const {
    return bs_beam_szas_;
  }
  void bs_beam_szas(const blitz::Array<double, 1>& bs_beam_szas_in) {
    bs_beam_szas_ = bs_beam_szas_in;
  }
  
  const int& bs_n_user_relazms() const {
    return *transfer_struct_c.bs_n_user_relazms_;
  }
  void bs_n_user_relazms(const int& bs_n_user_relazms_in) {
    *transfer_struct_c.bs_n_user_relazms_ = bs_n_user_relazms_in;
  }
  
  const blitz::Array<double, 1>& bs_user_relazms() const {
    return bs_user_relazms_;
  }
  void bs_user_relazms(const blitz::Array<double, 1>& bs_user_relazms_in) {
    bs_user_relazms_ = bs_user_relazms_in;
  }
  
  const int& bs_n_user_streams() const {
    return *transfer_struct_c.bs_n_user_streams_;
  }
  void bs_n_user_streams(const int& bs_n_user_streams_in) {
    *transfer_struct_c.bs_n_user_streams_ = bs_n_user_streams_in;
  }
  
  const blitz::Array<double, 1>& bs_user_angles_input() const {
    return bs_user_angles_input_;
  }
  void bs_user_angles_input(const blitz::Array<double, 1>& bs_user_angles_input_in) {
    bs_user_angles_input_ = bs_user_angles_input_in;
  }
  
  const int& bs_n_brdf_kernels() const {
    return *transfer_struct_c.bs_n_brdf_kernels_;
  }
  void bs_n_brdf_kernels(const int& bs_n_brdf_kernels_in) {
    *transfer_struct_c.bs_n_brdf_kernels_ = bs_n_brdf_kernels_in;
  }
  
  const std::vector< std::string > bs_brdf_names() const {
    std::vector< std::string > bs_brdf_names_ret;
    blitz::Array<char, 2> bs_brdf_names_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_brdf_names__f_shapes[0], transfer_struct_c.bs_brdf_names__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_sup_inputs_bs_brdf_names_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_brdf_names__f_shapes[0], &transfer_struct_c.bs_brdf_names__f_len, bs_brdf_names_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_brdf_names_lcl.extent(0); dim_0_idx++)
      bs_brdf_names_ret.push_back( std::string(std::string(bs_brdf_names_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_brdf_names_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_brdf_names_ret;
  }
  
  const blitz::Array<int, 1>& bs_which_brdf() const {
    return bs_which_brdf_;
  }
  void bs_which_brdf(const blitz::Array<int, 1>& bs_which_brdf_in) {
    bs_which_brdf_ = bs_which_brdf_in;
  }
  
  const blitz::Array<int, 1>& bs_n_brdf_parameters() const {
    return bs_n_brdf_parameters_;
  }
  void bs_n_brdf_parameters(const blitz::Array<int, 1>& bs_n_brdf_parameters_in) {
    bs_n_brdf_parameters_ = bs_n_brdf_parameters_in;
  }
  
  const blitz::Array<double, 2>& bs_brdf_parameters() const {
    return bs_brdf_parameters_;
  }
  void bs_brdf_parameters(const blitz::Array<double, 2>& bs_brdf_parameters_in) {
    bs_brdf_parameters_ = bs_brdf_parameters_in;
  }
  
  const blitz::Array<bool, 1> bs_lambertian_kernel_flag() const {
    blitz::Array<bool,1> as_bool(bs_lambertian_kernel_flag_.shape());
    as_bool = blitz::where(bs_lambertian_kernel_flag_ != 0, true, false);
    return as_bool;
  }
  void bs_lambertian_kernel_flag(const blitz::Array<bool, 1>& bs_lambertian_kernel_flag_in) {
    blitz::Array<int,1> as_int(bs_lambertian_kernel_flag_.shape());
    as_int = blitz::where(bs_lambertian_kernel_flag_in == true, FORTRAN_TRUE_INT, 0);
    bs_lambertian_kernel_flag_ = as_int;
  }
  
  const blitz::Array<double, 1>& bs_brdf_factors() const {
    return bs_brdf_factors_;
  }
  void bs_brdf_factors(const blitz::Array<double, 1>& bs_brdf_factors_in) {
    bs_brdf_factors_ = bs_brdf_factors_in;
  }
  
  const int& bs_nstreams_brdf() const {
    return *transfer_struct_c.bs_nstreams_brdf_;
  }
  void bs_nstreams_brdf(const int& bs_nstreams_brdf_in) {
    *transfer_struct_c.bs_nstreams_brdf_ = bs_nstreams_brdf_in;
  }
  
  const bool bs_do_shadow_effect() const {
    return *transfer_struct_c.bs_do_shadow_effect_ != 0;
  }
  void bs_do_shadow_effect(const bool& bs_do_shadow_effect_in) {
    *transfer_struct_c.bs_do_shadow_effect_ = bs_do_shadow_effect_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool bs_do_exactonly() const {
    return *transfer_struct_c.bs_do_exactonly_ != 0;
  }
  void bs_do_exactonly(const bool& bs_do_exactonly_in) {
    *transfer_struct_c.bs_do_exactonly_ = bs_do_exactonly_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool bs_do_glitter_msrcorr() const {
    return *transfer_struct_c.bs_do_glitter_msrcorr_ != 0;
  }
  void bs_do_glitter_msrcorr(const bool& bs_do_glitter_msrcorr_in) {
    *transfer_struct_c.bs_do_glitter_msrcorr_ = bs_do_glitter_msrcorr_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool bs_do_glitter_msrcorr_exactonly() const {
    return *transfer_struct_c.bs_do_glitter_msrcorr_exactonly_ != 0;
  }
  void bs_do_glitter_msrcorr_exactonly(const bool& bs_do_glitter_msrcorr_exactonly_in) {
    *transfer_struct_c.bs_do_glitter_msrcorr_exactonly_ = bs_do_glitter_msrcorr_exactonly_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const int& bs_glitter_msrcorr_order() const {
    return *transfer_struct_c.bs_glitter_msrcorr_order_;
  }
  void bs_glitter_msrcorr_order(const int& bs_glitter_msrcorr_order_in) {
    *transfer_struct_c.bs_glitter_msrcorr_order_ = bs_glitter_msrcorr_order_in;
  }
  
  const int& bs_glitter_msrcorr_nmuquad() const {
    return *transfer_struct_c.bs_glitter_msrcorr_nmuquad_;
  }
  void bs_glitter_msrcorr_nmuquad(const int& bs_glitter_msrcorr_nmuquad_in) {
    *transfer_struct_c.bs_glitter_msrcorr_nmuquad_ = bs_glitter_msrcorr_nmuquad_in;
  }
  
  const int& bs_glitter_msrcorr_nphiquad() const {
    return *transfer_struct_c.bs_glitter_msrcorr_nphiquad_;
  }
  void bs_glitter_msrcorr_nphiquad(const int& bs_glitter_msrcorr_nphiquad_in) {
    *transfer_struct_c.bs_glitter_msrcorr_nphiquad_ = bs_glitter_msrcorr_nphiquad_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Sup_Inputs:" << std::endl
      << "             bs_do_user_streams: " << bs_do_user_streams()  << std::endl
      << "             bs_do_brdf_surface: " << bs_do_brdf_surface()  << std::endl
      << "         bs_do_surface_emission: " << bs_do_surface_emission()  << std::endl
      << "                    bs_nstreams: " << bs_nstreams()  << std::endl
      << "                      bs_nbeams: " << bs_nbeams()  << std::endl
      << "                   bs_beam_szas: " << std::endl << bs_beam_szas()  << std::endl
      << "              bs_n_user_relazms: " << bs_n_user_relazms()  << std::endl
      << "                bs_user_relazms: " << std::endl << bs_user_relazms()  << std::endl
      << "              bs_n_user_streams: " << bs_n_user_streams()  << std::endl
      << "           bs_user_angles_input: " << std::endl << bs_user_angles_input()  << std::endl
      << "              bs_n_brdf_kernels: " << bs_n_brdf_kernels()  << std::endl
      << "                  bs_brdf_names: " << std::endl;
    std::vector< std::string > bs_brdf_names_lcl = bs_brdf_names();
    for(unsigned int idx = 0; idx < bs_brdf_names_lcl.size(); idx++)
      if ( bs_brdf_names_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_brdf_names_lcl[idx] << "\"" << std::endl;
    output_stream
      << "                  bs_which_brdf: " << std::endl << bs_which_brdf()  << std::endl
      << "           bs_n_brdf_parameters: " << std::endl << bs_n_brdf_parameters()  << std::endl
      << "             bs_brdf_parameters: " << std::endl << bs_brdf_parameters()  << std::endl
      << "      bs_lambertian_kernel_flag: " << std::endl << bs_lambertian_kernel_flag()  << std::endl
      << "                bs_brdf_factors: " << std::endl << bs_brdf_factors()  << std::endl
      << "               bs_nstreams_brdf: " << bs_nstreams_brdf()  << std::endl
      << "            bs_do_shadow_effect: " << bs_do_shadow_effect()  << std::endl
      << "                bs_do_exactonly: " << bs_do_exactonly()  << std::endl
      << "          bs_do_glitter_msrcorr: " << bs_do_glitter_msrcorr()  << std::endl
      << "bs_do_glitter_msrcorr_exactonly: " << bs_do_glitter_msrcorr_exactonly()  << std::endl
      << "       bs_glitter_msrcorr_order: " << bs_glitter_msrcorr_order()  << std::endl
      << "     bs_glitter_msrcorr_nmuquad: " << bs_glitter_msrcorr_nmuquad()  << std::endl
      << "    bs_glitter_msrcorr_nphiquad: " << bs_glitter_msrcorr_nphiquad()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_do_user_streams_",sizeof(*transfer_struct_c.bs_do_user_streams_),transfer_struct_c.bs_do_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_brdf_surface_",sizeof(*transfer_struct_c.bs_do_brdf_surface_),transfer_struct_c.bs_do_brdf_surface__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_surface_emission_",sizeof(*transfer_struct_c.bs_do_surface_emission_),transfer_struct_c.bs_do_surface_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nstreams_",sizeof(*transfer_struct_c.bs_nstreams_),transfer_struct_c.bs_nstreams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nbeams_",sizeof(*transfer_struct_c.bs_nbeams_),transfer_struct_c.bs_nbeams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_beam_szas_",sizeof(*transfer_struct_c.bs_beam_szas_),transfer_struct_c.bs_beam_szas__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_relazms_",sizeof(*transfer_struct_c.bs_n_user_relazms_),transfer_struct_c.bs_n_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_relazms_",sizeof(*transfer_struct_c.bs_user_relazms_),transfer_struct_c.bs_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_streams_",sizeof(*transfer_struct_c.bs_n_user_streams_),transfer_struct_c.bs_n_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_angles_input_",sizeof(*transfer_struct_c.bs_user_angles_input_),transfer_struct_c.bs_user_angles_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_brdf_kernels_",sizeof(*transfer_struct_c.bs_n_brdf_kernels_),transfer_struct_c.bs_n_brdf_kernels__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_which_brdf_",sizeof(*transfer_struct_c.bs_which_brdf_),transfer_struct_c.bs_which_brdf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_brdf_parameters_",sizeof(*transfer_struct_c.bs_n_brdf_parameters_),transfer_struct_c.bs_n_brdf_parameters__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_parameters_",sizeof(*transfer_struct_c.bs_brdf_parameters_),transfer_struct_c.bs_brdf_parameters__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_lambertian_kernel_flag_",sizeof(*transfer_struct_c.bs_lambertian_kernel_flag_),transfer_struct_c.bs_lambertian_kernel_flag__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_factors_",sizeof(*transfer_struct_c.bs_brdf_factors_),transfer_struct_c.bs_brdf_factors__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nstreams_brdf_",sizeof(*transfer_struct_c.bs_nstreams_brdf_),transfer_struct_c.bs_nstreams_brdf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_shadow_effect_",sizeof(*transfer_struct_c.bs_do_shadow_effect_),transfer_struct_c.bs_do_shadow_effect__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_exactonly_",sizeof(*transfer_struct_c.bs_do_exactonly_),transfer_struct_c.bs_do_exactonly__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_glitter_msrcorr_",sizeof(*transfer_struct_c.bs_do_glitter_msrcorr_),transfer_struct_c.bs_do_glitter_msrcorr__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_glitter_msrcorr_exactonly_",sizeof(*transfer_struct_c.bs_do_glitter_msrcorr_exactonly_),transfer_struct_c.bs_do_glitter_msrcorr_exactonly__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_order_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_order_),transfer_struct_c.bs_glitter_msrcorr_order__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_nmuquad_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_nmuquad_),transfer_struct_c.bs_glitter_msrcorr_nmuquad__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_nphiquad_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_nphiquad_),transfer_struct_c.bs_glitter_msrcorr_nphiquad__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_beam_szas_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_beam_szas_,
      blitz::shape(transfer_struct_c.bs_beam_szas__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_user_relazms_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_user_relazms_,
      blitz::shape(transfer_struct_c.bs_user_relazms__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_user_angles_input_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_user_angles_input_,
      blitz::shape(transfer_struct_c.bs_user_angles_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_which_brdf_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_which_brdf_,
      blitz::shape(transfer_struct_c.bs_which_brdf__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_n_brdf_parameters_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_n_brdf_parameters_,
      blitz::shape(transfer_struct_c.bs_n_brdf_parameters__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_brdf_parameters_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_brdf_parameters_,
      blitz::shape(transfer_struct_c.bs_brdf_parameters__f_shapes[0],
                   transfer_struct_c.bs_brdf_parameters__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_lambertian_kernel_flag_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_lambertian_kernel_flag_,
      blitz::shape(transfer_struct_c.bs_lambertian_kernel_flag__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_brdf_factors_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_brdf_factors_,
      blitz::shape(transfer_struct_c.bs_brdf_factors__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_sup_inputs transfer_struct_c;

  blitz::Array<double, 1> bs_beam_szas_;
  blitz::Array<double, 1> bs_user_relazms_;
  blitz::Array<double, 1> bs_user_angles_input_;
  blitz::Array<int, 1> bs_which_brdf_;
  blitz::Array<int, 1> bs_n_brdf_parameters_;
  blitz::Array<double, 2> bs_brdf_parameters_;
  blitz::Array<int, 1> bs_lambertian_kernel_flag_;
  blitz::Array<double, 1> bs_brdf_factors_;
  
};

// Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
extern "C" {
  void brdf_sup_outputs_c_alloc_init(struct brdf_sup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_outputs_c_init_only(struct brdf_sup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_sup_outputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_sup_outputs {
  double* bs_exactdb_brdfunc_;
  int bs_exactdb_brdfunc__f_shapes[3];
  int bs_exactdb_brdfunc__f_byte_size;

  double* bs_brdf_f_0_;
  int bs_brdf_f_0__f_shapes[3];
  int bs_brdf_f_0__f_byte_size;

  double* bs_brdf_f_;
  int bs_brdf_f__f_shapes[3];
  int bs_brdf_f__f_byte_size;

  double* bs_user_brdf_f_0_;
  int bs_user_brdf_f_0__f_shapes[3];
  int bs_user_brdf_f_0__f_byte_size;

  double* bs_user_brdf_f_;
  int bs_user_brdf_f__f_shapes[3];
  int bs_user_brdf_f__f_byte_size;

  double* bs_emissivity_;
  int bs_emissivity__f_shapes[2];
  int bs_emissivity__f_byte_size;

  double* bs_user_emissivity_;
  int bs_user_emissivity__f_shapes[2];
  int bs_user_emissivity__f_byte_size;

  
};

// Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
class Brdf_Sup_Outputs : public Lidort_Structure {
public:
  // Allocating constructor
  Brdf_Sup_Outputs() : Lidort_Structure() {
    brdf_sup_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Sup_Outputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    brdf_sup_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Sup_Outputs() {
    if (owns_pointer)
      brdf_sup_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& bs_exactdb_brdfunc() const {
    return bs_exactdb_brdfunc_;
  }
  void bs_exactdb_brdfunc(const blitz::Array<double, 3>& bs_exactdb_brdfunc_in) {
    bs_exactdb_brdfunc_ = bs_exactdb_brdfunc_in;
  }
  
  const blitz::Array<double, 3>& bs_brdf_f_0() const {
    return bs_brdf_f_0_;
  }
  void bs_brdf_f_0(const blitz::Array<double, 3>& bs_brdf_f_0_in) {
    bs_brdf_f_0_ = bs_brdf_f_0_in;
  }
  
  const blitz::Array<double, 3>& bs_brdf_f() const {
    return bs_brdf_f_;
  }
  void bs_brdf_f(const blitz::Array<double, 3>& bs_brdf_f_in) {
    bs_brdf_f_ = bs_brdf_f_in;
  }
  
  const blitz::Array<double, 3>& bs_user_brdf_f_0() const {
    return bs_user_brdf_f_0_;
  }
  void bs_user_brdf_f_0(const blitz::Array<double, 3>& bs_user_brdf_f_0_in) {
    bs_user_brdf_f_0_ = bs_user_brdf_f_0_in;
  }
  
  const blitz::Array<double, 3>& bs_user_brdf_f() const {
    return bs_user_brdf_f_;
  }
  void bs_user_brdf_f(const blitz::Array<double, 3>& bs_user_brdf_f_in) {
    bs_user_brdf_f_ = bs_user_brdf_f_in;
  }
  
  const blitz::Array<double, 2>& bs_emissivity() const {
    return bs_emissivity_;
  }
  void bs_emissivity(const blitz::Array<double, 2>& bs_emissivity_in) {
    bs_emissivity_ = bs_emissivity_in;
  }
  
  const blitz::Array<double, 2>& bs_user_emissivity() const {
    return bs_user_emissivity_;
  }
  void bs_user_emissivity(const blitz::Array<double, 2>& bs_user_emissivity_in) {
    bs_user_emissivity_ = bs_user_emissivity_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Sup_Outputs:" << std::endl
      << "bs_exactdb_brdfunc: " << std::endl << bs_exactdb_brdfunc()  << std::endl
      << "       bs_brdf_f_0: " << std::endl << bs_brdf_f_0()  << std::endl
      << "         bs_brdf_f: " << std::endl << bs_brdf_f()  << std::endl
      << "  bs_user_brdf_f_0: " << std::endl << bs_user_brdf_f_0()  << std::endl
      << "    bs_user_brdf_f: " << std::endl << bs_user_brdf_f()  << std::endl
      << "     bs_emissivity: " << std::endl << bs_emissivity()  << std::endl
      << "bs_user_emissivity: " << std::endl << bs_user_emissivity()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_exactdb_brdfunc_",sizeof(*transfer_struct_c.bs_exactdb_brdfunc_),transfer_struct_c.bs_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_f_0_",sizeof(*transfer_struct_c.bs_brdf_f_0_),transfer_struct_c.bs_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_f_",sizeof(*transfer_struct_c.bs_brdf_f_),transfer_struct_c.bs_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_brdf_f_0_",sizeof(*transfer_struct_c.bs_user_brdf_f_0_),transfer_struct_c.bs_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_brdf_f_",sizeof(*transfer_struct_c.bs_user_brdf_f_),transfer_struct_c.bs_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_emissivity_",sizeof(*transfer_struct_c.bs_emissivity_),transfer_struct_c.bs_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_emissivity_",sizeof(*transfer_struct_c.bs_user_emissivity_),transfer_struct_c.bs_user_emissivity__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_exactdb_brdfunc_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.bs_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.bs_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.bs_exactdb_brdfunc__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_brdf_f_,
      blitz::shape(transfer_struct_c.bs_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_user_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_user_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_user_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_user_brdf_f_,
      blitz::shape(transfer_struct_c.bs_user_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_user_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_user_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_emissivity_,
      blitz::shape(transfer_struct_c.bs_emissivity__f_shapes[0],
                   transfer_struct_c.bs_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_user_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_user_emissivity_,
      blitz::shape(transfer_struct_c.bs_user_emissivity__f_shapes[0],
                   transfer_struct_c.bs_user_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_sup_outputs transfer_struct_c;

  blitz::Array<double, 3> bs_exactdb_brdfunc_;
  blitz::Array<double, 3> bs_brdf_f_0_;
  blitz::Array<double, 3> bs_brdf_f_;
  blitz::Array<double, 3> bs_user_brdf_f_0_;
  blitz::Array<double, 3> bs_user_brdf_f_;
  blitz::Array<double, 2> bs_emissivity_;
  blitz::Array<double, 2> bs_user_emissivity_;
  
};

// Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
extern "C" {
  void brdf_input_exception_handling_c_alloc_init(struct brdf_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_input_exception_handling_c_init_only(struct brdf_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_input_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_input_exception_handling_c_destroy(void **fortran_type_c);
  void brdf_input_exception_handling_bs_inputmessages_get(void **fortran_type_c, const int* bs_inputmessages_in_shape_1, const int* bs_inputmessages_in_len, const char* bs_inputmessages_in);
  void brdf_input_exception_handling_bs_inputactions_get(void **fortran_type_c, const int* bs_inputactions_in_shape_1, const int* bs_inputactions_in_len, const char* bs_inputactions_in);
  
}

struct brdf_input_exception_handling {
  int* bs_status_inputread_;
  int bs_status_inputread__f_byte_size;

  int* bs_ninputmessages_;
  int bs_ninputmessages__f_byte_size;

  
  int bs_inputmessages__f_shapes[1];
  int bs_inputmessages__f_len;

  
  int bs_inputactions__f_shapes[1];
  int bs_inputactions__f_len;

  
};

// Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
class Brdf_Input_Exception_Handling : public Lidort_Structure {
public:
  // Allocating constructor
  Brdf_Input_Exception_Handling() : Lidort_Structure() {
    brdf_input_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Input_Exception_Handling(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    brdf_input_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Input_Exception_Handling() {
    if (owns_pointer)
      brdf_input_exception_handling_c_destroy(&fortran_type_c);
  }

  const int& bs_status_inputread() const {
    return *transfer_struct_c.bs_status_inputread_;
  }
  void bs_status_inputread(const int& bs_status_inputread_in) {
    *transfer_struct_c.bs_status_inputread_ = bs_status_inputread_in;
  }
  
  const int& bs_ninputmessages() const {
    return *transfer_struct_c.bs_ninputmessages_;
  }
  void bs_ninputmessages(const int& bs_ninputmessages_in) {
    *transfer_struct_c.bs_ninputmessages_ = bs_ninputmessages_in;
  }
  
  const std::vector< std::string > bs_inputmessages() const {
    std::vector< std::string > bs_inputmessages_ret;
    blitz::Array<char, 2> bs_inputmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_inputmessages__f_shapes[0], transfer_struct_c.bs_inputmessages__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_input_exception_handling_bs_inputmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_inputmessages__f_shapes[0], &transfer_struct_c.bs_inputmessages__f_len, bs_inputmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_inputmessages_lcl.extent(0); dim_0_idx++)
      bs_inputmessages_ret.push_back( std::string(std::string(bs_inputmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_inputmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_inputmessages_ret;
  }
  
  const std::vector< std::string > bs_inputactions() const {
    std::vector< std::string > bs_inputactions_ret;
    blitz::Array<char, 2> bs_inputactions_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_inputactions__f_shapes[0], transfer_struct_c.bs_inputactions__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_input_exception_handling_bs_inputactions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_inputactions__f_shapes[0], &transfer_struct_c.bs_inputactions__f_len, bs_inputactions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_inputactions_lcl.extent(0); dim_0_idx++)
      bs_inputactions_ret.push_back( std::string(std::string(bs_inputactions_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_inputactions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_inputactions_ret;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Input_Exception_Handling:" << std::endl
      << "bs_status_inputread: " << bs_status_inputread()  << std::endl
      << "  bs_ninputmessages: " << bs_ninputmessages()  << std::endl
      << "   bs_inputmessages: " << std::endl;
    std::vector< std::string > bs_inputmessages_lcl = bs_inputmessages();
    for(unsigned int idx = 0; idx < bs_inputmessages_lcl.size(); idx++)
      if ( bs_inputmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_inputmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "    bs_inputactions: " << std::endl;
    std::vector< std::string > bs_inputactions_lcl = bs_inputactions();
    for(unsigned int idx = 0; idx < bs_inputactions_lcl.size(); idx++)
      if ( bs_inputactions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_inputactions_lcl[idx] << "\"" << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_status_inputread_",sizeof(*transfer_struct_c.bs_status_inputread_),transfer_struct_c.bs_status_inputread__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ninputmessages_",sizeof(*transfer_struct_c.bs_ninputmessages_),transfer_struct_c.bs_ninputmessages__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct brdf_input_exception_handling transfer_struct_c;

  
};

// Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
extern "C" {
  void lidort_fixed_lincontrol_c_alloc_init(struct lidort_fixed_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lincontrol_c_init_only(struct lidort_fixed_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lincontrol_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_lincontrol_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_lincontrol {
  int* ts_do_column_linearization_;
  int ts_do_column_linearization__f_byte_size;

  int* ts_do_profile_linearization_;
  int ts_do_profile_linearization__f_byte_size;

  int* ts_do_surface_linearization_;
  int ts_do_surface_linearization__f_byte_size;

  int* ts_do_sleave_wfs_;
  int ts_do_sleave_wfs__f_byte_size;

  int* ts_layer_vary_flag_;
  int ts_layer_vary_flag__f_shapes[1];
  int ts_layer_vary_flag__f_byte_size;

  int* ts_layer_vary_number_;
  int ts_layer_vary_number__f_shapes[1];
  int ts_layer_vary_number__f_byte_size;

  int* ts_n_totalcolumn_wfs_;
  int ts_n_totalcolumn_wfs__f_byte_size;

  int* ts_n_surface_wfs_;
  int ts_n_surface_wfs__f_byte_size;

  int* ts_n_sleave_wfs_;
  int ts_n_sleave_wfs__f_byte_size;

  
};

// Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
class Lidort_Fixed_Lincontrol : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Lincontrol() : Lidort_Structure() {
    lidort_fixed_lincontrol_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Lincontrol(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_lincontrol_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Lincontrol() {
    if (owns_pointer)
      lidort_fixed_lincontrol_c_destroy(&fortran_type_c);
  }

  const bool ts_do_column_linearization() const {
    return *transfer_struct_c.ts_do_column_linearization_ != 0;
  }
  void ts_do_column_linearization(const bool& ts_do_column_linearization_in) {
    *transfer_struct_c.ts_do_column_linearization_ = ts_do_column_linearization_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_profile_linearization() const {
    return *transfer_struct_c.ts_do_profile_linearization_ != 0;
  }
  void ts_do_profile_linearization(const bool& ts_do_profile_linearization_in) {
    *transfer_struct_c.ts_do_profile_linearization_ = ts_do_profile_linearization_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_surface_linearization() const {
    return *transfer_struct_c.ts_do_surface_linearization_ != 0;
  }
  void ts_do_surface_linearization(const bool& ts_do_surface_linearization_in) {
    *transfer_struct_c.ts_do_surface_linearization_ = ts_do_surface_linearization_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_sleave_wfs() const {
    return *transfer_struct_c.ts_do_sleave_wfs_ != 0;
  }
  void ts_do_sleave_wfs(const bool& ts_do_sleave_wfs_in) {
    *transfer_struct_c.ts_do_sleave_wfs_ = ts_do_sleave_wfs_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const blitz::Array<bool, 1> ts_layer_vary_flag() const {
    blitz::Array<bool,1> as_bool(ts_layer_vary_flag_.shape());
    as_bool = blitz::where(ts_layer_vary_flag_ != 0, true, false);
    return as_bool;
  }
  void ts_layer_vary_flag(const blitz::Array<bool, 1>& ts_layer_vary_flag_in) {
    blitz::Array<int,1> as_int(ts_layer_vary_flag_.shape());
    as_int = blitz::where(ts_layer_vary_flag_in == true, FORTRAN_TRUE_INT, 0);
    ts_layer_vary_flag_ = as_int;
  }
  
  const blitz::Array<int, 1>& ts_layer_vary_number() const {
    return ts_layer_vary_number_;
  }
  void ts_layer_vary_number(const blitz::Array<int, 1>& ts_layer_vary_number_in) {
    ts_layer_vary_number_ = ts_layer_vary_number_in;
  }
  
  const int& ts_n_totalcolumn_wfs() const {
    return *transfer_struct_c.ts_n_totalcolumn_wfs_;
  }
  void ts_n_totalcolumn_wfs(const int& ts_n_totalcolumn_wfs_in) {
    *transfer_struct_c.ts_n_totalcolumn_wfs_ = ts_n_totalcolumn_wfs_in;
  }
  
  const int& ts_n_surface_wfs() const {
    return *transfer_struct_c.ts_n_surface_wfs_;
  }
  void ts_n_surface_wfs(const int& ts_n_surface_wfs_in) {
    *transfer_struct_c.ts_n_surface_wfs_ = ts_n_surface_wfs_in;
  }
  
  const int& ts_n_sleave_wfs() const {
    return *transfer_struct_c.ts_n_sleave_wfs_;
  }
  void ts_n_sleave_wfs(const int& ts_n_sleave_wfs_in) {
    *transfer_struct_c.ts_n_sleave_wfs_ = ts_n_sleave_wfs_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Lincontrol:" << std::endl
      << " ts_do_column_linearization: " << ts_do_column_linearization()  << std::endl
      << "ts_do_profile_linearization: " << ts_do_profile_linearization()  << std::endl
      << "ts_do_surface_linearization: " << ts_do_surface_linearization()  << std::endl
      << "           ts_do_sleave_wfs: " << ts_do_sleave_wfs()  << std::endl
      << "         ts_layer_vary_flag: " << std::endl << ts_layer_vary_flag()  << std::endl
      << "       ts_layer_vary_number: " << std::endl << ts_layer_vary_number()  << std::endl
      << "       ts_n_totalcolumn_wfs: " << ts_n_totalcolumn_wfs()  << std::endl
      << "           ts_n_surface_wfs: " << ts_n_surface_wfs()  << std::endl
      << "            ts_n_sleave_wfs: " << ts_n_sleave_wfs()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_column_linearization_",sizeof(*transfer_struct_c.ts_do_column_linearization_),transfer_struct_c.ts_do_column_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_profile_linearization_",sizeof(*transfer_struct_c.ts_do_profile_linearization_),transfer_struct_c.ts_do_profile_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_linearization_",sizeof(*transfer_struct_c.ts_do_surface_linearization_),transfer_struct_c.ts_do_surface_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sleave_wfs_",sizeof(*transfer_struct_c.ts_do_sleave_wfs_),transfer_struct_c.ts_do_sleave_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_layer_vary_flag_",sizeof(*transfer_struct_c.ts_layer_vary_flag_),transfer_struct_c.ts_layer_vary_flag__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_layer_vary_number_",sizeof(*transfer_struct_c.ts_layer_vary_number_),transfer_struct_c.ts_layer_vary_number__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_totalcolumn_wfs_",sizeof(*transfer_struct_c.ts_n_totalcolumn_wfs_),transfer_struct_c.ts_n_totalcolumn_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_surface_wfs_",sizeof(*transfer_struct_c.ts_n_surface_wfs_),transfer_struct_c.ts_n_surface_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_sleave_wfs_",sizeof(*transfer_struct_c.ts_n_sleave_wfs_),transfer_struct_c.ts_n_sleave_wfs__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_layer_vary_flag_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_layer_vary_flag_,
      blitz::shape(transfer_struct_c.ts_layer_vary_flag__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_layer_vary_number_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_layer_vary_number_,
      blitz::shape(transfer_struct_c.ts_layer_vary_number__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_lincontrol transfer_struct_c;

  blitz::Array<int, 1> ts_layer_vary_flag_;
  blitz::Array<int, 1> ts_layer_vary_number_;
  
};

// Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
extern "C" {
  void lidort_fixed_linoptical_c_alloc_init(struct lidort_fixed_linoptical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_linoptical_c_init_only(struct lidort_fixed_linoptical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_linoptical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_linoptical_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_linoptical {
  double* ts_l_deltau_vert_input_;
  int ts_l_deltau_vert_input__f_shapes[3];
  int ts_l_deltau_vert_input__f_byte_size;

  double* ts_l_omega_total_input_;
  int ts_l_omega_total_input__f_shapes[3];
  int ts_l_omega_total_input__f_byte_size;

  double* ts_l_phasmoms_total_input_;
  int ts_l_phasmoms_total_input__f_shapes[4];
  int ts_l_phasmoms_total_input__f_byte_size;

  
};

// Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
class Lidort_Fixed_Linoptical : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Linoptical() : Lidort_Structure() {
    lidort_fixed_linoptical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Linoptical(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_linoptical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Linoptical() {
    if (owns_pointer)
      lidort_fixed_linoptical_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_l_deltau_vert_input() const {
    return ts_l_deltau_vert_input_;
  }
  void ts_l_deltau_vert_input(const blitz::Array<double, 3>& ts_l_deltau_vert_input_in) {
    ts_l_deltau_vert_input_ = ts_l_deltau_vert_input_in;
  }
  
  const blitz::Array<double, 3>& ts_l_omega_total_input() const {
    return ts_l_omega_total_input_;
  }
  void ts_l_omega_total_input(const blitz::Array<double, 3>& ts_l_omega_total_input_in) {
    ts_l_omega_total_input_ = ts_l_omega_total_input_in;
  }
  
  const blitz::Array<double, 4>& ts_l_phasmoms_total_input() const {
    return ts_l_phasmoms_total_input_;
  }
  void ts_l_phasmoms_total_input(const blitz::Array<double, 4>& ts_l_phasmoms_total_input_in) {
    ts_l_phasmoms_total_input_ = ts_l_phasmoms_total_input_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Linoptical:" << std::endl
      << "   ts_l_deltau_vert_input: " << std::endl << ts_l_deltau_vert_input()  << std::endl
      << "   ts_l_omega_total_input: " << std::endl << ts_l_omega_total_input()  << std::endl
      << "ts_l_phasmoms_total_input: " << std::endl << ts_l_phasmoms_total_input()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_l_deltau_vert_input_",sizeof(*transfer_struct_c.ts_l_deltau_vert_input_),transfer_struct_c.ts_l_deltau_vert_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_omega_total_input_",sizeof(*transfer_struct_c.ts_l_omega_total_input_),transfer_struct_c.ts_l_omega_total_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_phasmoms_total_input_",sizeof(*transfer_struct_c.ts_l_phasmoms_total_input_),transfer_struct_c.ts_l_phasmoms_total_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_l_deltau_vert_input_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_l_deltau_vert_input_,
      blitz::shape(transfer_struct_c.ts_l_deltau_vert_input__f_shapes[0],
                   transfer_struct_c.ts_l_deltau_vert_input__f_shapes[1],
                   transfer_struct_c.ts_l_deltau_vert_input__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_l_omega_total_input_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_l_omega_total_input_,
      blitz::shape(transfer_struct_c.ts_l_omega_total_input__f_shapes[0],
                   transfer_struct_c.ts_l_omega_total_input__f_shapes[1],
                   transfer_struct_c.ts_l_omega_total_input__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_l_phasmoms_total_input_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_l_phasmoms_total_input_,
      blitz::shape(transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[0],
                   transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[1],
                   transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[2],
                   transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_linoptical transfer_struct_c;

  blitz::Array<double, 3> ts_l_deltau_vert_input_;
  blitz::Array<double, 3> ts_l_omega_total_input_;
  blitz::Array<double, 4> ts_l_phasmoms_total_input_;
  
};

// Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
extern "C" {
  void lidort_fixed_lininputs_c_alloc_init(struct lidort_fixed_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lininputs_c_init_only(struct lidort_fixed_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lininputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_lininputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_lininputs {
  void* cont_;
  int cont__f_byte_size;

  void* optical_;
  int optical__f_byte_size;

  
};

// Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
class Lidort_Fixed_Lininputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Lininputs() : Lidort_Structure() {
    lidort_fixed_lininputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Lininputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_lininputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Lininputs() {
    if (owns_pointer)
      lidort_fixed_lininputs_c_destroy(&fortran_type_c);
  }

  Lidort_Fixed_Lincontrol& cont() {
    return *cont_;
  }
  const Lidort_Fixed_Lincontrol& cont() const {
    return *cont_;
  }
  void cont(Lidort_Fixed_Lincontrol& cont_in) {
    void* src_ptr = cont_in.fortran_type_ptr();
    void* dst_ptr = cont_->fortran_type_ptr();
    lidort_fixed_lincontrol_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Linoptical& optical() {
    return *optical_;
  }
  const Lidort_Fixed_Linoptical& optical() const {
    return *optical_;
  }
  void optical(Lidort_Fixed_Linoptical& optical_in) {
    void* src_ptr = optical_in.fortran_type_ptr();
    void* dst_ptr = optical_->fortran_type_ptr();
    lidort_fixed_linoptical_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Lininputs:" << std::endl
      << "   cont: " << cont()  << std::endl
      << "optical: " << optical()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    cont_.reset( new Lidort_Fixed_Lincontrol(transfer_struct_c.cont_) );
    optical_.reset( new Lidort_Fixed_Linoptical(transfer_struct_c.optical_) );
    
  }

  struct lidort_fixed_lininputs transfer_struct_c;

  boost::shared_ptr<Lidort_Fixed_Lincontrol> cont_;
  boost::shared_ptr<Lidort_Fixed_Linoptical> optical_;
  
};

// Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
extern "C" {
  void lidort_modified_lininputs_c_alloc_init(struct lidort_modified_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lininputs_c_init_only(struct lidort_modified_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lininputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_lininputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_lininputs {
  int* dummy_;
  int dummy__f_byte_size;

  
};

// Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
class Lidort_Modified_Lininputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Lininputs() : Lidort_Structure() {
    lidort_modified_lininputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Lininputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_lininputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Lininputs() {
    if (owns_pointer)
      lidort_modified_lininputs_c_destroy(&fortran_type_c);
  }

  const int& dummy() const {
    return *transfer_struct_c.dummy_;
  }
  void dummy(const int& dummy_in) {
    *transfer_struct_c.dummy_ = dummy_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Lininputs:" << std::endl
      << "dummy: " << dummy()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("dummy_",sizeof(*transfer_struct_c.dummy_),transfer_struct_c.dummy__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_lininputs transfer_struct_c;

  
};

// Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
extern "C" {
  void lidort_linatmos_c_alloc_init(struct lidort_linatmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linatmos_c_init_only(struct lidort_linatmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linatmos_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linatmos_c_destroy(void **fortran_type_c);
  
}

struct lidort_linatmos {
  double* ts_columnwf_;
  int ts_columnwf__f_shapes[5];
  int ts_columnwf__f_byte_size;

  double* ts_mint_columnwf_;
  int ts_mint_columnwf__f_shapes[5];
  int ts_mint_columnwf__f_byte_size;

  double* ts_flux_columnwf_;
  int ts_flux_columnwf__f_shapes[5];
  int ts_flux_columnwf__f_byte_size;

  double* ts_profilewf_;
  int ts_profilewf__f_shapes[6];
  int ts_profilewf__f_byte_size;

  double* ts_mint_profilewf_;
  int ts_mint_profilewf__f_shapes[6];
  int ts_mint_profilewf__f_byte_size;

  double* ts_flux_profilewf_;
  int ts_flux_profilewf__f_shapes[6];
  int ts_flux_profilewf__f_byte_size;

  
};

// Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
class Lidort_Linatmos : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linatmos() : Lidort_Structure() {
    lidort_linatmos_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linatmos(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linatmos_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linatmos() {
    if (owns_pointer)
      lidort_linatmos_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 5>& ts_columnwf() const {
    return ts_columnwf_;
  }
  void ts_columnwf(const blitz::Array<double, 5>& ts_columnwf_in) {
    ts_columnwf_ = ts_columnwf_in;
  }
  
  const blitz::Array<double, 5>& ts_mint_columnwf() const {
    return ts_mint_columnwf_;
  }
  void ts_mint_columnwf(const blitz::Array<double, 5>& ts_mint_columnwf_in) {
    ts_mint_columnwf_ = ts_mint_columnwf_in;
  }
  
  const blitz::Array<double, 5>& ts_flux_columnwf() const {
    return ts_flux_columnwf_;
  }
  void ts_flux_columnwf(const blitz::Array<double, 5>& ts_flux_columnwf_in) {
    ts_flux_columnwf_ = ts_flux_columnwf_in;
  }
  
  const blitz::Array<double, 6>& ts_profilewf() const {
    return ts_profilewf_;
  }
  void ts_profilewf(const blitz::Array<double, 6>& ts_profilewf_in) {
    ts_profilewf_ = ts_profilewf_in;
  }
  
  const blitz::Array<double, 6>& ts_mint_profilewf() const {
    return ts_mint_profilewf_;
  }
  void ts_mint_profilewf(const blitz::Array<double, 6>& ts_mint_profilewf_in) {
    ts_mint_profilewf_ = ts_mint_profilewf_in;
  }
  
  const blitz::Array<double, 6>& ts_flux_profilewf() const {
    return ts_flux_profilewf_;
  }
  void ts_flux_profilewf(const blitz::Array<double, 6>& ts_flux_profilewf_in) {
    ts_flux_profilewf_ = ts_flux_profilewf_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linatmos:" << std::endl
      << "      ts_columnwf: " << std::endl << ts_columnwf()  << std::endl
      << " ts_mint_columnwf: " << std::endl << ts_mint_columnwf()  << std::endl
      << " ts_flux_columnwf: " << std::endl << ts_flux_columnwf()  << std::endl
      << "     ts_profilewf: " << std::endl << ts_profilewf()  << std::endl
      << "ts_mint_profilewf: " << std::endl << ts_mint_profilewf()  << std::endl
      << "ts_flux_profilewf: " << std::endl << ts_flux_profilewf()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_",sizeof(*transfer_struct_c.ts_columnwf_),transfer_struct_c.ts_columnwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_mint_columnwf_",sizeof(*transfer_struct_c.ts_mint_columnwf_),transfer_struct_c.ts_mint_columnwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_columnwf_",sizeof(*transfer_struct_c.ts_flux_columnwf_),transfer_struct_c.ts_flux_columnwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_",sizeof(*transfer_struct_c.ts_profilewf_),transfer_struct_c.ts_profilewf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_mint_profilewf_",sizeof(*transfer_struct_c.ts_mint_profilewf_),transfer_struct_c.ts_mint_profilewf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_profilewf_",sizeof(*transfer_struct_c.ts_flux_profilewf_),transfer_struct_c.ts_flux_profilewf__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_columnwf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_columnwf_,
      blitz::shape(transfer_struct_c.ts_columnwf__f_shapes[0],
                   transfer_struct_c.ts_columnwf__f_shapes[1],
                   transfer_struct_c.ts_columnwf__f_shapes[2],
                   transfer_struct_c.ts_columnwf__f_shapes[3],
                   transfer_struct_c.ts_columnwf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_mint_columnwf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_mint_columnwf_,
      blitz::shape(transfer_struct_c.ts_mint_columnwf__f_shapes[0],
                   transfer_struct_c.ts_mint_columnwf__f_shapes[1],
                   transfer_struct_c.ts_mint_columnwf__f_shapes[2],
                   transfer_struct_c.ts_mint_columnwf__f_shapes[3],
                   transfer_struct_c.ts_mint_columnwf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_flux_columnwf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_flux_columnwf_,
      blitz::shape(transfer_struct_c.ts_flux_columnwf__f_shapes[0],
                   transfer_struct_c.ts_flux_columnwf__f_shapes[1],
                   transfer_struct_c.ts_flux_columnwf__f_shapes[2],
                   transfer_struct_c.ts_flux_columnwf__f_shapes[3],
                   transfer_struct_c.ts_flux_columnwf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_profilewf_.reference(blitz::Array<double, 6>(transfer_struct_c.ts_profilewf_,
      blitz::shape(transfer_struct_c.ts_profilewf__f_shapes[0],
                   transfer_struct_c.ts_profilewf__f_shapes[1],
                   transfer_struct_c.ts_profilewf__f_shapes[2],
                   transfer_struct_c.ts_profilewf__f_shapes[3],
                   transfer_struct_c.ts_profilewf__f_shapes[4],
                   transfer_struct_c.ts_profilewf__f_shapes[5]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<6>()));
    ts_mint_profilewf_.reference(blitz::Array<double, 6>(transfer_struct_c.ts_mint_profilewf_,
      blitz::shape(transfer_struct_c.ts_mint_profilewf__f_shapes[0],
                   transfer_struct_c.ts_mint_profilewf__f_shapes[1],
                   transfer_struct_c.ts_mint_profilewf__f_shapes[2],
                   transfer_struct_c.ts_mint_profilewf__f_shapes[3],
                   transfer_struct_c.ts_mint_profilewf__f_shapes[4],
                   transfer_struct_c.ts_mint_profilewf__f_shapes[5]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<6>()));
    ts_flux_profilewf_.reference(blitz::Array<double, 6>(transfer_struct_c.ts_flux_profilewf_,
      blitz::shape(transfer_struct_c.ts_flux_profilewf__f_shapes[0],
                   transfer_struct_c.ts_flux_profilewf__f_shapes[1],
                   transfer_struct_c.ts_flux_profilewf__f_shapes[2],
                   transfer_struct_c.ts_flux_profilewf__f_shapes[3],
                   transfer_struct_c.ts_flux_profilewf__f_shapes[4],
                   transfer_struct_c.ts_flux_profilewf__f_shapes[5]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<6>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linatmos transfer_struct_c;

  blitz::Array<double, 5> ts_columnwf_;
  blitz::Array<double, 5> ts_mint_columnwf_;
  blitz::Array<double, 5> ts_flux_columnwf_;
  blitz::Array<double, 6> ts_profilewf_;
  blitz::Array<double, 6> ts_mint_profilewf_;
  blitz::Array<double, 6> ts_flux_profilewf_;
  
};

// Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
extern "C" {
  void lidort_linsurf_c_alloc_init(struct lidort_linsurf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsurf_c_init_only(struct lidort_linsurf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsurf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsurf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsurf {
  double* ts_surfacewf_;
  int ts_surfacewf__f_shapes[5];
  int ts_surfacewf__f_byte_size;

  double* ts_mint_surfacewf_;
  int ts_mint_surfacewf__f_shapes[5];
  int ts_mint_surfacewf__f_byte_size;

  double* ts_flux_surfacewf_;
  int ts_flux_surfacewf__f_shapes[5];
  int ts_flux_surfacewf__f_byte_size;

  
};

// Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
class Lidort_Linsurf : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsurf() : Lidort_Structure() {
    lidort_linsurf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsurf(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsurf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsurf() {
    if (owns_pointer)
      lidort_linsurf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 5>& ts_surfacewf() const {
    return ts_surfacewf_;
  }
  void ts_surfacewf(const blitz::Array<double, 5>& ts_surfacewf_in) {
    ts_surfacewf_ = ts_surfacewf_in;
  }
  
  const blitz::Array<double, 5>& ts_mint_surfacewf() const {
    return ts_mint_surfacewf_;
  }
  void ts_mint_surfacewf(const blitz::Array<double, 5>& ts_mint_surfacewf_in) {
    ts_mint_surfacewf_ = ts_mint_surfacewf_in;
  }
  
  const blitz::Array<double, 5>& ts_flux_surfacewf() const {
    return ts_flux_surfacewf_;
  }
  void ts_flux_surfacewf(const blitz::Array<double, 5>& ts_flux_surfacewf_in) {
    ts_flux_surfacewf_ = ts_flux_surfacewf_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsurf:" << std::endl
      << "     ts_surfacewf: " << std::endl << ts_surfacewf()  << std::endl
      << "ts_mint_surfacewf: " << std::endl << ts_mint_surfacewf()  << std::endl
      << "ts_flux_surfacewf: " << std::endl << ts_flux_surfacewf()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_surfacewf_",sizeof(*transfer_struct_c.ts_surfacewf_),transfer_struct_c.ts_surfacewf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_mint_surfacewf_",sizeof(*transfer_struct_c.ts_mint_surfacewf_),transfer_struct_c.ts_mint_surfacewf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_surfacewf_",sizeof(*transfer_struct_c.ts_flux_surfacewf_),transfer_struct_c.ts_flux_surfacewf__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_surfacewf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_surfacewf_,
      blitz::shape(transfer_struct_c.ts_surfacewf__f_shapes[0],
                   transfer_struct_c.ts_surfacewf__f_shapes[1],
                   transfer_struct_c.ts_surfacewf__f_shapes[2],
                   transfer_struct_c.ts_surfacewf__f_shapes[3],
                   transfer_struct_c.ts_surfacewf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_mint_surfacewf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_mint_surfacewf_,
      blitz::shape(transfer_struct_c.ts_mint_surfacewf__f_shapes[0],
                   transfer_struct_c.ts_mint_surfacewf__f_shapes[1],
                   transfer_struct_c.ts_mint_surfacewf__f_shapes[2],
                   transfer_struct_c.ts_mint_surfacewf__f_shapes[3],
                   transfer_struct_c.ts_mint_surfacewf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_flux_surfacewf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_flux_surfacewf_,
      blitz::shape(transfer_struct_c.ts_flux_surfacewf__f_shapes[0],
                   transfer_struct_c.ts_flux_surfacewf__f_shapes[1],
                   transfer_struct_c.ts_flux_surfacewf__f_shapes[2],
                   transfer_struct_c.ts_flux_surfacewf__f_shapes[3],
                   transfer_struct_c.ts_flux_surfacewf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsurf transfer_struct_c;

  blitz::Array<double, 5> ts_surfacewf_;
  blitz::Array<double, 5> ts_mint_surfacewf_;
  blitz::Array<double, 5> ts_flux_surfacewf_;
  
};

// Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
extern "C" {
  void lidort_linoutputs_c_alloc_init(struct lidort_linoutputs *transfer_struct_c, void **fortran_type_c);
  void lidort_linoutputs_c_init_only(struct lidort_linoutputs *transfer_struct_c, void **fortran_type_c);
  void lidort_linoutputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linoutputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_linoutputs {
  void* atmos_;
  int atmos__f_byte_size;

  void* surf_;
  int surf__f_byte_size;

  
};

// Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
class Lidort_Linoutputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linoutputs() : Lidort_Structure() {
    lidort_linoutputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linoutputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linoutputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linoutputs() {
    if (owns_pointer)
      lidort_linoutputs_c_destroy(&fortran_type_c);
  }

  Lidort_Linatmos& atmos() {
    return *atmos_;
  }
  const Lidort_Linatmos& atmos() const {
    return *atmos_;
  }
  void atmos(Lidort_Linatmos& atmos_in) {
    void* src_ptr = atmos_in.fortran_type_ptr();
    void* dst_ptr = atmos_->fortran_type_ptr();
    lidort_linatmos_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Linsurf& surf() {
    return *surf_;
  }
  const Lidort_Linsurf& surf() const {
    return *surf_;
  }
  void surf(Lidort_Linsurf& surf_in) {
    void* src_ptr = surf_in.fortran_type_ptr();
    void* dst_ptr = surf_->fortran_type_ptr();
    lidort_linsurf_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linoutputs:" << std::endl
      << "atmos: " << atmos()  << std::endl
      << " surf: " << surf()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    atmos_.reset( new Lidort_Linatmos(transfer_struct_c.atmos_) );
    surf_.reset( new Lidort_Linsurf(transfer_struct_c.surf_) );
    
  }

  struct lidort_linoutputs transfer_struct_c;

  boost::shared_ptr<Lidort_Linatmos> atmos_;
  boost::shared_ptr<Lidort_Linsurf> surf_;
  
};

// Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
extern "C" {
  void lidort_linsup_brdf_c_alloc_init(struct lidort_linsup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_brdf_c_init_only(struct lidort_linsup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_brdf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_brdf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_brdf {
  double* ts_ls_exactdb_brdfunc_;
  int ts_ls_exactdb_brdfunc__f_shapes[4];
  int ts_ls_exactdb_brdfunc__f_byte_size;

  double* ts_ls_brdf_f_0_;
  int ts_ls_brdf_f_0__f_shapes[4];
  int ts_ls_brdf_f_0__f_byte_size;

  double* ts_ls_brdf_f_;
  int ts_ls_brdf_f__f_shapes[4];
  int ts_ls_brdf_f__f_byte_size;

  double* ts_ls_user_brdf_f_0_;
  int ts_ls_user_brdf_f_0__f_shapes[4];
  int ts_ls_user_brdf_f_0__f_byte_size;

  double* ts_ls_user_brdf_f_;
  int ts_ls_user_brdf_f__f_shapes[4];
  int ts_ls_user_brdf_f__f_byte_size;

  double* ts_ls_emissivity_;
  int ts_ls_emissivity__f_shapes[3];
  int ts_ls_emissivity__f_byte_size;

  double* ts_ls_user_emissivity_;
  int ts_ls_user_emissivity__f_shapes[3];
  int ts_ls_user_emissivity__f_byte_size;

  
};

// Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
class Lidort_Linsup_Brdf : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Brdf() : Lidort_Structure() {
    lidort_linsup_brdf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Brdf(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_brdf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Brdf() {
    if (owns_pointer)
      lidort_linsup_brdf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc() const {
    return ts_ls_exactdb_brdfunc_;
  }
  void ts_ls_exactdb_brdfunc(const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc_in) {
    ts_ls_exactdb_brdfunc_ = ts_ls_exactdb_brdfunc_in;
  }
  
  const blitz::Array<double, 4>& ts_ls_brdf_f_0() const {
    return ts_ls_brdf_f_0_;
  }
  void ts_ls_brdf_f_0(const blitz::Array<double, 4>& ts_ls_brdf_f_0_in) {
    ts_ls_brdf_f_0_ = ts_ls_brdf_f_0_in;
  }
  
  const blitz::Array<double, 4>& ts_ls_brdf_f() const {
    return ts_ls_brdf_f_;
  }
  void ts_ls_brdf_f(const blitz::Array<double, 4>& ts_ls_brdf_f_in) {
    ts_ls_brdf_f_ = ts_ls_brdf_f_in;
  }
  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f_0() const {
    return ts_ls_user_brdf_f_0_;
  }
  void ts_ls_user_brdf_f_0(const blitz::Array<double, 4>& ts_ls_user_brdf_f_0_in) {
    ts_ls_user_brdf_f_0_ = ts_ls_user_brdf_f_0_in;
  }
  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f() const {
    return ts_ls_user_brdf_f_;
  }
  void ts_ls_user_brdf_f(const blitz::Array<double, 4>& ts_ls_user_brdf_f_in) {
    ts_ls_user_brdf_f_ = ts_ls_user_brdf_f_in;
  }
  
  const blitz::Array<double, 3>& ts_ls_emissivity() const {
    return ts_ls_emissivity_;
  }
  void ts_ls_emissivity(const blitz::Array<double, 3>& ts_ls_emissivity_in) {
    ts_ls_emissivity_ = ts_ls_emissivity_in;
  }
  
  const blitz::Array<double, 3>& ts_ls_user_emissivity() const {
    return ts_ls_user_emissivity_;
  }
  void ts_ls_user_emissivity(const blitz::Array<double, 3>& ts_ls_user_emissivity_in) {
    ts_ls_user_emissivity_ = ts_ls_user_emissivity_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Brdf:" << std::endl
      << "ts_ls_exactdb_brdfunc: " << std::endl << ts_ls_exactdb_brdfunc()  << std::endl
      << "       ts_ls_brdf_f_0: " << std::endl << ts_ls_brdf_f_0()  << std::endl
      << "         ts_ls_brdf_f: " << std::endl << ts_ls_brdf_f()  << std::endl
      << "  ts_ls_user_brdf_f_0: " << std::endl << ts_ls_user_brdf_f_0()  << std::endl
      << "    ts_ls_user_brdf_f: " << std::endl << ts_ls_user_brdf_f()  << std::endl
      << "     ts_ls_emissivity: " << std::endl << ts_ls_emissivity()  << std::endl
      << "ts_ls_user_emissivity: " << std::endl << ts_ls_user_emissivity()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_ls_exactdb_brdfunc_",sizeof(*transfer_struct_c.ts_ls_exactdb_brdfunc_),transfer_struct_c.ts_ls_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_brdf_f_0_",sizeof(*transfer_struct_c.ts_ls_brdf_f_0_),transfer_struct_c.ts_ls_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_brdf_f_",sizeof(*transfer_struct_c.ts_ls_brdf_f_),transfer_struct_c.ts_ls_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_brdf_f_0_",sizeof(*transfer_struct_c.ts_ls_user_brdf_f_0_),transfer_struct_c.ts_ls_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_brdf_f_",sizeof(*transfer_struct_c.ts_ls_user_brdf_f_),transfer_struct_c.ts_ls_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_emissivity_",sizeof(*transfer_struct_c.ts_ls_emissivity_),transfer_struct_c.ts_ls_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_emissivity_",sizeof(*transfer_struct_c.ts_ls_user_emissivity_),transfer_struct_c.ts_ls_user_emissivity__f_byte_size);
    
  }
  
  void copy_from_sup(Brdf_Linsup_Outputs& supp_obj) { 
    // This is only safe for LIDORT 3.6 because the BRDF and LIDORT
    // sup structures are the exact same structure, but with different
    // names. This MUST be reevaluated on future LIDORT versions
    void* sup_ptr = supp_obj.fortran_type_ptr();
    lidort_linsup_brdf_c_copy(&sup_ptr, &fortran_type_c);
  }

private:
  void link_blitz_arrays() {
    ts_ls_exactdb_brdfunc_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[2],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_ls_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[2],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_brdf_f_,
      blitz::shape(transfer_struct_c.ts_ls_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[2],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_user_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[2],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_user_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_user_brdf_f_,
      blitz::shape(transfer_struct_c.ts_ls_user_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[2],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_emissivity_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_ls_emissivity_,
      blitz::shape(transfer_struct_c.ts_ls_emissivity__f_shapes[0],
                   transfer_struct_c.ts_ls_emissivity__f_shapes[1],
                   transfer_struct_c.ts_ls_emissivity__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_ls_user_emissivity_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_ls_user_emissivity_,
      blitz::shape(transfer_struct_c.ts_ls_user_emissivity__f_shapes[0],
                   transfer_struct_c.ts_ls_user_emissivity__f_shapes[1],
                   transfer_struct_c.ts_ls_user_emissivity__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_brdf transfer_struct_c;

  blitz::Array<double, 4> ts_ls_exactdb_brdfunc_;
  blitz::Array<double, 4> ts_ls_brdf_f_0_;
  blitz::Array<double, 4> ts_ls_brdf_f_;
  blitz::Array<double, 4> ts_ls_user_brdf_f_0_;
  blitz::Array<double, 4> ts_ls_user_brdf_f_;
  blitz::Array<double, 3> ts_ls_emissivity_;
  blitz::Array<double, 3> ts_ls_user_emissivity_;
  
};

// Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
extern "C" {
  void lidort_linsup_ss_atmos_c_alloc_init(struct lidort_linsup_ss_atmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_atmos_c_init_only(struct lidort_linsup_ss_atmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_atmos_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_atmos_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss_atmos {
  double* ts_columnwf_ss_;
  int ts_columnwf_ss__f_shapes[4];
  int ts_columnwf_ss__f_byte_size;

  double* ts_columnwf_db_;
  int ts_columnwf_db__f_shapes[3];
  int ts_columnwf_db__f_byte_size;

  double* ts_profilewf_ss_;
  int ts_profilewf_ss__f_shapes[5];
  int ts_profilewf_ss__f_byte_size;

  double* ts_profilewf_db_;
  int ts_profilewf_db__f_shapes[4];
  int ts_profilewf_db__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
class Lidort_Linsup_Ss_Atmos : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Ss_Atmos() : Lidort_Structure() {
    lidort_linsup_ss_atmos_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss_Atmos(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_ss_atmos_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss_Atmos() {
    if (owns_pointer)
      lidort_linsup_ss_atmos_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_columnwf_ss() const {
    return ts_columnwf_ss_;
  }
  void ts_columnwf_ss(const blitz::Array<double, 4>& ts_columnwf_ss_in) {
    ts_columnwf_ss_ = ts_columnwf_ss_in;
  }
  
  const blitz::Array<double, 3>& ts_columnwf_db() const {
    return ts_columnwf_db_;
  }
  void ts_columnwf_db(const blitz::Array<double, 3>& ts_columnwf_db_in) {
    ts_columnwf_db_ = ts_columnwf_db_in;
  }
  
  const blitz::Array<double, 5>& ts_profilewf_ss() const {
    return ts_profilewf_ss_;
  }
  void ts_profilewf_ss(const blitz::Array<double, 5>& ts_profilewf_ss_in) {
    ts_profilewf_ss_ = ts_profilewf_ss_in;
  }
  
  const blitz::Array<double, 4>& ts_profilewf_db() const {
    return ts_profilewf_db_;
  }
  void ts_profilewf_db(const blitz::Array<double, 4>& ts_profilewf_db_in) {
    ts_profilewf_db_ = ts_profilewf_db_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss_Atmos:" << std::endl
      << " ts_columnwf_ss: " << std::endl << ts_columnwf_ss()  << std::endl
      << " ts_columnwf_db: " << std::endl << ts_columnwf_db()  << std::endl
      << "ts_profilewf_ss: " << std::endl << ts_profilewf_ss()  << std::endl
      << "ts_profilewf_db: " << std::endl << ts_profilewf_db()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_ss_",sizeof(*transfer_struct_c.ts_columnwf_ss_),transfer_struct_c.ts_columnwf_ss__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_db_",sizeof(*transfer_struct_c.ts_columnwf_db_),transfer_struct_c.ts_columnwf_db__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_ss_",sizeof(*transfer_struct_c.ts_profilewf_ss_),transfer_struct_c.ts_profilewf_ss__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_db_",sizeof(*transfer_struct_c.ts_profilewf_db_),transfer_struct_c.ts_profilewf_db__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_columnwf_ss_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_columnwf_ss_,
      blitz::shape(transfer_struct_c.ts_columnwf_ss__f_shapes[0],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[1],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[2],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_columnwf_db_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_columnwf_db_,
      blitz::shape(transfer_struct_c.ts_columnwf_db__f_shapes[0],
                   transfer_struct_c.ts_columnwf_db__f_shapes[1],
                   transfer_struct_c.ts_columnwf_db__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_profilewf_ss_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_profilewf_ss_,
      blitz::shape(transfer_struct_c.ts_profilewf_ss__f_shapes[0],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[1],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[2],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[3],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_profilewf_db_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_profilewf_db_,
      blitz::shape(transfer_struct_c.ts_profilewf_db__f_shapes[0],
                   transfer_struct_c.ts_profilewf_db__f_shapes[1],
                   transfer_struct_c.ts_profilewf_db__f_shapes[2],
                   transfer_struct_c.ts_profilewf_db__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_ss_atmos transfer_struct_c;

  blitz::Array<double, 4> ts_columnwf_ss_;
  blitz::Array<double, 3> ts_columnwf_db_;
  blitz::Array<double, 5> ts_profilewf_ss_;
  blitz::Array<double, 4> ts_profilewf_db_;
  
};

// Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
extern "C" {
  void lidort_linsup_ss_surf_c_alloc_init(struct lidort_linsup_ss_surf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_surf_c_init_only(struct lidort_linsup_ss_surf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_surf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_surf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss_surf {
  double* ts_surfacewf_db_;
  int ts_surfacewf_db__f_shapes[3];
  int ts_surfacewf_db__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
class Lidort_Linsup_Ss_Surf : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Ss_Surf() : Lidort_Structure() {
    lidort_linsup_ss_surf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss_Surf(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_ss_surf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss_Surf() {
    if (owns_pointer)
      lidort_linsup_ss_surf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_surfacewf_db() const {
    return ts_surfacewf_db_;
  }
  void ts_surfacewf_db(const blitz::Array<double, 3>& ts_surfacewf_db_in) {
    ts_surfacewf_db_ = ts_surfacewf_db_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss_Surf:" << std::endl
      << "ts_surfacewf_db: " << std::endl << ts_surfacewf_db()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_surfacewf_db_",sizeof(*transfer_struct_c.ts_surfacewf_db_),transfer_struct_c.ts_surfacewf_db__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_surfacewf_db_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_surfacewf_db_,
      blitz::shape(transfer_struct_c.ts_surfacewf_db__f_shapes[0],
                   transfer_struct_c.ts_surfacewf_db__f_shapes[1],
                   transfer_struct_c.ts_surfacewf_db__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_ss_surf transfer_struct_c;

  blitz::Array<double, 3> ts_surfacewf_db_;
  
};

// Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
extern "C" {
  void lidort_linsup_ss_c_alloc_init(struct lidort_linsup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_c_init_only(struct lidort_linsup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss {
  void* atmos_;
  int atmos__f_byte_size;

  void* surf_;
  int surf__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
class Lidort_Linsup_Ss : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Ss() : Lidort_Structure() {
    lidort_linsup_ss_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_ss_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss() {
    if (owns_pointer)
      lidort_linsup_ss_c_destroy(&fortran_type_c);
  }

  Lidort_Linsup_Ss_Atmos& atmos() {
    return *atmos_;
  }
  const Lidort_Linsup_Ss_Atmos& atmos() const {
    return *atmos_;
  }
  void atmos(Lidort_Linsup_Ss_Atmos& atmos_in) {
    void* src_ptr = atmos_in.fortran_type_ptr();
    void* dst_ptr = atmos_->fortran_type_ptr();
    lidort_linsup_ss_atmos_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Linsup_Ss_Surf& surf() {
    return *surf_;
  }
  const Lidort_Linsup_Ss_Surf& surf() const {
    return *surf_;
  }
  void surf(Lidort_Linsup_Ss_Surf& surf_in) {
    void* src_ptr = surf_in.fortran_type_ptr();
    void* dst_ptr = surf_->fortran_type_ptr();
    lidort_linsup_ss_surf_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss:" << std::endl
      << "atmos: " << atmos()  << std::endl
      << " surf: " << surf()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    atmos_.reset( new Lidort_Linsup_Ss_Atmos(transfer_struct_c.atmos_) );
    surf_.reset( new Lidort_Linsup_Ss_Surf(transfer_struct_c.surf_) );
    
  }

  struct lidort_linsup_ss transfer_struct_c;

  boost::shared_ptr<Lidort_Linsup_Ss_Atmos> atmos_;
  boost::shared_ptr<Lidort_Linsup_Ss_Surf> surf_;
  
};

// Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
extern "C" {
  void lidort_linsup_sleave_c_alloc_init(struct lidort_linsup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_sleave_c_init_only(struct lidort_linsup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_sleave_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_sleave_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_sleave {
  double* ts_lssl_slterm_isotropic_;
  int ts_lssl_slterm_isotropic__f_shapes[2];
  int ts_lssl_slterm_isotropic__f_byte_size;

  double* ts_lssl_slterm_userangles_;
  int ts_lssl_slterm_userangles__f_shapes[4];
  int ts_lssl_slterm_userangles__f_byte_size;

  double* ts_lssl_slterm_f_0_;
  int ts_lssl_slterm_f_0__f_shapes[4];
  int ts_lssl_slterm_f_0__f_byte_size;

  double* ts_lssl_user_slterm_f_0_;
  int ts_lssl_user_slterm_f_0__f_shapes[4];
  int ts_lssl_user_slterm_f_0__f_byte_size;

  
};

// Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
class Lidort_Linsup_Sleave : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Sleave() : Lidort_Structure() {
    lidort_linsup_sleave_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Sleave(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_sleave_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Sleave() {
    if (owns_pointer)
      lidort_linsup_sleave_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 2>& ts_lssl_slterm_isotropic() const {
    return ts_lssl_slterm_isotropic_;
  }
  void ts_lssl_slterm_isotropic(const blitz::Array<double, 2>& ts_lssl_slterm_isotropic_in) {
    ts_lssl_slterm_isotropic_ = ts_lssl_slterm_isotropic_in;
  }
  
  const blitz::Array<double, 4>& ts_lssl_slterm_userangles() const {
    return ts_lssl_slterm_userangles_;
  }
  void ts_lssl_slterm_userangles(const blitz::Array<double, 4>& ts_lssl_slterm_userangles_in) {
    ts_lssl_slterm_userangles_ = ts_lssl_slterm_userangles_in;
  }
  
  const blitz::Array<double, 4>& ts_lssl_slterm_f_0() const {
    return ts_lssl_slterm_f_0_;
  }
  void ts_lssl_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_slterm_f_0_in) {
    ts_lssl_slterm_f_0_ = ts_lssl_slterm_f_0_in;
  }
  
  const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0() const {
    return ts_lssl_user_slterm_f_0_;
  }
  void ts_lssl_user_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0_in) {
    ts_lssl_user_slterm_f_0_ = ts_lssl_user_slterm_f_0_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Sleave:" << std::endl
      << " ts_lssl_slterm_isotropic: " << std::endl << ts_lssl_slterm_isotropic()  << std::endl
      << "ts_lssl_slterm_userangles: " << std::endl << ts_lssl_slterm_userangles()  << std::endl
      << "       ts_lssl_slterm_f_0: " << std::endl << ts_lssl_slterm_f_0()  << std::endl
      << "  ts_lssl_user_slterm_f_0: " << std::endl << ts_lssl_user_slterm_f_0()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_isotropic_",sizeof(*transfer_struct_c.ts_lssl_slterm_isotropic_),transfer_struct_c.ts_lssl_slterm_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_userangles_",sizeof(*transfer_struct_c.ts_lssl_slterm_userangles_),transfer_struct_c.ts_lssl_slterm_userangles__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_f_0_",sizeof(*transfer_struct_c.ts_lssl_slterm_f_0_),transfer_struct_c.ts_lssl_slterm_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_user_slterm_f_0_",sizeof(*transfer_struct_c.ts_lssl_user_slterm_f_0_),transfer_struct_c.ts_lssl_user_slterm_f_0__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_lssl_slterm_isotropic_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_lssl_slterm_isotropic_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_isotropic__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_isotropic__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_lssl_slterm_userangles_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_slterm_userangles_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[1],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[2],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_lssl_slterm_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[2],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_lssl_user_slterm_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_user_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[2],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_sleave transfer_struct_c;

  blitz::Array<double, 2> ts_lssl_slterm_isotropic_;
  blitz::Array<double, 4> ts_lssl_slterm_userangles_;
  blitz::Array<double, 4> ts_lssl_slterm_f_0_;
  blitz::Array<double, 4> ts_lssl_user_slterm_f_0_;
  
};

// Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
extern "C" {
  void lidort_linsup_inout_c_alloc_init(struct lidort_linsup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_inout_c_init_only(struct lidort_linsup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_inout_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_inout_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_inout {
  void* brdf_;
  int brdf__f_byte_size;

  void* ss_;
  int ss__f_byte_size;

  void* sleave_;
  int sleave__f_byte_size;

  
};

// Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
class Lidort_Linsup_Inout : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Inout() : Lidort_Structure() {
    lidort_linsup_inout_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Inout(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_linsup_inout_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Inout() {
    if (owns_pointer)
      lidort_linsup_inout_c_destroy(&fortran_type_c);
  }

  Lidort_Linsup_Brdf& brdf() {
    return *brdf_;
  }
  const Lidort_Linsup_Brdf& brdf() const {
    return *brdf_;
  }
  void brdf(Lidort_Linsup_Brdf& brdf_in) {
    void* src_ptr = brdf_in.fortran_type_ptr();
    void* dst_ptr = brdf_->fortran_type_ptr();
    lidort_linsup_brdf_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Linsup_Ss& ss() {
    return *ss_;
  }
  const Lidort_Linsup_Ss& ss() const {
    return *ss_;
  }
  void ss(Lidort_Linsup_Ss& ss_in) {
    void* src_ptr = ss_in.fortran_type_ptr();
    void* dst_ptr = ss_->fortran_type_ptr();
    lidort_linsup_ss_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Linsup_Sleave& sleave() {
    return *sleave_;
  }
  const Lidort_Linsup_Sleave& sleave() const {
    return *sleave_;
  }
  void sleave(Lidort_Linsup_Sleave& sleave_in) {
    void* src_ptr = sleave_in.fortran_type_ptr();
    void* dst_ptr = sleave_->fortran_type_ptr();
    lidort_linsup_sleave_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Inout:" << std::endl
      << "  brdf: " << brdf()  << std::endl
      << "    ss: " << ss()  << std::endl
      << "sleave: " << sleave()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    brdf_.reset( new Lidort_Linsup_Brdf(transfer_struct_c.brdf_) );
    ss_.reset( new Lidort_Linsup_Ss(transfer_struct_c.ss_) );
    sleave_.reset( new Lidort_Linsup_Sleave(transfer_struct_c.sleave_) );
    
  }

  struct lidort_linsup_inout transfer_struct_c;

  boost::shared_ptr<Lidort_Linsup_Brdf> brdf_;
  boost::shared_ptr<Lidort_Linsup_Ss> ss_;
  boost::shared_ptr<Lidort_Linsup_Sleave> sleave_;
  
};

// Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
extern "C" {
  void lidort_main_outputs_c_alloc_init(struct lidort_main_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_main_outputs_c_init_only(struct lidort_main_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_main_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_main_outputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_main_outputs {
  double* ts_intensity_;
  int ts_intensity__f_shapes[4];
  int ts_intensity__f_byte_size;

  double* ts_mean_intensity_;
  int ts_mean_intensity__f_shapes[4];
  int ts_mean_intensity__f_byte_size;

  double* ts_flux_integral_;
  int ts_flux_integral__f_shapes[4];
  int ts_flux_integral__f_byte_size;

  double* ts_dnflux_direct_;
  int ts_dnflux_direct__f_shapes[3];
  int ts_dnflux_direct__f_byte_size;

  double* ts_dnmean_direct_;
  int ts_dnmean_direct__f_shapes[3];
  int ts_dnmean_direct__f_byte_size;

  int* ts_fourier_saved_;
  int ts_fourier_saved__f_shapes[2];
  int ts_fourier_saved__f_byte_size;

  int* ts_n_geometries_;
  int ts_n_geometries__f_byte_size;

  
};

// Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
class Lidort_Main_Outputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Main_Outputs() : Lidort_Structure() {
    lidort_main_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Main_Outputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_main_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Main_Outputs() {
    if (owns_pointer)
      lidort_main_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_intensity() const {
    return ts_intensity_;
  }
  void ts_intensity(const blitz::Array<double, 4>& ts_intensity_in) {
    ts_intensity_ = ts_intensity_in;
  }
  
  const blitz::Array<double, 4>& ts_mean_intensity() const {
    return ts_mean_intensity_;
  }
  void ts_mean_intensity(const blitz::Array<double, 4>& ts_mean_intensity_in) {
    ts_mean_intensity_ = ts_mean_intensity_in;
  }
  
  const blitz::Array<double, 4>& ts_flux_integral() const {
    return ts_flux_integral_;
  }
  void ts_flux_integral(const blitz::Array<double, 4>& ts_flux_integral_in) {
    ts_flux_integral_ = ts_flux_integral_in;
  }
  
  const blitz::Array<double, 3>& ts_dnflux_direct() const {
    return ts_dnflux_direct_;
  }
  void ts_dnflux_direct(const blitz::Array<double, 3>& ts_dnflux_direct_in) {
    ts_dnflux_direct_ = ts_dnflux_direct_in;
  }
  
  const blitz::Array<double, 3>& ts_dnmean_direct() const {
    return ts_dnmean_direct_;
  }
  void ts_dnmean_direct(const blitz::Array<double, 3>& ts_dnmean_direct_in) {
    ts_dnmean_direct_ = ts_dnmean_direct_in;
  }
  
  const blitz::Array<int, 2>& ts_fourier_saved() const {
    return ts_fourier_saved_;
  }
  void ts_fourier_saved(const blitz::Array<int, 2>& ts_fourier_saved_in) {
    ts_fourier_saved_ = ts_fourier_saved_in;
  }
  
  const int& ts_n_geometries() const {
    return *transfer_struct_c.ts_n_geometries_;
  }
  void ts_n_geometries(const int& ts_n_geometries_in) {
    *transfer_struct_c.ts_n_geometries_ = ts_n_geometries_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Main_Outputs:" << std::endl
      << "     ts_intensity: " << std::endl << ts_intensity()  << std::endl
      << "ts_mean_intensity: " << std::endl << ts_mean_intensity()  << std::endl
      << " ts_flux_integral: " << std::endl << ts_flux_integral()  << std::endl
      << " ts_dnflux_direct: " << std::endl << ts_dnflux_direct()  << std::endl
      << " ts_dnmean_direct: " << std::endl << ts_dnmean_direct()  << std::endl
      << " ts_fourier_saved: " << std::endl << ts_fourier_saved()  << std::endl
      << "  ts_n_geometries: " << ts_n_geometries()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_intensity_",sizeof(*transfer_struct_c.ts_intensity_),transfer_struct_c.ts_intensity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_mean_intensity_",sizeof(*transfer_struct_c.ts_mean_intensity_),transfer_struct_c.ts_mean_intensity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_integral_",sizeof(*transfer_struct_c.ts_flux_integral_),transfer_struct_c.ts_flux_integral__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnflux_direct_",sizeof(*transfer_struct_c.ts_dnflux_direct_),transfer_struct_c.ts_dnflux_direct__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnmean_direct_",sizeof(*transfer_struct_c.ts_dnmean_direct_),transfer_struct_c.ts_dnmean_direct__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_fourier_saved_",sizeof(*transfer_struct_c.ts_fourier_saved_),transfer_struct_c.ts_fourier_saved__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_geometries_",sizeof(*transfer_struct_c.ts_n_geometries_),transfer_struct_c.ts_n_geometries__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_intensity_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_intensity_,
      blitz::shape(transfer_struct_c.ts_intensity__f_shapes[0],
                   transfer_struct_c.ts_intensity__f_shapes[1],
                   transfer_struct_c.ts_intensity__f_shapes[2],
                   transfer_struct_c.ts_intensity__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_mean_intensity_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_mean_intensity_,
      blitz::shape(transfer_struct_c.ts_mean_intensity__f_shapes[0],
                   transfer_struct_c.ts_mean_intensity__f_shapes[1],
                   transfer_struct_c.ts_mean_intensity__f_shapes[2],
                   transfer_struct_c.ts_mean_intensity__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_flux_integral_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_flux_integral_,
      blitz::shape(transfer_struct_c.ts_flux_integral__f_shapes[0],
                   transfer_struct_c.ts_flux_integral__f_shapes[1],
                   transfer_struct_c.ts_flux_integral__f_shapes[2],
                   transfer_struct_c.ts_flux_integral__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_dnflux_direct_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_dnflux_direct_,
      blitz::shape(transfer_struct_c.ts_dnflux_direct__f_shapes[0],
                   transfer_struct_c.ts_dnflux_direct__f_shapes[1],
                   transfer_struct_c.ts_dnflux_direct__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_dnmean_direct_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_dnmean_direct_,
      blitz::shape(transfer_struct_c.ts_dnmean_direct__f_shapes[0],
                   transfer_struct_c.ts_dnmean_direct__f_shapes[1],
                   transfer_struct_c.ts_dnmean_direct__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_fourier_saved_.reference(blitz::Array<int, 2>(transfer_struct_c.ts_fourier_saved_,
      blitz::shape(transfer_struct_c.ts_fourier_saved__f_shapes[0],
                   transfer_struct_c.ts_fourier_saved__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_main_outputs transfer_struct_c;

  blitz::Array<double, 4> ts_intensity_;
  blitz::Array<double, 4> ts_mean_intensity_;
  blitz::Array<double, 4> ts_flux_integral_;
  blitz::Array<double, 3> ts_dnflux_direct_;
  blitz::Array<double, 3> ts_dnmean_direct_;
  blitz::Array<int, 2> ts_fourier_saved_;
  
};

// Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
extern "C" {
  void lidort_exception_handling_c_alloc_init(struct lidort_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_exception_handling_c_init_only(struct lidort_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_exception_handling_c_destroy(void **fortran_type_c);
  void exception_handling_ts_checkmessages_get(void **fortran_type_c, const int* ts_checkmessages_in_shape_1, const int* ts_checkmessages_in_len, const char* ts_checkmessages_in);
  void exception_handling_ts_actions_get(void **fortran_type_c, const int* ts_actions_in_shape_1, const int* ts_actions_in_len, const char* ts_actions_in);
  void exception_handling_ts_message_get(void **fortran_type_c, const int* ts_message_in_len, const char* ts_message_in);
  void exception_handling_ts_trace_1_get(void **fortran_type_c, const int* ts_trace_1_in_len, const char* ts_trace_1_in);
  void exception_handling_ts_trace_2_get(void **fortran_type_c, const int* ts_trace_2_in_len, const char* ts_trace_2_in);
  void exception_handling_ts_trace_3_get(void **fortran_type_c, const int* ts_trace_3_in_len, const char* ts_trace_3_in);
  
}

struct lidort_exception_handling {
  int* ts_status_inputcheck_;
  int ts_status_inputcheck__f_byte_size;

  int* ts_ncheckmessages_;
  int ts_ncheckmessages__f_byte_size;

  
  int ts_checkmessages__f_shapes[1];
  int ts_checkmessages__f_len;

  
  int ts_actions__f_shapes[1];
  int ts_actions__f_len;

  int* ts_status_calculation_;
  int ts_status_calculation__f_byte_size;

  
  int ts_message__f_len;

  
  int ts_trace_1__f_len;

  
  int ts_trace_2__f_len;

  
  int ts_trace_3__f_len;

  
};

// Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
class Lidort_Exception_Handling : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Exception_Handling() : Lidort_Structure() {
    lidort_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Exception_Handling(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Exception_Handling() {
    if (owns_pointer)
      lidort_exception_handling_c_destroy(&fortran_type_c);
  }

  const int& ts_status_inputcheck() const {
    return *transfer_struct_c.ts_status_inputcheck_;
  }
  void ts_status_inputcheck(const int& ts_status_inputcheck_in) {
    *transfer_struct_c.ts_status_inputcheck_ = ts_status_inputcheck_in;
  }
  
  const int& ts_ncheckmessages() const {
    return *transfer_struct_c.ts_ncheckmessages_;
  }
  void ts_ncheckmessages(const int& ts_ncheckmessages_in) {
    *transfer_struct_c.ts_ncheckmessages_ = ts_ncheckmessages_in;
  }
  
  const std::vector< std::string > ts_checkmessages() const {
    std::vector< std::string > ts_checkmessages_ret;
    blitz::Array<char, 2> ts_checkmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_checkmessages__f_shapes[0], transfer_struct_c.ts_checkmessages__f_len+1, blitz::ColumnMajorArray<2>());
    exception_handling_ts_checkmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_checkmessages__f_shapes[0], &transfer_struct_c.ts_checkmessages__f_len, ts_checkmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_checkmessages_lcl.extent(0); dim_0_idx++)
      ts_checkmessages_ret.push_back( std::string(std::string(ts_checkmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_checkmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_checkmessages_ret;
  }
  
  const std::vector< std::string > ts_actions() const {
    std::vector< std::string > ts_actions_ret;
    blitz::Array<char, 2> ts_actions_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_actions__f_shapes[0], transfer_struct_c.ts_actions__f_len+1, blitz::ColumnMajorArray<2>());
    exception_handling_ts_actions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_actions__f_shapes[0], &transfer_struct_c.ts_actions__f_len, ts_actions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_actions_lcl.extent(0); dim_0_idx++)
      ts_actions_ret.push_back( std::string(std::string(ts_actions_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_actions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_actions_ret;
  }
  
  const int& ts_status_calculation() const {
    return *transfer_struct_c.ts_status_calculation_;
  }
  void ts_status_calculation(const int& ts_status_calculation_in) {
    *transfer_struct_c.ts_status_calculation_ = ts_status_calculation_in;
  }
  
  const std::string ts_message() const {
    std::string ts_message_ret;
    blitz::Array<char, 1> ts_message_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_message__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_message_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_message__f_len, ts_message_lcl.dataFirst());
    ts_message_ret = ( std::string(std::string(ts_message_lcl(blitz::Range::all()).begin(), ts_message_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_message_ret;
  }
  
  const std::string ts_trace_1() const {
    std::string ts_trace_1_ret;
    blitz::Array<char, 1> ts_trace_1_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_1__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_1_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_1__f_len, ts_trace_1_lcl.dataFirst());
    ts_trace_1_ret = ( std::string(std::string(ts_trace_1_lcl(blitz::Range::all()).begin(), ts_trace_1_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_1_ret;
  }
  
  const std::string ts_trace_2() const {
    std::string ts_trace_2_ret;
    blitz::Array<char, 1> ts_trace_2_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_2__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_2_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_2__f_len, ts_trace_2_lcl.dataFirst());
    ts_trace_2_ret = ( std::string(std::string(ts_trace_2_lcl(blitz::Range::all()).begin(), ts_trace_2_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_2_ret;
  }
  
  const std::string ts_trace_3() const {
    std::string ts_trace_3_ret;
    blitz::Array<char, 1> ts_trace_3_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_3__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_3_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_3__f_len, ts_trace_3_lcl.dataFirst());
    ts_trace_3_ret = ( std::string(std::string(ts_trace_3_lcl(blitz::Range::all()).begin(), ts_trace_3_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_3_ret;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Exception_Handling:" << std::endl
      << " ts_status_inputcheck: " << ts_status_inputcheck()  << std::endl
      << "    ts_ncheckmessages: " << ts_ncheckmessages()  << std::endl
      << "     ts_checkmessages: " << std::endl;
    std::vector< std::string > ts_checkmessages_lcl = ts_checkmessages();
    for(unsigned int idx = 0; idx < ts_checkmessages_lcl.size(); idx++)
      if ( ts_checkmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_checkmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "           ts_actions: " << std::endl;
    std::vector< std::string > ts_actions_lcl = ts_actions();
    for(unsigned int idx = 0; idx < ts_actions_lcl.size(); idx++)
      if ( ts_actions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_actions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "ts_status_calculation: " << ts_status_calculation()  << std::endl
      << "           ts_message: " << "\"" << ts_message() << "\"" << std::endl
      << "           ts_trace_1: " << "\"" << ts_trace_1() << "\"" << std::endl
      << "           ts_trace_2: " << "\"" << ts_trace_2() << "\"" << std::endl
      << "           ts_trace_3: " << "\"" << ts_trace_3() << "\"" << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_status_inputcheck_",sizeof(*transfer_struct_c.ts_status_inputcheck_),transfer_struct_c.ts_status_inputcheck__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ncheckmessages_",sizeof(*transfer_struct_c.ts_ncheckmessages_),transfer_struct_c.ts_ncheckmessages__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_status_calculation_",sizeof(*transfer_struct_c.ts_status_calculation_),transfer_struct_c.ts_status_calculation__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_exception_handling transfer_struct_c;

  
};

// Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
extern "C" {
  void lidort_input_exception_handling_c_alloc_init(struct lidort_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_input_exception_handling_c_init_only(struct lidort_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_input_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_input_exception_handling_c_destroy(void **fortran_type_c);
  void input_exception_handling_ts_inputmessages_get(void **fortran_type_c, const int* ts_inputmessages_in_shape_1, const int* ts_inputmessages_in_len, const char* ts_inputmessages_in);
  void input_exception_handling_ts_inputactions_get(void **fortran_type_c, const int* ts_inputactions_in_shape_1, const int* ts_inputactions_in_len, const char* ts_inputactions_in);
  
}

struct lidort_input_exception_handling {
  int* ts_status_inputread_;
  int ts_status_inputread__f_byte_size;

  int* ts_ninputmessages_;
  int ts_ninputmessages__f_byte_size;

  
  int ts_inputmessages__f_shapes[1];
  int ts_inputmessages__f_len;

  
  int ts_inputactions__f_shapes[1];
  int ts_inputactions__f_len;

  
};

// Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
class Lidort_Input_Exception_Handling : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Input_Exception_Handling() : Lidort_Structure() {
    lidort_input_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Input_Exception_Handling(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_input_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Input_Exception_Handling() {
    if (owns_pointer)
      lidort_input_exception_handling_c_destroy(&fortran_type_c);
  }

  const int& ts_status_inputread() const {
    return *transfer_struct_c.ts_status_inputread_;
  }
  void ts_status_inputread(const int& ts_status_inputread_in) {
    *transfer_struct_c.ts_status_inputread_ = ts_status_inputread_in;
  }
  
  const int& ts_ninputmessages() const {
    return *transfer_struct_c.ts_ninputmessages_;
  }
  void ts_ninputmessages(const int& ts_ninputmessages_in) {
    *transfer_struct_c.ts_ninputmessages_ = ts_ninputmessages_in;
  }
  
  const std::vector< std::string > ts_inputmessages() const {
    std::vector< std::string > ts_inputmessages_ret;
    blitz::Array<char, 2> ts_inputmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_inputmessages__f_shapes[0], transfer_struct_c.ts_inputmessages__f_len+1, blitz::ColumnMajorArray<2>());
    input_exception_handling_ts_inputmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_inputmessages__f_shapes[0], &transfer_struct_c.ts_inputmessages__f_len, ts_inputmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_inputmessages_lcl.extent(0); dim_0_idx++)
      ts_inputmessages_ret.push_back( std::string(std::string(ts_inputmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_inputmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_inputmessages_ret;
  }
  
  const std::vector< std::string > ts_inputactions() const {
    std::vector< std::string > ts_inputactions_ret;
    blitz::Array<char, 2> ts_inputactions_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_inputactions__f_shapes[0], transfer_struct_c.ts_inputactions__f_len+1, blitz::ColumnMajorArray<2>());
    input_exception_handling_ts_inputactions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_inputactions__f_shapes[0], &transfer_struct_c.ts_inputactions__f_len, ts_inputactions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_inputactions_lcl.extent(0); dim_0_idx++)
      ts_inputactions_ret.push_back( std::string(std::string(ts_inputactions_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_inputactions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_inputactions_ret;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Input_Exception_Handling:" << std::endl
      << "ts_status_inputread: " << ts_status_inputread()  << std::endl
      << "  ts_ninputmessages: " << ts_ninputmessages()  << std::endl
      << "   ts_inputmessages: " << std::endl;
    std::vector< std::string > ts_inputmessages_lcl = ts_inputmessages();
    for(unsigned int idx = 0; idx < ts_inputmessages_lcl.size(); idx++)
      if ( ts_inputmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_inputmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "    ts_inputactions: " << std::endl;
    std::vector< std::string > ts_inputactions_lcl = ts_inputactions();
    for(unsigned int idx = 0; idx < ts_inputactions_lcl.size(); idx++)
      if ( ts_inputactions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_inputactions_lcl[idx] << "\"" << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_status_inputread_",sizeof(*transfer_struct_c.ts_status_inputread_),transfer_struct_c.ts_status_inputread__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ninputmessages_",sizeof(*transfer_struct_c.ts_ninputmessages_),transfer_struct_c.ts_ninputmessages__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_input_exception_handling transfer_struct_c;

  
};

// Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
extern "C" {
  void lidort_outputs_c_alloc_init(struct lidort_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_outputs_c_init_only(struct lidort_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_outputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_outputs {
  void* main_;
  int main__f_byte_size;

  void* status_;
  int status__f_byte_size;

  
};

// Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
class Lidort_Outputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Outputs() : Lidort_Structure() {
    lidort_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Outputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Outputs() {
    if (owns_pointer)
      lidort_outputs_c_destroy(&fortran_type_c);
  }

  Lidort_Main_Outputs& main() {
    return *main_;
  }
  const Lidort_Main_Outputs& main() const {
    return *main_;
  }
  void main(Lidort_Main_Outputs& main_in) {
    void* src_ptr = main_in.fortran_type_ptr();
    void* dst_ptr = main_->fortran_type_ptr();
    lidort_main_outputs_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Exception_Handling& status() {
    return *status_;
  }
  const Lidort_Exception_Handling& status() const {
    return *status_;
  }
  void status(Lidort_Exception_Handling& status_in) {
    void* src_ptr = status_in.fortran_type_ptr();
    void* dst_ptr = status_->fortran_type_ptr();
    lidort_exception_handling_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Outputs:" << std::endl
      << "  main: " << main()  << std::endl
      << "status: " << status()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    main_.reset( new Lidort_Main_Outputs(transfer_struct_c.main_) );
    status_.reset( new Lidort_Exception_Handling(transfer_struct_c.status_) );
    
  }

  struct lidort_outputs transfer_struct_c;

  boost::shared_ptr<Lidort_Main_Outputs> main_;
  boost::shared_ptr<Lidort_Exception_Handling> status_;
  
};

// Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
extern "C" {
  void lidort_sup_brdf_c_alloc_init(struct lidort_sup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_brdf_c_init_only(struct lidort_sup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_brdf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_brdf_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_brdf {
  double* ts_exactdb_brdfunc_;
  int ts_exactdb_brdfunc__f_shapes[3];
  int ts_exactdb_brdfunc__f_byte_size;

  double* ts_brdf_f_0_;
  int ts_brdf_f_0__f_shapes[3];
  int ts_brdf_f_0__f_byte_size;

  double* ts_brdf_f_;
  int ts_brdf_f__f_shapes[3];
  int ts_brdf_f__f_byte_size;

  double* ts_user_brdf_f_0_;
  int ts_user_brdf_f_0__f_shapes[3];
  int ts_user_brdf_f_0__f_byte_size;

  double* ts_user_brdf_f_;
  int ts_user_brdf_f__f_shapes[3];
  int ts_user_brdf_f__f_byte_size;

  double* ts_emissivity_;
  int ts_emissivity__f_shapes[2];
  int ts_emissivity__f_byte_size;

  double* ts_user_emissivity_;
  int ts_user_emissivity__f_shapes[2];
  int ts_user_emissivity__f_byte_size;

  
};

// Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
class Lidort_Sup_Brdf : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Sup_Brdf() : Lidort_Structure() {
    lidort_sup_brdf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Brdf(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_sup_brdf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Brdf() {
    if (owns_pointer)
      lidort_sup_brdf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_exactdb_brdfunc() const {
    return ts_exactdb_brdfunc_;
  }
  void ts_exactdb_brdfunc(const blitz::Array<double, 3>& ts_exactdb_brdfunc_in) {
    ts_exactdb_brdfunc_ = ts_exactdb_brdfunc_in;
  }
  
  const blitz::Array<double, 3>& ts_brdf_f_0() const {
    return ts_brdf_f_0_;
  }
  void ts_brdf_f_0(const blitz::Array<double, 3>& ts_brdf_f_0_in) {
    ts_brdf_f_0_ = ts_brdf_f_0_in;
  }
  
  const blitz::Array<double, 3>& ts_brdf_f() const {
    return ts_brdf_f_;
  }
  void ts_brdf_f(const blitz::Array<double, 3>& ts_brdf_f_in) {
    ts_brdf_f_ = ts_brdf_f_in;
  }
  
  const blitz::Array<double, 3>& ts_user_brdf_f_0() const {
    return ts_user_brdf_f_0_;
  }
  void ts_user_brdf_f_0(const blitz::Array<double, 3>& ts_user_brdf_f_0_in) {
    ts_user_brdf_f_0_ = ts_user_brdf_f_0_in;
  }
  
  const blitz::Array<double, 3>& ts_user_brdf_f() const {
    return ts_user_brdf_f_;
  }
  void ts_user_brdf_f(const blitz::Array<double, 3>& ts_user_brdf_f_in) {
    ts_user_brdf_f_ = ts_user_brdf_f_in;
  }
  
  const blitz::Array<double, 2>& ts_emissivity() const {
    return ts_emissivity_;
  }
  void ts_emissivity(const blitz::Array<double, 2>& ts_emissivity_in) {
    ts_emissivity_ = ts_emissivity_in;
  }
  
  const blitz::Array<double, 2>& ts_user_emissivity() const {
    return ts_user_emissivity_;
  }
  void ts_user_emissivity(const blitz::Array<double, 2>& ts_user_emissivity_in) {
    ts_user_emissivity_ = ts_user_emissivity_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Brdf:" << std::endl
      << "ts_exactdb_brdfunc: " << std::endl << ts_exactdb_brdfunc()  << std::endl
      << "       ts_brdf_f_0: " << std::endl << ts_brdf_f_0()  << std::endl
      << "         ts_brdf_f: " << std::endl << ts_brdf_f()  << std::endl
      << "  ts_user_brdf_f_0: " << std::endl << ts_user_brdf_f_0()  << std::endl
      << "    ts_user_brdf_f: " << std::endl << ts_user_brdf_f()  << std::endl
      << "     ts_emissivity: " << std::endl << ts_emissivity()  << std::endl
      << "ts_user_emissivity: " << std::endl << ts_user_emissivity()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_exactdb_brdfunc_",sizeof(*transfer_struct_c.ts_exactdb_brdfunc_),transfer_struct_c.ts_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_brdf_f_0_",sizeof(*transfer_struct_c.ts_brdf_f_0_),transfer_struct_c.ts_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_brdf_f_",sizeof(*transfer_struct_c.ts_brdf_f_),transfer_struct_c.ts_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_brdf_f_0_",sizeof(*transfer_struct_c.ts_user_brdf_f_0_),transfer_struct_c.ts_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_brdf_f_",sizeof(*transfer_struct_c.ts_user_brdf_f_),transfer_struct_c.ts_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_emissivity_",sizeof(*transfer_struct_c.ts_emissivity_),transfer_struct_c.ts_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_emissivity_",sizeof(*transfer_struct_c.ts_user_emissivity_),transfer_struct_c.ts_user_emissivity__f_byte_size);
    
  }
  
  void copy_from_sup(Brdf_Sup_Outputs& supp_obj) { 
    // This is only safe for LIDORT 3.6 because the BRDF and LIDORT
    // sup structures are the exact same structure, but with different
    // names. This MUST be reevaluated on future LIDORT versions
    void* sup_ptr = supp_obj.fortran_type_ptr();
    lidort_sup_brdf_c_copy(&sup_ptr, &fortran_type_c);
  }

private:
  void link_blitz_arrays() {
    ts_exactdb_brdfunc_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.ts_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.ts_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.ts_exactdb_brdfunc__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_brdf_f_,
      blitz::shape(transfer_struct_c.ts_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_user_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_brdf_f_,
      blitz::shape(transfer_struct_c.ts_user_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_user_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_user_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_emissivity_,
      blitz::shape(transfer_struct_c.ts_emissivity__f_shapes[0],
                   transfer_struct_c.ts_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_user_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_user_emissivity_,
      blitz::shape(transfer_struct_c.ts_user_emissivity__f_shapes[0],
                   transfer_struct_c.ts_user_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_brdf transfer_struct_c;

  blitz::Array<double, 3> ts_exactdb_brdfunc_;
  blitz::Array<double, 3> ts_brdf_f_0_;
  blitz::Array<double, 3> ts_brdf_f_;
  blitz::Array<double, 3> ts_user_brdf_f_0_;
  blitz::Array<double, 3> ts_user_brdf_f_;
  blitz::Array<double, 2> ts_emissivity_;
  blitz::Array<double, 2> ts_user_emissivity_;
  
};

// Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
extern "C" {
  void lidort_sup_sleave_c_alloc_init(struct lidort_sup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_sleave_c_init_only(struct lidort_sup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_sleave_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_sleave_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_sleave {
  double* ts_slterm_isotropic_;
  int ts_slterm_isotropic__f_shapes[1];
  int ts_slterm_isotropic__f_byte_size;

  double* ts_slterm_userangles_;
  int ts_slterm_userangles__f_shapes[3];
  int ts_slterm_userangles__f_byte_size;

  double* ts_slterm_f_0_;
  int ts_slterm_f_0__f_shapes[3];
  int ts_slterm_f_0__f_byte_size;

  double* ts_user_slterm_f_0_;
  int ts_user_slterm_f_0__f_shapes[3];
  int ts_user_slterm_f_0__f_byte_size;

  
};

// Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
class Lidort_Sup_Sleave : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Sup_Sleave() : Lidort_Structure() {
    lidort_sup_sleave_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Sleave(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_sup_sleave_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Sleave() {
    if (owns_pointer)
      lidort_sup_sleave_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& ts_slterm_isotropic() const {
    return ts_slterm_isotropic_;
  }
  void ts_slterm_isotropic(const blitz::Array<double, 1>& ts_slterm_isotropic_in) {
    ts_slterm_isotropic_ = ts_slterm_isotropic_in;
  }
  
  const blitz::Array<double, 3>& ts_slterm_userangles() const {
    return ts_slterm_userangles_;
  }
  void ts_slterm_userangles(const blitz::Array<double, 3>& ts_slterm_userangles_in) {
    ts_slterm_userangles_ = ts_slterm_userangles_in;
  }
  
  const blitz::Array<double, 3>& ts_slterm_f_0() const {
    return ts_slterm_f_0_;
  }
  void ts_slterm_f_0(const blitz::Array<double, 3>& ts_slterm_f_0_in) {
    ts_slterm_f_0_ = ts_slterm_f_0_in;
  }
  
  const blitz::Array<double, 3>& ts_user_slterm_f_0() const {
    return ts_user_slterm_f_0_;
  }
  void ts_user_slterm_f_0(const blitz::Array<double, 3>& ts_user_slterm_f_0_in) {
    ts_user_slterm_f_0_ = ts_user_slterm_f_0_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Sleave:" << std::endl
      << " ts_slterm_isotropic: " << std::endl << ts_slterm_isotropic()  << std::endl
      << "ts_slterm_userangles: " << std::endl << ts_slterm_userangles()  << std::endl
      << "       ts_slterm_f_0: " << std::endl << ts_slterm_f_0()  << std::endl
      << "  ts_user_slterm_f_0: " << std::endl << ts_user_slterm_f_0()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_slterm_isotropic_",sizeof(*transfer_struct_c.ts_slterm_isotropic_),transfer_struct_c.ts_slterm_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_slterm_userangles_",sizeof(*transfer_struct_c.ts_slterm_userangles_),transfer_struct_c.ts_slterm_userangles__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_slterm_f_0_",sizeof(*transfer_struct_c.ts_slterm_f_0_),transfer_struct_c.ts_slterm_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_slterm_f_0_",sizeof(*transfer_struct_c.ts_user_slterm_f_0_),transfer_struct_c.ts_user_slterm_f_0__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_slterm_isotropic_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_slterm_isotropic_,
      blitz::shape(transfer_struct_c.ts_slterm_isotropic__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_slterm_userangles_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_slterm_userangles_,
      blitz::shape(transfer_struct_c.ts_slterm_userangles__f_shapes[0],
                   transfer_struct_c.ts_slterm_userangles__f_shapes[1],
                   transfer_struct_c.ts_slterm_userangles__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_slterm_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_slterm_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_slterm_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_user_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_user_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_user_slterm_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_sleave transfer_struct_c;

  blitz::Array<double, 1> ts_slterm_isotropic_;
  blitz::Array<double, 3> ts_slterm_userangles_;
  blitz::Array<double, 3> ts_slterm_f_0_;
  blitz::Array<double, 3> ts_user_slterm_f_0_;
  
};

// Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
extern "C" {
  void lidort_sup_ss_c_alloc_init(struct lidort_sup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_ss_c_init_only(struct lidort_sup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_ss_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_ss_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_ss {
  double* ts_intensity_ss_;
  int ts_intensity_ss__f_shapes[3];
  int ts_intensity_ss__f_byte_size;

  double* ts_intensity_db_;
  int ts_intensity_db__f_shapes[2];
  int ts_intensity_db__f_byte_size;

  
};

// Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
class Lidort_Sup_Ss : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Sup_Ss() : Lidort_Structure() {
    lidort_sup_ss_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Ss(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_sup_ss_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Ss() {
    if (owns_pointer)
      lidort_sup_ss_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_intensity_ss() const {
    return ts_intensity_ss_;
  }
  void ts_intensity_ss(const blitz::Array<double, 3>& ts_intensity_ss_in) {
    ts_intensity_ss_ = ts_intensity_ss_in;
  }
  
  const blitz::Array<double, 2>& ts_intensity_db() const {
    return ts_intensity_db_;
  }
  void ts_intensity_db(const blitz::Array<double, 2>& ts_intensity_db_in) {
    ts_intensity_db_ = ts_intensity_db_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Ss:" << std::endl
      << "ts_intensity_ss: " << std::endl << ts_intensity_ss()  << std::endl
      << "ts_intensity_db: " << std::endl << ts_intensity_db()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_intensity_ss_",sizeof(*transfer_struct_c.ts_intensity_ss_),transfer_struct_c.ts_intensity_ss__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_intensity_db_",sizeof(*transfer_struct_c.ts_intensity_db_),transfer_struct_c.ts_intensity_db__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_intensity_ss_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_intensity_ss_,
      blitz::shape(transfer_struct_c.ts_intensity_ss__f_shapes[0],
                   transfer_struct_c.ts_intensity_ss__f_shapes[1],
                   transfer_struct_c.ts_intensity_ss__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_intensity_db_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_intensity_db_,
      blitz::shape(transfer_struct_c.ts_intensity_db__f_shapes[0],
                   transfer_struct_c.ts_intensity_db__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_ss transfer_struct_c;

  blitz::Array<double, 3> ts_intensity_ss_;
  blitz::Array<double, 2> ts_intensity_db_;
  
};

// Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
extern "C" {
  void lidort_sup_inout_c_alloc_init(struct lidort_sup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_inout_c_init_only(struct lidort_sup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_inout_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_inout_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_inout {
  void* brdf_;
  int brdf__f_byte_size;

  void* ss_;
  int ss__f_byte_size;

  void* sleave_;
  int sleave__f_byte_size;

  
};

// Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
class Lidort_Sup_Inout : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Sup_Inout() : Lidort_Structure() {
    lidort_sup_inout_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Inout(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_sup_inout_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Inout() {
    if (owns_pointer)
      lidort_sup_inout_c_destroy(&fortran_type_c);
  }

  Lidort_Sup_Brdf& brdf() {
    return *brdf_;
  }
  const Lidort_Sup_Brdf& brdf() const {
    return *brdf_;
  }
  void brdf(Lidort_Sup_Brdf& brdf_in) {
    void* src_ptr = brdf_in.fortran_type_ptr();
    void* dst_ptr = brdf_->fortran_type_ptr();
    lidort_sup_brdf_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Sup_Ss& ss() {
    return *ss_;
  }
  const Lidort_Sup_Ss& ss() const {
    return *ss_;
  }
  void ss(Lidort_Sup_Ss& ss_in) {
    void* src_ptr = ss_in.fortran_type_ptr();
    void* dst_ptr = ss_->fortran_type_ptr();
    lidort_sup_ss_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Sup_Sleave& sleave() {
    return *sleave_;
  }
  const Lidort_Sup_Sleave& sleave() const {
    return *sleave_;
  }
  void sleave(Lidort_Sup_Sleave& sleave_in) {
    void* src_ptr = sleave_in.fortran_type_ptr();
    void* dst_ptr = sleave_->fortran_type_ptr();
    lidort_sup_sleave_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Inout:" << std::endl
      << "  brdf: " << brdf()  << std::endl
      << "    ss: " << ss()  << std::endl
      << "sleave: " << sleave()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    brdf_.reset( new Lidort_Sup_Brdf(transfer_struct_c.brdf_) );
    ss_.reset( new Lidort_Sup_Ss(transfer_struct_c.ss_) );
    sleave_.reset( new Lidort_Sup_Sleave(transfer_struct_c.sleave_) );
    
  }

  struct lidort_sup_inout transfer_struct_c;

  boost::shared_ptr<Lidort_Sup_Brdf> brdf_;
  boost::shared_ptr<Lidort_Sup_Ss> ss_;
  boost::shared_ptr<Lidort_Sup_Sleave> sleave_;
  
};

// Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_boolean_c_alloc_init(struct lidort_fixed_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_boolean_c_init_only(struct lidort_fixed_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_boolean_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_boolean_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_boolean {
  int* ts_do_fullrad_mode_;
  int ts_do_fullrad_mode__f_byte_size;

  int* ts_do_sscorr_truncation_;
  int ts_do_sscorr_truncation__f_byte_size;

  int* ts_do_ss_external_;
  int ts_do_ss_external__f_byte_size;

  int* ts_do_ssfull_;
  int ts_do_ssfull__f_byte_size;

  int* ts_do_thermal_emission_;
  int ts_do_thermal_emission__f_byte_size;

  int* ts_do_surface_emission_;
  int ts_do_surface_emission__f_byte_size;

  int* ts_do_plane_parallel_;
  int ts_do_plane_parallel__f_byte_size;

  int* ts_do_brdf_surface_;
  int ts_do_brdf_surface__f_byte_size;

  int* ts_do_upwelling_;
  int ts_do_upwelling__f_byte_size;

  int* ts_do_dnwelling_;
  int ts_do_dnwelling__f_byte_size;

  int* ts_do_surface_leaving_;
  int ts_do_surface_leaving__f_byte_size;

  int* ts_do_sl_isotropic_;
  int ts_do_sl_isotropic__f_byte_size;

  
};

// Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Boolean : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Boolean() : Lidort_Structure() {
    lidort_fixed_boolean_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Boolean(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_boolean_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Boolean() {
    if (owns_pointer)
      lidort_fixed_boolean_c_destroy(&fortran_type_c);
  }

  const bool ts_do_fullrad_mode() const {
    return *transfer_struct_c.ts_do_fullrad_mode_ != 0;
  }
  void ts_do_fullrad_mode(const bool& ts_do_fullrad_mode_in) {
    *transfer_struct_c.ts_do_fullrad_mode_ = ts_do_fullrad_mode_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_sscorr_truncation() const {
    return *transfer_struct_c.ts_do_sscorr_truncation_ != 0;
  }
  void ts_do_sscorr_truncation(const bool& ts_do_sscorr_truncation_in) {
    *transfer_struct_c.ts_do_sscorr_truncation_ = ts_do_sscorr_truncation_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_ss_external() const {
    return *transfer_struct_c.ts_do_ss_external_ != 0;
  }
  void ts_do_ss_external(const bool& ts_do_ss_external_in) {
    *transfer_struct_c.ts_do_ss_external_ = ts_do_ss_external_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_ssfull() const {
    return *transfer_struct_c.ts_do_ssfull_ != 0;
  }
  void ts_do_ssfull(const bool& ts_do_ssfull_in) {
    *transfer_struct_c.ts_do_ssfull_ = ts_do_ssfull_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_thermal_emission() const {
    return *transfer_struct_c.ts_do_thermal_emission_ != 0;
  }
  void ts_do_thermal_emission(const bool& ts_do_thermal_emission_in) {
    *transfer_struct_c.ts_do_thermal_emission_ = ts_do_thermal_emission_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_surface_emission() const {
    return *transfer_struct_c.ts_do_surface_emission_ != 0;
  }
  void ts_do_surface_emission(const bool& ts_do_surface_emission_in) {
    *transfer_struct_c.ts_do_surface_emission_ = ts_do_surface_emission_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_plane_parallel() const {
    return *transfer_struct_c.ts_do_plane_parallel_ != 0;
  }
  void ts_do_plane_parallel(const bool& ts_do_plane_parallel_in) {
    *transfer_struct_c.ts_do_plane_parallel_ = ts_do_plane_parallel_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_brdf_surface() const {
    return *transfer_struct_c.ts_do_brdf_surface_ != 0;
  }
  void ts_do_brdf_surface(const bool& ts_do_brdf_surface_in) {
    *transfer_struct_c.ts_do_brdf_surface_ = ts_do_brdf_surface_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_upwelling() const {
    return *transfer_struct_c.ts_do_upwelling_ != 0;
  }
  void ts_do_upwelling(const bool& ts_do_upwelling_in) {
    *transfer_struct_c.ts_do_upwelling_ = ts_do_upwelling_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_dnwelling() const {
    return *transfer_struct_c.ts_do_dnwelling_ != 0;
  }
  void ts_do_dnwelling(const bool& ts_do_dnwelling_in) {
    *transfer_struct_c.ts_do_dnwelling_ = ts_do_dnwelling_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_surface_leaving() const {
    return *transfer_struct_c.ts_do_surface_leaving_ != 0;
  }
  void ts_do_surface_leaving(const bool& ts_do_surface_leaving_in) {
    *transfer_struct_c.ts_do_surface_leaving_ = ts_do_surface_leaving_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_sl_isotropic() const {
    return *transfer_struct_c.ts_do_sl_isotropic_ != 0;
  }
  void ts_do_sl_isotropic(const bool& ts_do_sl_isotropic_in) {
    *transfer_struct_c.ts_do_sl_isotropic_ = ts_do_sl_isotropic_in ? FORTRAN_TRUE_INT : 0;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Boolean:" << std::endl
      << "     ts_do_fullrad_mode: " << ts_do_fullrad_mode()  << std::endl
      << "ts_do_sscorr_truncation: " << ts_do_sscorr_truncation()  << std::endl
      << "      ts_do_ss_external: " << ts_do_ss_external()  << std::endl
      << "           ts_do_ssfull: " << ts_do_ssfull()  << std::endl
      << " ts_do_thermal_emission: " << ts_do_thermal_emission()  << std::endl
      << " ts_do_surface_emission: " << ts_do_surface_emission()  << std::endl
      << "   ts_do_plane_parallel: " << ts_do_plane_parallel()  << std::endl
      << "     ts_do_brdf_surface: " << ts_do_brdf_surface()  << std::endl
      << "        ts_do_upwelling: " << ts_do_upwelling()  << std::endl
      << "        ts_do_dnwelling: " << ts_do_dnwelling()  << std::endl
      << "  ts_do_surface_leaving: " << ts_do_surface_leaving()  << std::endl
      << "     ts_do_sl_isotropic: " << ts_do_sl_isotropic()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_fullrad_mode_",sizeof(*transfer_struct_c.ts_do_fullrad_mode_),transfer_struct_c.ts_do_fullrad_mode__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sscorr_truncation_",sizeof(*transfer_struct_c.ts_do_sscorr_truncation_),transfer_struct_c.ts_do_sscorr_truncation__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_ss_external_",sizeof(*transfer_struct_c.ts_do_ss_external_),transfer_struct_c.ts_do_ss_external__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_ssfull_",sizeof(*transfer_struct_c.ts_do_ssfull_),transfer_struct_c.ts_do_ssfull__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_thermal_emission_",sizeof(*transfer_struct_c.ts_do_thermal_emission_),transfer_struct_c.ts_do_thermal_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_emission_",sizeof(*transfer_struct_c.ts_do_surface_emission_),transfer_struct_c.ts_do_surface_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_plane_parallel_",sizeof(*transfer_struct_c.ts_do_plane_parallel_),transfer_struct_c.ts_do_plane_parallel__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_brdf_surface_",sizeof(*transfer_struct_c.ts_do_brdf_surface_),transfer_struct_c.ts_do_brdf_surface__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_upwelling_",sizeof(*transfer_struct_c.ts_do_upwelling_),transfer_struct_c.ts_do_upwelling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_dnwelling_",sizeof(*transfer_struct_c.ts_do_dnwelling_),transfer_struct_c.ts_do_dnwelling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_leaving_",sizeof(*transfer_struct_c.ts_do_surface_leaving_),transfer_struct_c.ts_do_surface_leaving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sl_isotropic_",sizeof(*transfer_struct_c.ts_do_sl_isotropic_),transfer_struct_c.ts_do_sl_isotropic__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_boolean transfer_struct_c;

  
};

// Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_control_c_alloc_init(struct lidort_fixed_control *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_control_c_init_only(struct lidort_fixed_control *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_control_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_control_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_control {
  int* ts_nstreams_;
  int ts_nstreams__f_byte_size;

  int* ts_nlayers_;
  int ts_nlayers__f_byte_size;

  int* ts_nfinelayers_;
  int ts_nfinelayers__f_byte_size;

  int* ts_n_thermal_coeffs_;
  int ts_n_thermal_coeffs__f_byte_size;

  double* ts_lidort_accuracy_;
  int ts_lidort_accuracy__f_byte_size;

  
};

// Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Control : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Control() : Lidort_Structure() {
    lidort_fixed_control_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Control(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_control_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Control() {
    if (owns_pointer)
      lidort_fixed_control_c_destroy(&fortran_type_c);
  }

  const int& ts_nstreams() const {
    return *transfer_struct_c.ts_nstreams_;
  }
  void ts_nstreams(const int& ts_nstreams_in) {
    *transfer_struct_c.ts_nstreams_ = ts_nstreams_in;
  }
  
  const int& ts_nlayers() const {
    return *transfer_struct_c.ts_nlayers_;
  }
  void ts_nlayers(const int& ts_nlayers_in) {
    *transfer_struct_c.ts_nlayers_ = ts_nlayers_in;
  }
  
  const int& ts_nfinelayers() const {
    return *transfer_struct_c.ts_nfinelayers_;
  }
  void ts_nfinelayers(const int& ts_nfinelayers_in) {
    *transfer_struct_c.ts_nfinelayers_ = ts_nfinelayers_in;
  }
  
  const int& ts_n_thermal_coeffs() const {
    return *transfer_struct_c.ts_n_thermal_coeffs_;
  }
  void ts_n_thermal_coeffs(const int& ts_n_thermal_coeffs_in) {
    *transfer_struct_c.ts_n_thermal_coeffs_ = ts_n_thermal_coeffs_in;
  }
  
  const double& ts_lidort_accuracy() const {
    return *transfer_struct_c.ts_lidort_accuracy_;
  }
  void ts_lidort_accuracy(const double& ts_lidort_accuracy_in) {
    *transfer_struct_c.ts_lidort_accuracy_ = ts_lidort_accuracy_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Control:" << std::endl
      << "        ts_nstreams: " << ts_nstreams()  << std::endl
      << "         ts_nlayers: " << ts_nlayers()  << std::endl
      << "     ts_nfinelayers: " << ts_nfinelayers()  << std::endl
      << "ts_n_thermal_coeffs: " << ts_n_thermal_coeffs()  << std::endl
      << " ts_lidort_accuracy: " << ts_lidort_accuracy()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_nstreams_",sizeof(*transfer_struct_c.ts_nstreams_),transfer_struct_c.ts_nstreams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nlayers_",sizeof(*transfer_struct_c.ts_nlayers_),transfer_struct_c.ts_nlayers__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nfinelayers_",sizeof(*transfer_struct_c.ts_nfinelayers_),transfer_struct_c.ts_nfinelayers__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_thermal_coeffs_",sizeof(*transfer_struct_c.ts_n_thermal_coeffs_),transfer_struct_c.ts_n_thermal_coeffs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lidort_accuracy_",sizeof(*transfer_struct_c.ts_lidort_accuracy_),transfer_struct_c.ts_lidort_accuracy__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_control transfer_struct_c;

  
};

// Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_sunrays_c_alloc_init(struct lidort_fixed_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_sunrays_c_init_only(struct lidort_fixed_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_sunrays_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_sunrays_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_sunrays {
  double* ts_flux_factor_;
  int ts_flux_factor__f_byte_size;

  
};

// Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Sunrays : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Sunrays() : Lidort_Structure() {
    lidort_fixed_sunrays_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Sunrays(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_sunrays_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Sunrays() {
    if (owns_pointer)
      lidort_fixed_sunrays_c_destroy(&fortran_type_c);
  }

  const double& ts_flux_factor() const {
    return *transfer_struct_c.ts_flux_factor_;
  }
  void ts_flux_factor(const double& ts_flux_factor_in) {
    *transfer_struct_c.ts_flux_factor_ = ts_flux_factor_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Sunrays:" << std::endl
      << "ts_flux_factor: " << ts_flux_factor()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_flux_factor_",sizeof(*transfer_struct_c.ts_flux_factor_),transfer_struct_c.ts_flux_factor__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_sunrays transfer_struct_c;

  
};

// Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_uservalues_c_alloc_init(struct lidort_fixed_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_uservalues_c_init_only(struct lidort_fixed_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_uservalues_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_uservalues_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_uservalues {
  int* ts_n_user_streams_;
  int ts_n_user_streams__f_byte_size;

  int* ts_n_user_levels_;
  int ts_n_user_levels__f_byte_size;

  
};

// Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Uservalues : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Uservalues() : Lidort_Structure() {
    lidort_fixed_uservalues_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Uservalues(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_uservalues_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Uservalues() {
    if (owns_pointer)
      lidort_fixed_uservalues_c_destroy(&fortran_type_c);
  }

  const int& ts_n_user_streams() const {
    return *transfer_struct_c.ts_n_user_streams_;
  }
  void ts_n_user_streams(const int& ts_n_user_streams_in) {
    *transfer_struct_c.ts_n_user_streams_ = ts_n_user_streams_in;
  }
  
  const int& ts_n_user_levels() const {
    return *transfer_struct_c.ts_n_user_levels_;
  }
  void ts_n_user_levels(const int& ts_n_user_levels_in) {
    *transfer_struct_c.ts_n_user_levels_ = ts_n_user_levels_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Uservalues:" << std::endl
      << "ts_n_user_streams: " << ts_n_user_streams()  << std::endl
      << " ts_n_user_levels: " << ts_n_user_levels()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_n_user_streams_",sizeof(*transfer_struct_c.ts_n_user_streams_),transfer_struct_c.ts_n_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_user_levels_",sizeof(*transfer_struct_c.ts_n_user_levels_),transfer_struct_c.ts_n_user_levels__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_uservalues transfer_struct_c;

  
};

// Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_chapman_c_alloc_init(struct lidort_fixed_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_chapman_c_init_only(struct lidort_fixed_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_chapman_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_chapman_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_chapman {
  double* ts_height_grid_;
  int ts_height_grid__f_shapes[1];
  int ts_height_grid__f_byte_size;

  double* ts_pressure_grid_;
  int ts_pressure_grid__f_shapes[1];
  int ts_pressure_grid__f_byte_size;

  double* ts_temperature_grid_;
  int ts_temperature_grid__f_shapes[1];
  int ts_temperature_grid__f_byte_size;

  int* ts_finegrid_;
  int ts_finegrid__f_shapes[1];
  int ts_finegrid__f_byte_size;

  double* ts_rfindex_parameter_;
  int ts_rfindex_parameter__f_byte_size;

  
};

// Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Chapman : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Chapman() : Lidort_Structure() {
    lidort_fixed_chapman_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Chapman(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_chapman_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Chapman() {
    if (owns_pointer)
      lidort_fixed_chapman_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& ts_height_grid() const {
    return ts_height_grid_;
  }
  void ts_height_grid(const blitz::Array<double, 1>& ts_height_grid_in) {
    ts_height_grid_ = ts_height_grid_in;
  }
  
  const blitz::Array<double, 1>& ts_pressure_grid() const {
    return ts_pressure_grid_;
  }
  void ts_pressure_grid(const blitz::Array<double, 1>& ts_pressure_grid_in) {
    ts_pressure_grid_ = ts_pressure_grid_in;
  }
  
  const blitz::Array<double, 1>& ts_temperature_grid() const {
    return ts_temperature_grid_;
  }
  void ts_temperature_grid(const blitz::Array<double, 1>& ts_temperature_grid_in) {
    ts_temperature_grid_ = ts_temperature_grid_in;
  }
  
  const blitz::Array<int, 1>& ts_finegrid() const {
    return ts_finegrid_;
  }
  void ts_finegrid(const blitz::Array<int, 1>& ts_finegrid_in) {
    ts_finegrid_ = ts_finegrid_in;
  }
  
  const double& ts_rfindex_parameter() const {
    return *transfer_struct_c.ts_rfindex_parameter_;
  }
  void ts_rfindex_parameter(const double& ts_rfindex_parameter_in) {
    *transfer_struct_c.ts_rfindex_parameter_ = ts_rfindex_parameter_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Chapman:" << std::endl
      << "      ts_height_grid: " << std::endl << ts_height_grid()  << std::endl
      << "    ts_pressure_grid: " << std::endl << ts_pressure_grid()  << std::endl
      << " ts_temperature_grid: " << std::endl << ts_temperature_grid()  << std::endl
      << "         ts_finegrid: " << std::endl << ts_finegrid()  << std::endl
      << "ts_rfindex_parameter: " << ts_rfindex_parameter()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_height_grid_",sizeof(*transfer_struct_c.ts_height_grid_),transfer_struct_c.ts_height_grid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_pressure_grid_",sizeof(*transfer_struct_c.ts_pressure_grid_),transfer_struct_c.ts_pressure_grid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_temperature_grid_",sizeof(*transfer_struct_c.ts_temperature_grid_),transfer_struct_c.ts_temperature_grid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_finegrid_",sizeof(*transfer_struct_c.ts_finegrid_),transfer_struct_c.ts_finegrid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_rfindex_parameter_",sizeof(*transfer_struct_c.ts_rfindex_parameter_),transfer_struct_c.ts_rfindex_parameter__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_height_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_height_grid_,
      blitz::shape(transfer_struct_c.ts_height_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_pressure_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_pressure_grid_,
      blitz::shape(transfer_struct_c.ts_pressure_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_temperature_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_temperature_grid_,
      blitz::shape(transfer_struct_c.ts_temperature_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_finegrid_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_finegrid_,
      blitz::shape(transfer_struct_c.ts_finegrid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_chapman transfer_struct_c;

  blitz::Array<double, 1> ts_height_grid_;
  blitz::Array<double, 1> ts_pressure_grid_;
  blitz::Array<double, 1> ts_temperature_grid_;
  blitz::Array<int, 1> ts_finegrid_;
  
};

// Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_optical_c_alloc_init(struct lidort_fixed_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_optical_c_init_only(struct lidort_fixed_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_optical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_optical_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_optical {
  double* ts_deltau_vert_input_;
  int ts_deltau_vert_input__f_shapes[2];
  int ts_deltau_vert_input__f_byte_size;

  double* ts_phasmoms_total_input_;
  int ts_phasmoms_total_input__f_shapes[3];
  int ts_phasmoms_total_input__f_byte_size;

  double* ts_thermal_bb_input_;
  int ts_thermal_bb_input__f_shapes[2];
  int ts_thermal_bb_input__f_byte_size;

  double* ts_lambertian_albedo_;
  int ts_lambertian_albedo__f_shapes[1];
  int ts_lambertian_albedo__f_byte_size;

  double* ts_surface_bb_input_;
  int ts_surface_bb_input__f_shapes[1];
  int ts_surface_bb_input__f_byte_size;

  
};

// Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Optical : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Optical() : Lidort_Structure() {
    lidort_fixed_optical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Optical(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_optical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Optical() {
    if (owns_pointer)
      lidort_fixed_optical_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 2>& ts_deltau_vert_input() const {
    return ts_deltau_vert_input_;
  }
  void ts_deltau_vert_input(const blitz::Array<double, 2>& ts_deltau_vert_input_in) {
    ts_deltau_vert_input_ = ts_deltau_vert_input_in;
  }
  
  const blitz::Array<double, 3>& ts_phasmoms_total_input() const {
    return ts_phasmoms_total_input_;
  }
  void ts_phasmoms_total_input(const blitz::Array<double, 3>& ts_phasmoms_total_input_in) {
    ts_phasmoms_total_input_ = ts_phasmoms_total_input_in;
  }
  
  const blitz::Array<double, 2>& ts_thermal_bb_input() const {
    return ts_thermal_bb_input_;
  }
  void ts_thermal_bb_input(const blitz::Array<double, 2>& ts_thermal_bb_input_in) {
    ts_thermal_bb_input_ = ts_thermal_bb_input_in;
  }
  
  const blitz::Array<double, 1>& ts_lambertian_albedo() const {
    return ts_lambertian_albedo_;
  }
  void ts_lambertian_albedo(const blitz::Array<double, 1>& ts_lambertian_albedo_in) {
    ts_lambertian_albedo_ = ts_lambertian_albedo_in;
  }
  
  const blitz::Array<double, 1>& ts_surface_bb_input() const {
    return ts_surface_bb_input_;
  }
  void ts_surface_bb_input(const blitz::Array<double, 1>& ts_surface_bb_input_in) {
    ts_surface_bb_input_ = ts_surface_bb_input_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Optical:" << std::endl
      << "   ts_deltau_vert_input: " << std::endl << ts_deltau_vert_input()  << std::endl
      << "ts_phasmoms_total_input: " << std::endl << ts_phasmoms_total_input()  << std::endl
      << "    ts_thermal_bb_input: " << std::endl << ts_thermal_bb_input()  << std::endl
      << "   ts_lambertian_albedo: " << std::endl << ts_lambertian_albedo()  << std::endl
      << "    ts_surface_bb_input: " << std::endl << ts_surface_bb_input()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_deltau_vert_input_",sizeof(*transfer_struct_c.ts_deltau_vert_input_),transfer_struct_c.ts_deltau_vert_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_phasmoms_total_input_",sizeof(*transfer_struct_c.ts_phasmoms_total_input_),transfer_struct_c.ts_phasmoms_total_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_thermal_bb_input_",sizeof(*transfer_struct_c.ts_thermal_bb_input_),transfer_struct_c.ts_thermal_bb_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lambertian_albedo_",sizeof(*transfer_struct_c.ts_lambertian_albedo_),transfer_struct_c.ts_lambertian_albedo__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_surface_bb_input_",sizeof(*transfer_struct_c.ts_surface_bb_input_),transfer_struct_c.ts_surface_bb_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_deltau_vert_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_deltau_vert_input_,
      blitz::shape(transfer_struct_c.ts_deltau_vert_input__f_shapes[0],
                   transfer_struct_c.ts_deltau_vert_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_phasmoms_total_input_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_phasmoms_total_input_,
      blitz::shape(transfer_struct_c.ts_phasmoms_total_input__f_shapes[0],
                   transfer_struct_c.ts_phasmoms_total_input__f_shapes[1],
                   transfer_struct_c.ts_phasmoms_total_input__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_thermal_bb_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_thermal_bb_input_,
      blitz::shape(transfer_struct_c.ts_thermal_bb_input__f_shapes[0],
                   transfer_struct_c.ts_thermal_bb_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_lambertian_albedo_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_lambertian_albedo_,
      blitz::shape(transfer_struct_c.ts_lambertian_albedo__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_surface_bb_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_surface_bb_input_,
      blitz::shape(transfer_struct_c.ts_surface_bb_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_optical transfer_struct_c;

  blitz::Array<double, 2> ts_deltau_vert_input_;
  blitz::Array<double, 3> ts_phasmoms_total_input_;
  blitz::Array<double, 2> ts_thermal_bb_input_;
  blitz::Array<double, 1> ts_lambertian_albedo_;
  blitz::Array<double, 1> ts_surface_bb_input_;
  
};

// Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_fixed_inputs_c_alloc_init(struct lidort_fixed_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_inputs_c_init_only(struct lidort_fixed_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_inputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_inputs {
  void* bool_;
  int bool__f_byte_size;

  void* cont_;
  int cont__f_byte_size;

  void* sunrays_;
  int sunrays__f_byte_size;

  void* userval_;
  int userval__f_byte_size;

  void* chapman_;
  int chapman__f_byte_size;

  void* optical_;
  int optical__f_byte_size;

  
};

// Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Fixed_Inputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Fixed_Inputs() : Lidort_Structure() {
    lidort_fixed_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Inputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_fixed_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Inputs() {
    if (owns_pointer)
      lidort_fixed_inputs_c_destroy(&fortran_type_c);
  }

  Lidort_Fixed_Boolean& f_bool() {
    return *bool_;
  }
  const Lidort_Fixed_Boolean& f_bool() const {
    return *bool_;
  }
  void f_bool(Lidort_Fixed_Boolean& bool_in) {
    void* src_ptr = bool_in.fortran_type_ptr();
    void* dst_ptr = bool_->fortran_type_ptr();
    lidort_fixed_boolean_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Control& cont() {
    return *cont_;
  }
  const Lidort_Fixed_Control& cont() const {
    return *cont_;
  }
  void cont(Lidort_Fixed_Control& cont_in) {
    void* src_ptr = cont_in.fortran_type_ptr();
    void* dst_ptr = cont_->fortran_type_ptr();
    lidort_fixed_control_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Sunrays& sunrays() {
    return *sunrays_;
  }
  const Lidort_Fixed_Sunrays& sunrays() const {
    return *sunrays_;
  }
  void sunrays(Lidort_Fixed_Sunrays& sunrays_in) {
    void* src_ptr = sunrays_in.fortran_type_ptr();
    void* dst_ptr = sunrays_->fortran_type_ptr();
    lidort_fixed_sunrays_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Uservalues& userval() {
    return *userval_;
  }
  const Lidort_Fixed_Uservalues& userval() const {
    return *userval_;
  }
  void userval(Lidort_Fixed_Uservalues& userval_in) {
    void* src_ptr = userval_in.fortran_type_ptr();
    void* dst_ptr = userval_->fortran_type_ptr();
    lidort_fixed_uservalues_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Chapman& chapman() {
    return *chapman_;
  }
  const Lidort_Fixed_Chapman& chapman() const {
    return *chapman_;
  }
  void chapman(Lidort_Fixed_Chapman& chapman_in) {
    void* src_ptr = chapman_in.fortran_type_ptr();
    void* dst_ptr = chapman_->fortran_type_ptr();
    lidort_fixed_chapman_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Fixed_Optical& optical() {
    return *optical_;
  }
  const Lidort_Fixed_Optical& optical() const {
    return *optical_;
  }
  void optical(Lidort_Fixed_Optical& optical_in) {
    void* src_ptr = optical_in.fortran_type_ptr();
    void* dst_ptr = optical_->fortran_type_ptr();
    lidort_fixed_optical_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Inputs:" << std::endl
      << "   bool: " << f_bool()  << std::endl
      << "   cont: " << cont()  << std::endl
      << "sunrays: " << sunrays()  << std::endl
      << "userval: " << userval()  << std::endl
      << "chapman: " << chapman()  << std::endl
      << "optical: " << optical()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    bool_.reset( new Lidort_Fixed_Boolean(transfer_struct_c.bool_) );
    cont_.reset( new Lidort_Fixed_Control(transfer_struct_c.cont_) );
    sunrays_.reset( new Lidort_Fixed_Sunrays(transfer_struct_c.sunrays_) );
    userval_.reset( new Lidort_Fixed_Uservalues(transfer_struct_c.userval_) );
    chapman_.reset( new Lidort_Fixed_Chapman(transfer_struct_c.chapman_) );
    optical_.reset( new Lidort_Fixed_Optical(transfer_struct_c.optical_) );
    
  }

  struct lidort_fixed_inputs transfer_struct_c;

  boost::shared_ptr<Lidort_Fixed_Boolean> bool_;
  boost::shared_ptr<Lidort_Fixed_Control> cont_;
  boost::shared_ptr<Lidort_Fixed_Sunrays> sunrays_;
  boost::shared_ptr<Lidort_Fixed_Uservalues> userval_;
  boost::shared_ptr<Lidort_Fixed_Chapman> chapman_;
  boost::shared_ptr<Lidort_Fixed_Optical> optical_;
  
};

// Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_boolean_c_alloc_init(struct lidort_modified_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_boolean_c_init_only(struct lidort_modified_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_boolean_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_boolean_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_boolean {
  int* ts_do_sscorr_nadir_;
  int ts_do_sscorr_nadir__f_byte_size;

  int* ts_do_sscorr_outgoing_;
  int ts_do_sscorr_outgoing__f_byte_size;

  int* ts_do_double_convtest_;
  int ts_do_double_convtest__f_byte_size;

  int* ts_do_solar_sources_;
  int ts_do_solar_sources__f_byte_size;

  int* ts_do_refractive_geometry_;
  int ts_do_refractive_geometry__f_byte_size;

  int* ts_do_chapman_function_;
  int ts_do_chapman_function__f_byte_size;

  int* ts_do_rayleigh_only_;
  int ts_do_rayleigh_only__f_byte_size;

  int* ts_do_isotropic_only_;
  int ts_do_isotropic_only__f_byte_size;

  int* ts_do_no_azimuth_;
  int ts_do_no_azimuth__f_byte_size;

  int* ts_do_all_fourier_;
  int ts_do_all_fourier__f_byte_size;

  int* ts_do_deltam_scaling_;
  int ts_do_deltam_scaling__f_byte_size;

  int* ts_do_solution_saving_;
  int ts_do_solution_saving__f_byte_size;

  int* ts_do_bvp_telescoping_;
  int ts_do_bvp_telescoping__f_byte_size;

  int* ts_do_user_streams_;
  int ts_do_user_streams__f_byte_size;

  int* ts_do_additional_mvout_;
  int ts_do_additional_mvout__f_byte_size;

  int* ts_do_mvout_only_;
  int ts_do_mvout_only__f_byte_size;

  int* ts_do_thermal_transonly_;
  int ts_do_thermal_transonly__f_byte_size;

  
};

// Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Boolean : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Boolean() : Lidort_Structure() {
    lidort_modified_boolean_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Boolean(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_boolean_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Boolean() {
    if (owns_pointer)
      lidort_modified_boolean_c_destroy(&fortran_type_c);
  }

  const bool ts_do_sscorr_nadir() const {
    return *transfer_struct_c.ts_do_sscorr_nadir_ != 0;
  }
  void ts_do_sscorr_nadir(const bool& ts_do_sscorr_nadir_in) {
    *transfer_struct_c.ts_do_sscorr_nadir_ = ts_do_sscorr_nadir_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_sscorr_outgoing() const {
    return *transfer_struct_c.ts_do_sscorr_outgoing_ != 0;
  }
  void ts_do_sscorr_outgoing(const bool& ts_do_sscorr_outgoing_in) {
    *transfer_struct_c.ts_do_sscorr_outgoing_ = ts_do_sscorr_outgoing_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_double_convtest() const {
    return *transfer_struct_c.ts_do_double_convtest_ != 0;
  }
  void ts_do_double_convtest(const bool& ts_do_double_convtest_in) {
    *transfer_struct_c.ts_do_double_convtest_ = ts_do_double_convtest_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_solar_sources() const {
    return *transfer_struct_c.ts_do_solar_sources_ != 0;
  }
  void ts_do_solar_sources(const bool& ts_do_solar_sources_in) {
    *transfer_struct_c.ts_do_solar_sources_ = ts_do_solar_sources_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_refractive_geometry() const {
    return *transfer_struct_c.ts_do_refractive_geometry_ != 0;
  }
  void ts_do_refractive_geometry(const bool& ts_do_refractive_geometry_in) {
    *transfer_struct_c.ts_do_refractive_geometry_ = ts_do_refractive_geometry_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_chapman_function() const {
    return *transfer_struct_c.ts_do_chapman_function_ != 0;
  }
  void ts_do_chapman_function(const bool& ts_do_chapman_function_in) {
    *transfer_struct_c.ts_do_chapman_function_ = ts_do_chapman_function_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_rayleigh_only() const {
    return *transfer_struct_c.ts_do_rayleigh_only_ != 0;
  }
  void ts_do_rayleigh_only(const bool& ts_do_rayleigh_only_in) {
    *transfer_struct_c.ts_do_rayleigh_only_ = ts_do_rayleigh_only_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_isotropic_only() const {
    return *transfer_struct_c.ts_do_isotropic_only_ != 0;
  }
  void ts_do_isotropic_only(const bool& ts_do_isotropic_only_in) {
    *transfer_struct_c.ts_do_isotropic_only_ = ts_do_isotropic_only_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_no_azimuth() const {
    return *transfer_struct_c.ts_do_no_azimuth_ != 0;
  }
  void ts_do_no_azimuth(const bool& ts_do_no_azimuth_in) {
    *transfer_struct_c.ts_do_no_azimuth_ = ts_do_no_azimuth_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_all_fourier() const {
    return *transfer_struct_c.ts_do_all_fourier_ != 0;
  }
  void ts_do_all_fourier(const bool& ts_do_all_fourier_in) {
    *transfer_struct_c.ts_do_all_fourier_ = ts_do_all_fourier_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_deltam_scaling() const {
    return *transfer_struct_c.ts_do_deltam_scaling_ != 0;
  }
  void ts_do_deltam_scaling(const bool& ts_do_deltam_scaling_in) {
    *transfer_struct_c.ts_do_deltam_scaling_ = ts_do_deltam_scaling_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_solution_saving() const {
    return *transfer_struct_c.ts_do_solution_saving_ != 0;
  }
  void ts_do_solution_saving(const bool& ts_do_solution_saving_in) {
    *transfer_struct_c.ts_do_solution_saving_ = ts_do_solution_saving_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_bvp_telescoping() const {
    return *transfer_struct_c.ts_do_bvp_telescoping_ != 0;
  }
  void ts_do_bvp_telescoping(const bool& ts_do_bvp_telescoping_in) {
    *transfer_struct_c.ts_do_bvp_telescoping_ = ts_do_bvp_telescoping_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_user_streams() const {
    return *transfer_struct_c.ts_do_user_streams_ != 0;
  }
  void ts_do_user_streams(const bool& ts_do_user_streams_in) {
    *transfer_struct_c.ts_do_user_streams_ = ts_do_user_streams_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_additional_mvout() const {
    return *transfer_struct_c.ts_do_additional_mvout_ != 0;
  }
  void ts_do_additional_mvout(const bool& ts_do_additional_mvout_in) {
    *transfer_struct_c.ts_do_additional_mvout_ = ts_do_additional_mvout_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_mvout_only() const {
    return *transfer_struct_c.ts_do_mvout_only_ != 0;
  }
  void ts_do_mvout_only(const bool& ts_do_mvout_only_in) {
    *transfer_struct_c.ts_do_mvout_only_ = ts_do_mvout_only_in ? FORTRAN_TRUE_INT : 0;
  }
  
  const bool ts_do_thermal_transonly() const {
    return *transfer_struct_c.ts_do_thermal_transonly_ != 0;
  }
  void ts_do_thermal_transonly(const bool& ts_do_thermal_transonly_in) {
    *transfer_struct_c.ts_do_thermal_transonly_ = ts_do_thermal_transonly_in ? FORTRAN_TRUE_INT : 0;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Boolean:" << std::endl
      << "       ts_do_sscorr_nadir: " << ts_do_sscorr_nadir()  << std::endl
      << "    ts_do_sscorr_outgoing: " << ts_do_sscorr_outgoing()  << std::endl
      << "    ts_do_double_convtest: " << ts_do_double_convtest()  << std::endl
      << "      ts_do_solar_sources: " << ts_do_solar_sources()  << std::endl
      << "ts_do_refractive_geometry: " << ts_do_refractive_geometry()  << std::endl
      << "   ts_do_chapman_function: " << ts_do_chapman_function()  << std::endl
      << "      ts_do_rayleigh_only: " << ts_do_rayleigh_only()  << std::endl
      << "     ts_do_isotropic_only: " << ts_do_isotropic_only()  << std::endl
      << "         ts_do_no_azimuth: " << ts_do_no_azimuth()  << std::endl
      << "        ts_do_all_fourier: " << ts_do_all_fourier()  << std::endl
      << "     ts_do_deltam_scaling: " << ts_do_deltam_scaling()  << std::endl
      << "    ts_do_solution_saving: " << ts_do_solution_saving()  << std::endl
      << "    ts_do_bvp_telescoping: " << ts_do_bvp_telescoping()  << std::endl
      << "       ts_do_user_streams: " << ts_do_user_streams()  << std::endl
      << "   ts_do_additional_mvout: " << ts_do_additional_mvout()  << std::endl
      << "         ts_do_mvout_only: " << ts_do_mvout_only()  << std::endl
      << "  ts_do_thermal_transonly: " << ts_do_thermal_transonly()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_sscorr_nadir_",sizeof(*transfer_struct_c.ts_do_sscorr_nadir_),transfer_struct_c.ts_do_sscorr_nadir__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sscorr_outgoing_",sizeof(*transfer_struct_c.ts_do_sscorr_outgoing_),transfer_struct_c.ts_do_sscorr_outgoing__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_double_convtest_",sizeof(*transfer_struct_c.ts_do_double_convtest_),transfer_struct_c.ts_do_double_convtest__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_solar_sources_",sizeof(*transfer_struct_c.ts_do_solar_sources_),transfer_struct_c.ts_do_solar_sources__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_refractive_geometry_",sizeof(*transfer_struct_c.ts_do_refractive_geometry_),transfer_struct_c.ts_do_refractive_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_chapman_function_",sizeof(*transfer_struct_c.ts_do_chapman_function_),transfer_struct_c.ts_do_chapman_function__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_rayleigh_only_",sizeof(*transfer_struct_c.ts_do_rayleigh_only_),transfer_struct_c.ts_do_rayleigh_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_isotropic_only_",sizeof(*transfer_struct_c.ts_do_isotropic_only_),transfer_struct_c.ts_do_isotropic_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_no_azimuth_",sizeof(*transfer_struct_c.ts_do_no_azimuth_),transfer_struct_c.ts_do_no_azimuth__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_all_fourier_",sizeof(*transfer_struct_c.ts_do_all_fourier_),transfer_struct_c.ts_do_all_fourier__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_deltam_scaling_",sizeof(*transfer_struct_c.ts_do_deltam_scaling_),transfer_struct_c.ts_do_deltam_scaling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_solution_saving_",sizeof(*transfer_struct_c.ts_do_solution_saving_),transfer_struct_c.ts_do_solution_saving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_bvp_telescoping_",sizeof(*transfer_struct_c.ts_do_bvp_telescoping_),transfer_struct_c.ts_do_bvp_telescoping__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_user_streams_",sizeof(*transfer_struct_c.ts_do_user_streams_),transfer_struct_c.ts_do_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_additional_mvout_",sizeof(*transfer_struct_c.ts_do_additional_mvout_),transfer_struct_c.ts_do_additional_mvout__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_mvout_only_",sizeof(*transfer_struct_c.ts_do_mvout_only_),transfer_struct_c.ts_do_mvout_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_thermal_transonly_",sizeof(*transfer_struct_c.ts_do_thermal_transonly_),transfer_struct_c.ts_do_thermal_transonly__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_boolean transfer_struct_c;

  
};

// Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_control_c_alloc_init(struct lidort_modified_control *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_control_c_init_only(struct lidort_modified_control *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_control_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_control_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_control {
  int* ts_nmoments_input_;
  int ts_nmoments_input__f_byte_size;

  
};

// Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Control : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Control() : Lidort_Structure() {
    lidort_modified_control_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Control(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_control_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Control() {
    if (owns_pointer)
      lidort_modified_control_c_destroy(&fortran_type_c);
  }

  const int& ts_nmoments_input() const {
    return *transfer_struct_c.ts_nmoments_input_;
  }
  void ts_nmoments_input(const int& ts_nmoments_input_in) {
    *transfer_struct_c.ts_nmoments_input_ = ts_nmoments_input_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Control:" << std::endl
      << "ts_nmoments_input: " << ts_nmoments_input()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_nmoments_input_",sizeof(*transfer_struct_c.ts_nmoments_input_),transfer_struct_c.ts_nmoments_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_control transfer_struct_c;

  
};

// Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_sunrays_c_alloc_init(struct lidort_modified_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_sunrays_c_init_only(struct lidort_modified_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_sunrays_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_sunrays_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_sunrays {
  int* ts_nbeams_;
  int ts_nbeams__f_byte_size;

  double* ts_beam_szas_;
  int ts_beam_szas__f_shapes[1];
  int ts_beam_szas__f_byte_size;

  
};

// Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Sunrays : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Sunrays() : Lidort_Structure() {
    lidort_modified_sunrays_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Sunrays(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_sunrays_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Sunrays() {
    if (owns_pointer)
      lidort_modified_sunrays_c_destroy(&fortran_type_c);
  }

  const int& ts_nbeams() const {
    return *transfer_struct_c.ts_nbeams_;
  }
  void ts_nbeams(const int& ts_nbeams_in) {
    *transfer_struct_c.ts_nbeams_ = ts_nbeams_in;
  }
  
  const blitz::Array<double, 1>& ts_beam_szas() const {
    return ts_beam_szas_;
  }
  void ts_beam_szas(const blitz::Array<double, 1>& ts_beam_szas_in) {
    ts_beam_szas_ = ts_beam_szas_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Sunrays:" << std::endl
      << "   ts_nbeams: " << ts_nbeams()  << std::endl
      << "ts_beam_szas: " << std::endl << ts_beam_szas()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_nbeams_",sizeof(*transfer_struct_c.ts_nbeams_),transfer_struct_c.ts_nbeams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_beam_szas_",sizeof(*transfer_struct_c.ts_beam_szas_),transfer_struct_c.ts_beam_szas__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_beam_szas_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_beam_szas_,
      blitz::shape(transfer_struct_c.ts_beam_szas__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_sunrays transfer_struct_c;

  blitz::Array<double, 1> ts_beam_szas_;
  
};

// Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_uservalues_c_alloc_init(struct lidort_modified_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_uservalues_c_init_only(struct lidort_modified_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_uservalues_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_uservalues_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_uservalues {
  int* ts_n_user_relazms_;
  int ts_n_user_relazms__f_byte_size;

  double* ts_user_relazms_;
  int ts_user_relazms__f_shapes[1];
  int ts_user_relazms__f_byte_size;

  double* ts_user_angles_input_;
  int ts_user_angles_input__f_shapes[1];
  int ts_user_angles_input__f_byte_size;

  double* ts_user_levels_;
  int ts_user_levels__f_shapes[1];
  int ts_user_levels__f_byte_size;

  double* ts_geometry_specheight_;
  int ts_geometry_specheight__f_byte_size;

  
};

// Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Uservalues : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Uservalues() : Lidort_Structure() {
    lidort_modified_uservalues_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Uservalues(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_uservalues_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Uservalues() {
    if (owns_pointer)
      lidort_modified_uservalues_c_destroy(&fortran_type_c);
  }

  const int& ts_n_user_relazms() const {
    return *transfer_struct_c.ts_n_user_relazms_;
  }
  void ts_n_user_relazms(const int& ts_n_user_relazms_in) {
    *transfer_struct_c.ts_n_user_relazms_ = ts_n_user_relazms_in;
  }
  
  const blitz::Array<double, 1>& ts_user_relazms() const {
    return ts_user_relazms_;
  }
  void ts_user_relazms(const blitz::Array<double, 1>& ts_user_relazms_in) {
    ts_user_relazms_ = ts_user_relazms_in;
  }
  
  const blitz::Array<double, 1>& ts_user_angles_input() const {
    return ts_user_angles_input_;
  }
  void ts_user_angles_input(const blitz::Array<double, 1>& ts_user_angles_input_in) {
    ts_user_angles_input_ = ts_user_angles_input_in;
  }
  
  const blitz::Array<double, 1>& ts_user_levels() const {
    return ts_user_levels_;
  }
  void ts_user_levels(const blitz::Array<double, 1>& ts_user_levels_in) {
    ts_user_levels_ = ts_user_levels_in;
  }
  
  const double& ts_geometry_specheight() const {
    return *transfer_struct_c.ts_geometry_specheight_;
  }
  void ts_geometry_specheight(const double& ts_geometry_specheight_in) {
    *transfer_struct_c.ts_geometry_specheight_ = ts_geometry_specheight_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Uservalues:" << std::endl
      << "     ts_n_user_relazms: " << ts_n_user_relazms()  << std::endl
      << "       ts_user_relazms: " << std::endl << ts_user_relazms()  << std::endl
      << "  ts_user_angles_input: " << std::endl << ts_user_angles_input()  << std::endl
      << "        ts_user_levels: " << std::endl << ts_user_levels()  << std::endl
      << "ts_geometry_specheight: " << ts_geometry_specheight()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_n_user_relazms_",sizeof(*transfer_struct_c.ts_n_user_relazms_),transfer_struct_c.ts_n_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_relazms_",sizeof(*transfer_struct_c.ts_user_relazms_),transfer_struct_c.ts_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_angles_input_",sizeof(*transfer_struct_c.ts_user_angles_input_),transfer_struct_c.ts_user_angles_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_levels_",sizeof(*transfer_struct_c.ts_user_levels_),transfer_struct_c.ts_user_levels__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_geometry_specheight_",sizeof(*transfer_struct_c.ts_geometry_specheight_),transfer_struct_c.ts_geometry_specheight__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_user_relazms_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_relazms_,
      blitz::shape(transfer_struct_c.ts_user_relazms__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_user_angles_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_angles_input_,
      blitz::shape(transfer_struct_c.ts_user_angles_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_user_levels_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_levels_,
      blitz::shape(transfer_struct_c.ts_user_levels__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_uservalues transfer_struct_c;

  blitz::Array<double, 1> ts_user_relazms_;
  blitz::Array<double, 1> ts_user_angles_input_;
  blitz::Array<double, 1> ts_user_levels_;
  
};

// Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_chapman_c_alloc_init(struct lidort_modified_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_chapman_c_init_only(struct lidort_modified_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_chapman_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_chapman_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_chapman {
  double* ts_earth_radius_;
  int ts_earth_radius__f_byte_size;

  
};

// Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Chapman : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Chapman() : Lidort_Structure() {
    lidort_modified_chapman_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Chapman(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_chapman_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Chapman() {
    if (owns_pointer)
      lidort_modified_chapman_c_destroy(&fortran_type_c);
  }

  const double& ts_earth_radius() const {
    return *transfer_struct_c.ts_earth_radius_;
  }
  void ts_earth_radius(const double& ts_earth_radius_in) {
    *transfer_struct_c.ts_earth_radius_ = ts_earth_radius_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Chapman:" << std::endl
      << "ts_earth_radius: " << ts_earth_radius()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_earth_radius_",sizeof(*transfer_struct_c.ts_earth_radius_),transfer_struct_c.ts_earth_radius__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_chapman transfer_struct_c;

  
};

// Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_optical_c_alloc_init(struct lidort_modified_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_optical_c_init_only(struct lidort_modified_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_optical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_optical_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_optical {
  double* ts_omega_total_input_;
  int ts_omega_total_input__f_shapes[2];
  int ts_omega_total_input__f_byte_size;

  
};

// Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Optical : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Optical() : Lidort_Structure() {
    lidort_modified_optical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Optical(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_optical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Optical() {
    if (owns_pointer)
      lidort_modified_optical_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 2>& ts_omega_total_input() const {
    return ts_omega_total_input_;
  }
  void ts_omega_total_input(const blitz::Array<double, 2>& ts_omega_total_input_in) {
    ts_omega_total_input_ = ts_omega_total_input_in;
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Optical:" << std::endl
      << "ts_omega_total_input: " << std::endl << ts_omega_total_input()  << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_omega_total_input_",sizeof(*transfer_struct_c.ts_omega_total_input_),transfer_struct_c.ts_omega_total_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_omega_total_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_omega_total_input_,
      blitz::shape(transfer_struct_c.ts_omega_total_input__f_shapes[0],
                   transfer_struct_c.ts_omega_total_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_optical transfer_struct_c;

  blitz::Array<double, 2> ts_omega_total_input_;
  
};

// Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
extern "C" {
  void lidort_modified_inputs_c_alloc_init(struct lidort_modified_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_inputs_c_init_only(struct lidort_modified_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_inputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_inputs {
  void* mbool_;
  int mbool__f_byte_size;

  void* mcont_;
  int mcont__f_byte_size;

  void* msunrays_;
  int msunrays__f_byte_size;

  void* muserval_;
  int muserval__f_byte_size;

  void* mchapman_;
  int mchapman__f_byte_size;

  void* moptical_;
  int moptical__f_byte_size;

  
};

// Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
class Lidort_Modified_Inputs : public Lidort_Structure {
public:
  // Allocating constructor
  Lidort_Modified_Inputs() : Lidort_Structure() {
    lidort_modified_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Inputs(void* allocated_f_type_c) :  Lidort_Structure(allocated_f_type_c) {
    lidort_modified_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Inputs() {
    if (owns_pointer)
      lidort_modified_inputs_c_destroy(&fortran_type_c);
  }

  Lidort_Modified_Boolean& mbool() {
    return *mbool_;
  }
  const Lidort_Modified_Boolean& mbool() const {
    return *mbool_;
  }
  void mbool(Lidort_Modified_Boolean& mbool_in) {
    void* src_ptr = mbool_in.fortran_type_ptr();
    void* dst_ptr = mbool_->fortran_type_ptr();
    lidort_modified_boolean_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Modified_Control& mcont() {
    return *mcont_;
  }
  const Lidort_Modified_Control& mcont() const {
    return *mcont_;
  }
  void mcont(Lidort_Modified_Control& mcont_in) {
    void* src_ptr = mcont_in.fortran_type_ptr();
    void* dst_ptr = mcont_->fortran_type_ptr();
    lidort_modified_control_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Modified_Sunrays& msunrays() {
    return *msunrays_;
  }
  const Lidort_Modified_Sunrays& msunrays() const {
    return *msunrays_;
  }
  void msunrays(Lidort_Modified_Sunrays& msunrays_in) {
    void* src_ptr = msunrays_in.fortran_type_ptr();
    void* dst_ptr = msunrays_->fortran_type_ptr();
    lidort_modified_sunrays_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Modified_Uservalues& muserval() {
    return *muserval_;
  }
  const Lidort_Modified_Uservalues& muserval() const {
    return *muserval_;
  }
  void muserval(Lidort_Modified_Uservalues& muserval_in) {
    void* src_ptr = muserval_in.fortran_type_ptr();
    void* dst_ptr = muserval_->fortran_type_ptr();
    lidort_modified_uservalues_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Modified_Chapman& mchapman() {
    return *mchapman_;
  }
  const Lidort_Modified_Chapman& mchapman() const {
    return *mchapman_;
  }
  void mchapman(Lidort_Modified_Chapman& mchapman_in) {
    void* src_ptr = mchapman_in.fortran_type_ptr();
    void* dst_ptr = mchapman_->fortran_type_ptr();
    lidort_modified_chapman_c_copy(&src_ptr, &dst_ptr);
  }
  
  Lidort_Modified_Optical& moptical() {
    return *moptical_;
  }
  const Lidort_Modified_Optical& moptical() const {
    return *moptical_;
  }
  void moptical(Lidort_Modified_Optical& moptical_in) {
    void* src_ptr = moptical_in.fortran_type_ptr();
    void* dst_ptr = moptical_->fortran_type_ptr();
    lidort_modified_optical_c_copy(&src_ptr, &dst_ptr);
  }
  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Inputs:" << std::endl
      << "   mbool: " << mbool()  << std::endl
      << "   mcont: " << mcont()  << std::endl
      << "msunrays: " << msunrays()  << std::endl
      << "muserval: " << muserval()  << std::endl
      << "mchapman: " << mchapman()  << std::endl
      << "moptical: " << moptical()  << std::endl;
  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    mbool_.reset( new Lidort_Modified_Boolean(transfer_struct_c.mbool_) );
    mcont_.reset( new Lidort_Modified_Control(transfer_struct_c.mcont_) );
    msunrays_.reset( new Lidort_Modified_Sunrays(transfer_struct_c.msunrays_) );
    muserval_.reset( new Lidort_Modified_Uservalues(transfer_struct_c.muserval_) );
    mchapman_.reset( new Lidort_Modified_Chapman(transfer_struct_c.mchapman_) );
    moptical_.reset( new Lidort_Modified_Optical(transfer_struct_c.moptical_) );
    
  }

  struct lidort_modified_inputs transfer_struct_c;

  boost::shared_ptr<Lidort_Modified_Boolean> mbool_;
  boost::shared_ptr<Lidort_Modified_Control> mcont_;
  boost::shared_ptr<Lidort_Modified_Sunrays> msunrays_;
  boost::shared_ptr<Lidort_Modified_Uservalues> muserval_;
  boost::shared_ptr<Lidort_Modified_Chapman> mchapman_;
  boost::shared_ptr<Lidort_Modified_Optical> moptical_;
  
};



}
#endif