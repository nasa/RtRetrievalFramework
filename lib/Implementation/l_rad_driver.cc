#include "l_rad_driver.h"

#include "wgs84_constant.h"
#include "linear_algebra.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

// l_rad Fortran code functions
extern "C" {
    void lr_init(const int* nstream,
                 const int* nstokes, const int* surftype,
                 void** l_rad_struct_c);
    void lr_cleanup(void** l_rad_struct_c);
    void init_l_rad(void** l_rad_struct_c, const double* sza,
                    const double* zen, const double* azm, const double* alt,
                    const double* radius_km, const int* nlayer, 
                    const bool* regular_ps, const bool* enhanced_ps, const bool* pure_nadir);
    void l_rad_first_driver(void* const* l_rad_struct_c, const int* nlayer,
                            const int* natm_der,
                            const int* nstokes, const double* tau,
                            const double* L_tau, const double* omega,
                            const double* L_omega,
                            const double* Zmat, const double* L_Zmat,
                            const double* fscale, const double* L_fscale,
                            const double* spars, const int* nspars,
                            const int* need_jacobians_i,
                            double* R1, double* L_R1,
                            double* Ls_R1);

    void l_rad_second_driver(void* const* l_rad_struct_c, const int* n_layers,
                             const int* n_params, const int* nstokes, const int* n_streams,
                             const int*n_scatt, const double* tau, const double* L_tau,
                             const double* omega, const double* L_omega,
                             const double* spars, const int* nspars,
                             const double* coefs, const double* L_coefs, const int* n_eff_coefs,
                             const int* need_jacobians_i,
                             double* R2, double* L_R2, double* Ls_R2,
                             double* ICorr, double* L_ICorr, double* Ls_ICorr);

    void calc_z(void* const* l_rad_struct_c, const int* nlay, const int* nstokes,
                const int* nmom, const int* npar, const double* dcoef,
                const double* l_dcoef, double* zmat, double* l_zmat);
}

LRadDriver::LRadDriver(int Number_stream, int Number_stokes,
                       int Surface_type,
                       bool Tms_correction,
                       bool Pure_nadir, const PsMode ps_mode)
    : nstream(Number_stream),
      nstokes(Number_stokes),
      surface_type(Surface_type),
      use_tms_correction(Tms_correction),
      pure_nadir(Pure_nadir)
{
    // There seems to a bug in l_rad if nstokes = 4,
    // it will not compute the third stokes value,
    // so for now the maximum value is 3
    range_check(nstokes, 1, 4);

    initialize(ps_mode);
}

//-----------------------------------------------------------------------
/// Finish up initialization.
//-----------------------------------------------------------------------

void LRadDriver::initialize(const PsMode ps_mode)
{
    range_min_check(nstream, 1);

    // Determine which pseudo-spherical mode to use based on the zenith angle
    // ie, how close are we to nadir, or use preset value
    set_ps_mode(ps_mode);

    // Initially we are not using jacobians until linear inputs are setup
    need_jacobians_i = 0;

    lr_init(&nstream, &nstokes, &surface_type, &l_rad_struct);
}

//-----------------------------------------------------------------------
/// Destructor
//-----------------------------------------------------------------------

LRadDriver::~LRadDriver()
{
    lr_cleanup(&l_rad_struct);
}

//-----------------------------------------------------------------------
/// Update the altitude information. This can change the number of
/// layers if desired.
//-----------------------------------------------------------------------

void LRadDriver::setup_geometry(Array<double, 1> alt, double sza, double zen, double azm) const
{
    int nlayer = alt.rows() - 1;

    double radius_km = OldConstant::wgs84_a.convert(units::km).value;

    init_l_rad(&l_rad_struct, &sza, &zen, &azm, alt.dataFirst(),
               &radius_km, &nlayer, &regular_ps, &enhanced_ps, &pure_nadir);
}

//-----------------------------------------------------------------------
/// Calculate the z matrix. This does the full calculation, and is
/// somewhat expensive to call.
///
/// Must be called after setup_geometry since this calculation uses
/// geometry values setup from it. 
//-----------------------------------------------------------------------

ArrayAd<double, 2>
LRadDriver::z_matrix(const ArrayAd<double, 3>& pf) const
{
    int nlayer = pf.cols();
    int natm_jac = pf.number_variable();

    ArrayAd<double, 2> zmat(nlayer, nstokes, natm_jac);
    zmat.value().reference
    (Array<double, 2>(nlayer, nstokes, ColumnMajorArray<2>()));
    zmat.jacobian().reference
    (Array<double, 3>(nlayer, nstokes, natm_jac, ColumnMajorArray<3>()));
    Array<double, 3> pf_f(to_fortran(pf.value()));
    Array<double, 4> l_pf_f(to_fortran(pf.jacobian()));
    int nmom = pf.rows() - 1;
    calc_z(&l_rad_struct, &nlayer, &nstokes, &nmom, &natm_jac, pf_f.dataFirst(),
           l_pf_f.dataFirst(), zmat.value().dataFirst(),
           zmat.jacobian().dataFirst());
    return zmat;
}

void LRadDriver::setup_surface_params(const Array<double, 1>& surface_param)
{
    if (surface_param_f.rows() != surface_param.rows()) {
        surface_param_f.reference(Array<double, 1>(surface_param.rows(), ColumnMajorArray<1>()));
    }

    // A hack to deal with the different ordering.  There should be a better
    // place to do this.
    if (surface_type != 2)
        surface_param_f = surface_param;
    else {
        surface_param_f(0) = surface_param(1);
        surface_param_f(1) = surface_param(2);
        surface_param_f(2) = surface_param(3);
        surface_param_f(3) = surface_param(4);
        surface_param_f(4) = surface_param(0);
    }
}

void LRadDriver::setup_optical_inputs(const Array<double, 1>& od, 
                                      const Array<double, 1>& ssa,
                                      const Array<double, 3>& pf,
                                      const Array<double, 2>& zmat)
{

    if (tau_f.rows() != od.rows()) {
        tau_f.reference(Array<double, 1>(od.shape(), ColumnMajorArray<1>()));
    }
    tau_f = od;

    if (omega_f.rows() != ssa.rows()) {
        omega_f.reference(Array<double, 1>(ssa.shape(), ColumnMajorArray<1>()));
    }
    omega_f = ssa;

    if(pf_f.rows() != pf.rows() or pf_f.cols() != pf.cols() or pf_f.depth() != pf.depth()) {
        pf_f.reference(Array<double, 3>(pf.shape(), ColumnMajorArray<3>()));
    }
    pf_f = pf;

    if (zmat_f.rows() != zmat.rows() or zmat_f.cols() != zmat.cols()) {
        zmat_f.reference(Array<double, 2>(zmat.shape(), ColumnMajorArray<2>()));
    }
    zmat_f = zmat;

    // Calculate factor needed for TMS correction. This is the
    // correction needed when we are combining with a multi-scattering
    // Corrects for errors due to delta-M scaling in intensity RT.
    // Compute single scattering using ALL themoments of the phase function
    // but using the truncated phase function to compute the multiple scattering
    // contribution.
    //
    // If we aren't using an underlying multiscattering RT, then this
    // should be zero.
    int nmom = 2 * nstream;

    if(fscale_f.rows() != od.rows()) {
        fscale_f.reference(Array<double, 1>(od.shape(), ColumnMajorArray<1>()));
    }

    if (use_tms_correction) {
        fscale_f = pf(nmom, Range::all(), 0) / (2 * nmom + 1);
    } else {
        fscale_f = 0;
    }

}

void LRadDriver::clear_linear_inputs()
{
    // We no longer need jacobians
    need_jacobians_i = 0;
        
    // These need to be resized to what the fortran interface is expecting
    // Also zero them out
    if(jac_atm_f.rows() != number_stokes() or jac_atm_f.cols() != tau_f.rows()) {
        jac_atm_f.resize(number_stokes(), tau_f.rows(), 1);
    }
    jac_atm_f = 0;

    if (jac_surf_f.rows() != number_stokes() or jac_surf_f.cols() != surface_param_f.rows()) {
        jac_surf_f.resize(number_stokes(), surface_param_f.rows());
    }
    jac_surf_f = 0;

    if(l_tau_f.rows() != tau_f.rows()) {
        l_tau_f.resize(tau_f.rows(), 1);
    }
    l_tau_f = 0;

    if(l_omega_f.rows() != omega_f.rows()) {
        l_omega_f.resize(omega_f.rows(), 1);
    }
    l_omega_f = 0;

    if(l_pf_f.rows() != pf_f.rows() or l_pf_f.cols() != pf_f.cols() or l_pf_f.depth() != pf_f.depth()) {
        l_pf_f.resize(pf_f.rows(), pf_f.cols(), pf_f.depth(), 1);
    }
    l_pf_f = 0;

    if(l_zmat_f.rows() != zmat_f.rows() or l_zmat_f.cols() != zmat_f.cols()) {
        l_zmat_f.resize(zmat_f.rows(), zmat_f.cols(), 1);
    }
    l_zmat_f = 0;

    if(l_fscale_f.rows() != fscale_f.rows()) {
        l_fscale_f.resize(fscale_f.rows(), 1);
    }
    l_fscale_f = 0;
}

void LRadDriver::setup_linear_inputs(const ArrayAd<double, 1>& od,
                                     const ArrayAd<double, 1>& ssa,
                                     const ArrayAd<double, 3>& pf,
                                     const ArrayAd<double, 2>& zmat)
{
    // Flag that we will be requesting jacobians
    need_jacobians_i = 1;

    int natm_jac = std::max(pf.number_variable(),
                            std::max(od.number_variable(),
                                     ssa.number_variable()));

    // Initialize output variables for jacobians
    if(jac_atm_f.rows() != number_stokes() or jac_atm_f.cols() != od.rows() or jac_atm_f.depth() != natm_jac) {
        jac_atm_f.reference(Array<double, 3>(number_stokes(), od.rows(), natm_jac, ColumnMajorArray<3>()));
    }

    if (jac_surf_f.rows() != number_stokes() or jac_surf_f.cols() != surface_param_f.rows()) {
        jac_surf_f.reference(Array<double, 2>(number_stokes(), surface_param_f.rows(), ColumnMajorArray<2>()));
    }

    // Copy over input jacobian arrays
    if(l_tau_f.rows() != od.jacobian().rows() or l_tau_f.cols() != od.jacobian().cols()) {
        l_tau_f.reference(Array<double, 2>(od.jacobian().shape(), ColumnMajorArray<2>()));
    }
    l_tau_f = od.jacobian();

    if(l_omega_f.rows() != ssa.jacobian().rows() or l_omega_f.cols() != ssa.jacobian().cols()) {
        l_omega_f.reference(Array<double, 2>(ssa.jacobian().shape(), ColumnMajorArray<2>()));
    }
    l_omega_f = ssa.jacobian();

    if(l_pf_f.rows() != pf.jacobian().rows() or l_pf_f.cols() != pf.jacobian().cols() or l_pf_f.depth() != pf.jacobian().extent(fourthDim) or l_pf_f.extent(fourthDim)) {
        l_pf_f.reference(Array<double, 4>(pf.jacobian().shape(), ColumnMajorArray<4>()));
    }
    l_pf_f = pf.jacobian();

    if(l_zmat_f.rows() != zmat.jacobian().rows() or l_zmat_f.cols() != zmat.jacobian().cols() or l_zmat_f.depth() != zmat.jacobian().depth()) {
        l_zmat_f.reference(Array<double, 3>(zmat.jacobian().shape(), ColumnMajorArray<3>()));
    }
    l_zmat_f = zmat.jacobian();

    // Calculate factor needed for TMS correction. This is the
    int nmom = 2 * nstream;
    Range ra(Range::all());

    if(l_fscale_f.rows() != od.rows() or l_fscale_f.cols() != natm_jac) {
        l_fscale_f.reference(Array<double, 2>(od.rows(), natm_jac, ColumnMajorArray<2>()));
    }

    if (use_tms_correction) {
        l_fscale_f = pf.jacobian()(nmom, ra, 0, ra) / (2 * nmom + 1);
    } else {
        l_fscale_f = 0;
    }

}

/// Check that inputs to RT methods are allocated and sized correctly
void LRadDriver::check_rt_inputs()
{
    if(surface_param_f.rows() == 0) {
       throw Exception("surface_param_f not allocated"); 
    }

    if(tau_f.rows() == 0) {
        throw Exception("tau_f not allocated");
    }

    if(omega_f.rows() == 0) {
       throw Exception("omega_f not allocated");
    }

    if(zmat_f.rows() == 0) {
        throw Exception("zmat_f not allocated");
    }

    if(fscale_f.rows() == 0) {
       throw Exception("fscale_f not allocated");
    }

    if (need_jacobians_i == 1) {
        if(l_tau_f.rows() == 0) {
            throw Exception("l_tau_f not allocated");
        }

        if(l_omega_f.rows() == 0) {
            throw Exception("l_omega_f not allocated");
        }

        if(l_zmat_f.rows() == 0) {
            throw Exception("l_zmat_f not allocated");
        }

        if(l_fscale_f.rows() == 0) {
            throw Exception("l_fscale_f not allocated");
        }

        if(jac_atm_f.rows() == 0) {
            throw Exception("jac_atm_f not allocated");
        }

        if(jac_surf_f.rows() == 0) {
            throw Exception("jac_surf_f not allocated");
        }
     }

}

void LRadDriver::calculate_first_order()
{
    check_rt_inputs();

    if(stokes_val_f.rows() != number_stokes()) {
        stokes_val_f.reference(Array<double, 1>(number_stokes(), ColumnMajorArray<1>()));
    }

    // Dump r1 value directly into stokes_val array
    Array<double, 1> r1(stokes_val_f);

    int nlayer = tau_f.rows();
    int nspars = surface_param_f.rows();
    int natm_jac = jac_atm_f.depth();

    l_rad_first_driver(&l_rad_struct, &nlayer, &natm_jac, &nstokes,
                       tau_f.dataFirst(), l_tau_f.dataFirst(),
                       omega_f.dataFirst(),
                       l_omega_f.dataFirst(), zmat_f.dataFirst(),
                       l_zmat_f.dataFirst(), fscale_f.dataFirst(),
                       l_fscale_f.dataFirst(), surface_param_f.dataFirst(),
                       &nspars,
                       &need_jacobians_i,
                       r1.dataFirst(), jac_atm_f.dataFirst(), jac_surf_f.dataFirst());
}

void LRadDriver::calculate_second_order()
{
    check_rt_inputs();

    // Only require pf_f to not be size of 0 for second order calculations
    if(pf_f.rows() == 0) {
       throw Exception("pf_f not allocated");
    }

    if (need_jacobians_i == 1) {
        if(l_pf_f.rows() == 0) {
            throw Exception("l_pf_f not allocated");
        }
    }

    if(stokes_val_f.rows() != number_stokes()) {
        stokes_val_f.reference(Array<double, 1>(number_stokes(), ColumnMajorArray<1>()));
    }

    int nscatt = pf_f.depth();
    int nlayer = tau_f.rows();
    int nmom = 2 * nstream;
    int nspars = surface_param_f.rows();
    int natm_jac = jac_atm_f.depth();

    Range ra(Range::all());
    Range rpol(1, number_stokes() - 1);

    // We only supply nstream coefficients to l_rad_second, so that is
    // essentially the number of effective coefficients
    Array<int, 1> neff_coefs(nlayer, ColumnMajorArray<1>());
    neff_coefs = nmom - 1;

    Array<double, 1> r2(number_stokes(), ColumnMajorArray<1>());
    Array<double, 2> ls_r2(number_stokes(), nspars, ColumnMajorArray<2>());
    Array<double, 3> l_r2(number_stokes(), nlayer,
                          natm_jac, ColumnMajorArray<3>());
    double icorr;
    Array<double, 1> ls_icorr(nspars);
    Array<double, 2> l_icorr(nlayer, natm_jac,
                             ColumnMajorArray<2>());

    // Use only the parts of the phase function expected by l_rad_second
    // for coeffs variables. We need to use a to_fortran command here
    // to ensure we have an array with column first ordering that is also 
    // contiguous. Using a reference to pf will not work.
    Array<double, 3> coefs_f(to_fortran(pf_f(Range(0, nmom - 1), ra, ra)));
    Array<double, 4> l_coefs_f(to_fortran(l_pf_f(Range(0, nmom - 1), ra, ra)));

    l_rad_second_driver(&l_rad_struct, &nlayer, &natm_jac, &nstokes,
                        &nmom, &nscatt,
                        tau_f.dataFirst(), l_tau_f.dataFirst(),
                        omega_f.dataFirst(), l_omega_f.dataFirst(),
                        surface_param_f.dataFirst(), &nspars,
                        coefs_f.dataFirst(), l_coefs_f.dataFirst(),
                        neff_coefs.dataFirst(),
                        &need_jacobians_i,
                        r2.dataFirst(), l_r2.dataFirst(), ls_r2.dataFirst(),
                        &icorr, l_icorr.dataFirst(), ls_icorr.dataFirst());

    // Second order calculation is split into parts, combine into how
    // they influence the stokes array and jacobians
    stokes_val_f(0) += icorr;
    stokes_val_f(rpol) += r2(rpol);

    if (need_jacobians_i == 1) {
        jac_surf_f(0, ra) += ls_icorr;
        jac_surf_f(rpol, ra) += ls_r2(rpol, ra);
        jac_atm_f(0, ra, ra) += l_icorr;
        jac_atm_f(rpol, ra, ra) += l_r2(rpol, ra, ra);  
    }
}

Array<double, 1> LRadDriver::stokes() const
{
    if (stokes_val_f.rows() == 0) {
        throw Exception("Stokes have not yet been allocated, use an RT operation first");
    }

    return stokes_val_f;
}

Array<double, 3> LRadDriver::atmospheric_jacobian() const
{
    if (jac_atm_f.rows() == 0) {
        throw Exception("Atmospheric jacobians array has not yet been allocated. Set up linear inputs first.");
    }

    return jac_atm_f;
}

Array<double, 2> LRadDriver::surface_jacobian() const
{
    if (jac_surf_f.rows() == 0) {
        throw Exception("Surface jacobians array has not yet been allocated. Set up linear inputs first.");
    }

    // A hack to deal with the different ordering.  There should be a better
    // place to do this.
    if (surface_type != 2)
        return jac_surf_f;
    else {
        Array<double, 2> temp(3,5);

        temp(Range::all(),0) = jac_surf_f(Range::all(),4);
        temp(Range::all(),1) = jac_surf_f(Range::all(),0);
        temp(Range::all(),2) = jac_surf_f(Range::all(),1);
        temp(Range::all(),3) = jac_surf_f(Range::all(),2);
        temp(Range::all(),4) = jac_surf_f(Range::all(),3);

        return temp;
    }
}

void LRadDriver::print(std::ostream& Os, bool Short_form) const
{
    Os << "LRadDriver:\n";
    Os << "  Number stream: " << nstream << "\n"
       << "  Number stokes: " << nstokes << "\n"
       << "  Surface type: " << surface_type << "\n"
       << "  Use TMS Correction: "
       << (use_tms_correction ? "True\n" : "False\n")
       << "  Pure nadir mode: "
       << (pure_nadir ? "True" : "False") << "\n";
}
