#ifndef L_RAD_DRIVER_H
#define L_RAD_DRIVER_H

#include <boost/noncopyable.hpp>
#include <blitz/array.h>

#include "linear_interpolate.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class drives the LRAD code, which gives a polarization
  correction to scalar intensity and jacobians.

  The correction used here is described in the paper "A fast
  linearized pseudo-spherical two orders of scattering model to
  account for polarization in vertically inhomogeneous
  scattering-absorbing media" by Vijah Natrah and Robert Spurr in
  Journal of Quantitative Spectroscopy & Radiative Transfer 107 (2007)
  263-293 (a copy of this paper can be found in the source tree at
  'doc/LRAD Paper.pdf').
*******************************************************************/

class LRadDriver : public Printable<LRadDriver>, public boost::noncopyable {

public:
    enum PsMode {REGULAR, ENHANCED, PLANE_PARALLEL, DETECT};

    LRadDriver(int Number_stream, int Number_stokes,
               int surface_type,
               bool Tms_Correction = false,
               bool Pure_nadir = false,
               const PsMode ps_mode = DETECT);

    virtual ~LRadDriver();

    virtual int number_stokes() const {return nstokes;}
    virtual int number_stream() const {return nstream;}

    ArrayAd<double, 2> z_matrix(const ArrayAd<double, 3>& pf) const;

    /// Setup viewing geometry, should only be called once per instance or if
    /// the viewing geometry changes
    virtual void setup_geometry(blitz::Array<double, 1> alt, double sza, double zen, double azm) const;

    /// Set up surface parameters for spectral point
    virtual void setup_surface_params(const blitz::Array<double, 1>& surface_param);

    /// Set up optical depth, single scattering albedo and scattering matrix
    /// Should be called per spectral point
    virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                      const blitz::Array<double, 1>& ssa,
                                      const blitz::Array<double, 3>& pf,
                                      const blitz::Array<double, 2>& zmat);

    /// Mark that we are not retrieving weighting functions
    virtual void clear_linear_inputs();

    /// Set up linearization, weighting functions
    virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                     const ArrayAd<double, 1>& ssa,
                                     const ArrayAd<double, 3>& pf,
                                     const ArrayAd<double, 2>& zmat);

    /// Perform radiative transfer calculation with the values
    /// setup by setup_optical_inputs and setup_linear_inputs
    /// for the first order of scattering
    virtual void calculate_first_order();

    /// Perform radiative transfer calculation with the values
    /// setup by setup_optical_inputs and setup_linear_inputs
    /// for the second order of scattering
    virtual void calculate_second_order();

    /// Retrieve the stokes values calculated
    virtual blitz::Array<double, 1> stokes() const;

    /// Atmospheric jacobian from last calculation
    virtual blitz::Array<double, 3> atmospheric_jacobian() const;
    
    /// Surface jacobian 
    virtual blitz::Array<double, 2> surface_jacobian() const;

    virtual void print(std::ostream& Os, bool Short_form = false) const;

private:

    /// Set the pseudo-spherical mode used by l_rad. Normally only used to force
    /// a mode for testing
    void set_ps_mode(const PsMode& mode)
    {
        switch(mode) {
        case REGULAR:
            regular_ps = true;
            enhanced_ps = false;
            break;

        case ENHANCED:
            regular_ps = false;
            enhanced_ps = true;
            break;

        case PLANE_PARALLEL:
            regular_ps = false;
            enhanced_ps = false;
            break;

        case DETECT:
        default:

            // Determine which pseudo-spherical mode to use based on the zenith angle
            // ie, how close are we to nadir
            if (pure_nadir) {
                set_ps_mode(REGULAR);
            } else {
                set_ps_mode(ENHANCED);
            }

            break;
        }
    }

    void initialize(const PsMode ps_mode);
    void check_rt_inputs();

    // Options configured at construction
    int nstream, nstokes;
    int surface_type;
    bool use_tms_correction;
    bool pure_nadir;
    bool regular_ps, enhanced_ps;

    // Integer boolean for l_rad indicating if jacobians are needed, false until linear
    // values are initialized
    int need_jacobians_i;

    // Fortran ordered Blitz arrays used when calling l_rad. They are initialized on first use
    blitz::Array<double, 1> surface_param_f; 
    blitz::Array<double, 1> tau_f; 
    blitz::Array<double, 1> omega_f; 
    blitz::Array<double, 3> pf_f; 
    blitz::Array<double, 2> zmat_f; 
    blitz::Array<double, 1> fscale_f; 
    blitz::Array<double, 3> jac_atm_f; 
    blitz::Array<double, 2> jac_surf_f; 
    blitz::Array<double, 2> l_tau_f; 
    blitz::Array<double, 2> l_omega_f; 
    blitz::Array<double, 4> l_pf_f; 
    blitz::Array<double, 3> l_zmat_f; 
    blitz::Array<double, 2> l_fscale_f; 
    blitz::Array<double, 1> stokes_val_f; 

    /// This points to Fortran 90 structure. We have a copy up here in
    /// C++ because we maintain the lifetime of this object here.
    mutable void *l_rad_struct;

};
}
#endif
