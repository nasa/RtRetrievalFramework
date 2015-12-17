// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "l_rad_driver.h"
%}

%import "array_ad.i"
%fp_shared_ptr(FullPhysics::LRadDriver);

namespace FullPhysics {
class LRadDriver {
public:
    enum PsMode {REGULAR, ENHANCED, PLANE_PARALLEL, DETECT};

    LRadDriver(int Number_stream, int Number_stokes,
               int surface_type,
               bool Tms_Correction = false,
               bool Pure_nadir = false,
               const PsMode ps_mode = DETECT);
  
    %python_attribute(number_stokes, virtual int)
    %python_attribute(number_stream, virtual int)

    ArrayAd<double, 2> z_matrix(const ArrayAd<double, 3>& pf) const;

    virtual void setup_geometry(blitz::Array<double, 1> alt, double sza, double zen, double azm) const;

    virtual void setup_surface_params(const blitz::Array<double, 1>& surface_param);

    virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                      const blitz::Array<double, 1>& ssa,
                                      const blitz::Array<double, 3>& pf,
                                      const blitz::Array<double, 2>& zmat);

    virtual void clear_linear_inputs();

    virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                     const ArrayAd<double, 1>& ssa,
                                     const ArrayAd<double, 3>& pf,
                                     const ArrayAd<double, 2>& zmat);

    virtual void calculate_first_order();
    virtual void calculate_second_order();

    virtual blitz::Array<double, 1> stokes() const;
    virtual blitz::Array<double, 3> atmospheric_jacobian() const;
    virtual blitz::Array<double, 2> surface_jacobian() const;

    virtual void print(std::ostream& Os, bool Short_form = false) const;

};
}
