#ifndef ABSORBER_VMR_SHAPE_H
#define ABSORBER_VMR_SHAPE_H
#include "absorber_vmr_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the absorbver VMR portion of the state. The
  absorber is retrieved as a set of scalings applied to a set of
  shape profiles combined with an initial guess VMR profile.

  The shape is computed as follows for each level i
  absorber_final(i) = absorber_base(i) + scaling1 * shape1(i) + scaling2 * shape2(i) + ... scalingN * shapeN(i)

  Where scalingN are the scale factors and shapeN are the shape
  profile values.
*******************************************************************/
class AbsorberVmrShape : public AbsorberVmrImpBase {
public:
    AbsorberVmrShape(const blitz::Array<double, 1> VMR_base,
                     const blitz::Array<double, 2> Shape_profile,
                     const blitz::Array<double, 1> Shape_scaling,
                     const boost::shared_ptr<Pressure>& Press,
                     const blitz::Array<bool, 1>& Shape_flag,
                     const std::string& Gas_name,
                     const bool Log_profiles = false);

    virtual ~AbsorberVmrShape() = default;

    virtual void print(std::ostream& Os) const;

    virtual boost::shared_ptr<AbsorberVmr> clone() const
    {
        return clone(press);
    }

    virtual boost::shared_ptr<AbsorberVmr> clone(const boost::shared_ptr<Pressure>& Press) const;

    virtual std::string sub_state_identifier() const
    {
        return "absorber_vmr_shape";
    }

    virtual std::string state_vector_name_i(int i) const;

    //-----------------------------------------------------------------------
    /// Base VMR profile associated with the pressure profile, values are in 
    /// Kelvin
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> vmr_base() const
    { return vmr_base_; }

    //-----------------------------------------------------------------------
    /// Shape profiles combined with base VMR dimenison are:
    /// N_level x N_scaling
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 2> shape_profile() const
    { return shape_prof; }

    //-----------------------------------------------------------------------
    /// Scale factors combined with the shape profiles with dimension:
    /// N_scaling 
    //-----------------------------------------------------------------------

    virtual blitz::Array<double, 1> shape_scaling() const
    { return coeff.value(); }

    //-----------------------------------------------------------------------
    /// Pressure levels that serve as the grid for the VMR values in
    /// units of Pascals
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> pressure_profile() const
    { return press->pressure_grid().value.value(); }

protected:
    virtual void calc_vmr() const;
private:
    blitz::Array<double, 1> vmr_base_;
    blitz::Array<double, 2> shape_prof;
    bool log_profiles;
};
}
#endif
