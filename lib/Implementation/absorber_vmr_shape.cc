#include "absorber_vmr_shape.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrShape, AbsorberVmr)
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 1>&,
                          const boost::shared_ptr<Pressure>&,
                          const blitz::Array<bool, 1>&,
                          const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrShape::AbsorberVmrShape(const blitz::Array<double, 1> VMR_base,
                                   const blitz::Array<double, 2> Shape_profile,
                                   const blitz::Array<double, 1> Shape_scaling,
                                   const boost::shared_ptr<Pressure>& Press,
                                   const blitz::Array<bool, 1>& Shape_flag,
                                   const std::string& Gas_name)

: AbsorberVmrImpBase(Gas_name, Shape_scaling, Shape_flag, Press, false),
  vmr_base_(VMR_base), shape_prof(Shape_profile)
{
}

//-----------------------------------------------------------------------
/// This calculates VMR grid to use for layer retrieval.
//-----------------------------------------------------------------------

void AbsorberVmrShape::calc_vmr() const
{
    blitz::Array<double, 1> press_profile( pressure_profile() );

    if (press_profile.rows() != vmr_base_.rows()) {
        std::stringstream err_msg;
        err_msg << "Size of pressure grid: "
                << press_profile.rows()
                << " != size of VMR levels: "
                << vmr_base_.rows();
        throw Exception(err_msg.str());
    }

    if (press_profile.rows() != shape_prof.rows()) {
        std::stringstream err_msg;
        err_msg << "Size of pressure grid: "
                << press_profile.rows()
                << " != size of shape profile levels: "
                << shape_prof.rows();
        throw Exception(err_msg.str());
    }

    std::vector<AutoDerivative<double> > plist;
    std::vector<AutoDerivative<double> > vmrlist;

    for(int lev_idx = 0; lev_idx < press_profile.rows(); ++lev_idx) {
        plist.push_back(press_profile(lev_idx));

        AutoDerivative<double> vmr_val = vmr_base_(lev_idx);
        for(int shape_idx = 0; shape_idx < shape_prof.cols(); shape_idx++) {
            vmr_val += coeff(shape_idx) * shape_prof(lev_idx, shape_idx);
        }
        vmrlist.push_back(vmr_val);
    }

    typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > lin_type;
    boost::shared_ptr<lin_type> lin(new lin_type(plist.begin(), plist.end(), vmrlist.begin()));
    vmr = boost::bind(&lin_type::operator(), lin, _1);
}

// See base class for description of this function.
std::string AbsorberVmrShape::state_vector_name_i(int i) const
{
    return gas_name() + " VMR Shape Scaling #" +
           boost::lexical_cast<std::string>(i + 1);
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrShape::clone
(const boost::shared_ptr<Pressure>& Press) const
{
    return boost::shared_ptr<AbsorberVmr>
        (new AbsorberVmrShape(vmr_base_, shape_prof, coeff.value(), Press->clone(), used_flag, gas_name()));
}

void AbsorberVmrShape::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "AbsorberVmrShape:\n"
       << "  Gas name:       " << gas_name() << "\n";
    opad.strict_sync();
}
