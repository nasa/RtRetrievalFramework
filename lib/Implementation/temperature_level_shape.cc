#include "temperature_level_shape.h"

#include "linear_interpolate.h"

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureLevelShape, Temperature)
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 1>&,
                          const boost::shared_ptr<Pressure>&,
                          const blitz::Array<bool, 1>&>())
REGISTER_LUA_END()
#endif


//-----------------------------------------------------------------------
/// Set up Temperature
//-----------------------------------------------------------------------

TemperatureLevelShape::TemperatureLevelShape(const blitz::Array<double, 1> Temp_base,
                                             const blitz::Array<double, 2> Shape_profile,
                                             const blitz::Array<double, 1> Shape_scaling,
                                             const boost::shared_ptr<Pressure>& Press,
                                             const blitz::Array<bool, 1>& Shape_flag)
: pressure(Press), temp_base(Temp_base), shape_prof(Shape_profile) 
{
    init(Shape_scaling, Shape_flag, Press);
}

//-----------------------------------------------------------------------
/// This calculates temperature grid to use for layer retrieval.
//-----------------------------------------------------------------------

void TemperatureLevelShape::calc_temperature_grid() const
{
    blitz::Array<double, 1> press_profile( pressure_profile() );

    if (press_profile.rows() != temp_base.rows()) {
        std::stringstream err_msg;
        err_msg << "Size of pressure grid: "
                << press_profile.rows()
                << " != size of temperature levels: "
                << temp_base.rows();
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
    std::vector<AutoDerivative<double> > tlist;

    for(int lev_idx = 0; lev_idx < press_profile.rows(); ++lev_idx) {
        plist.push_back(press_profile(lev_idx));

        AutoDerivative<double> temp_val = temp_base(lev_idx);
        for(int shape_idx = 0; shape_idx < shape_prof.cols(); shape_idx++) {
            temp_val += coeff(shape_idx) * shape_prof(lev_idx, shape_idx);
        }
        tlist.push_back(temp_val);
    }

    typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > lin_type;
    boost::shared_ptr<lin_type> lin(new lin_type(plist.begin(), plist.end(), tlist.begin()));
    tgrid = boost::bind(&lin_type::operator(), lin, _1);
}

// See base class for description of this function.
std::string TemperatureLevelShape::state_vector_name_i(int i) const
{
    return "Temperature Shape Scaling #" +
           boost::lexical_cast<std::string>(i + 1);
}

boost::shared_ptr<Temperature> TemperatureLevelShape::clone(const boost::shared_ptr<Pressure>& Press) const
{
    return boost::shared_ptr<TemperatureLevelShape>
        (new TemperatureLevelShape(temp_base, shape_prof, coeff.value(), Press->clone(), used_flag));
}

void TemperatureLevelShape::print(std::ostream& Os) const
{
    Os << "Temperature Level Shape\n";
}
