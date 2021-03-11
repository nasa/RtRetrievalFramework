%include "fp_common.i"
%{
#include "temperature_level.h"
%}
%base_import(temperature_imp_base)

%fp_shared_ptr(FullPhysics::TemperatureLevel);

namespace FullPhysics {
class TemperatureLevel: public TemperatureImpBase {
public:
    TemperatureLevel(const blitz::Array<double, 1> Temp,
                     const boost::shared_ptr<Pressure>& Press,
                     const blitz::Array<bool, 1>& Temp_flag);

    virtual ~TemperatureLevel();
    virtual void print(std::ostream& Os) const;
    virtual boost::shared_ptr<Temperature> clone(const boost::shared_ptr<Pressure>& Press) const;
    virtual std::string sub_state_identifier() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual blitz::Array<double, 1> temperature_profile() const;
    virtual blitz::Array<double, 1> pressure_profile() const;
    virtual ArrayWithUnit<double, 1> important_pressure_level() const;
protected:
    void calc_temperature_grid() const;
};
}
