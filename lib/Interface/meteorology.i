%include "common.i"

%{
#include "meteorology.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::Meteorology);
%nodefaultctor FullPhysics::Meteorology;

namespace FullPhysics {
class Meteorology : public GenericObject {
public:
    virtual ~Meteorology();
    std::string print_to_string() const;
    %python_attribute(pressure_levels, virtual blitz::Array<double, 1>);
    %python_attribute(specific_humidity, virtual blitz::Array<double, 1>);
    virtual blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const;
    blitz::Array<double, 1> vmr(const std::string& Species) const;
    blitz::Array<double, 1> vmr(const std::string& Species, const blitz::Array<double, 1>& Pressure_level) const;
    %python_attribute(temperature, virtual blitz::Array<double, 1>);
    virtual blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const;
    virtual double surface_pressure() const;
    virtual double windspeed() const;
};
}

