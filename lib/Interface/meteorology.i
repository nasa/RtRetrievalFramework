%include "common.i"

%{
#include "meteorology.h"
%}

namespace FullPhysics {
class Meteorology : public Printable<Meteorology> {
public:
    virtual ~Meteorology();
    std::string print_to_string() const;
    %python_attribute(pressure_levels, virtual blitz::Array<double, 1>);
    %python_attribute(specific_humidity, virtual blitz::Array<double, 1>);
    virtual blitz::Array<double, 1> specific_humidity(const blitz::Array<double, 1>& Pressure_level) const = 0;
    blitz::Array<double, 1> vmr(const std::string& Species) const = 0;
    blitz::Array<double, 1> vmr(const std::string& Species, const blitz::Array<double, 1>& Pressure_level) const = 0;
    %python_attribute(temperature, virtual blitz::Array<double, 1>);
    virtual blitz::Array<double, 1> temperature(const blitz::Array<double, 1>& Pressure_level) const = 0;
    virtual double surface_pressure() const = 0;
    virtual double windspeed() const = 0;
};
}

