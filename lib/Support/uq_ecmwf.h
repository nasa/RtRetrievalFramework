#ifndef UQ_ECMWF_H
#define UQ_ECMWF_H
#include "meteorology.h"

namespace FullPhysics {

/****************************************************************//**
  This class implements the Uncertainty Quantification specific
  ECMWF reading functionality.
*******************************************************************/

class UqEcmwf : public Meteorology {
public:
    UqEcmwf(const std::string& Fname);
    ~UqEcmwf() {}

    // Define how to read various items
    using Meteorology::pressure_levels;
    blitz::Array<double, 1> pressure_levels() const
        { return read_array("ECMWF/vector_pressure_levels_ecmwf"); }

    using Meteorology::specific_humidity;
    blitz::Array<double, 1> specific_humidity() const
        { return read_array("ECMWF/specific_humidity_profile_ecmwf"); }

    double surface_pressure() const
        { return read_scalar("ECMWF/surface_pressure_ecmwf"); }

    double windspeed_u() const
        { return read_scalar("ECMWF/windspeed_u_ecmwf"); }

    double windspeed_v() const
        { return read_scalar("ECMWF/windspeed_v_ecmwf"); }

    using Meteorology::temperature;
    blitz::Array<double, 1> temperature() const
        { return read_array("ECMWF/temperature_profile_ecmwf"); }

    void print(std::ostream& Os) const
        { Os << "UqEcmwf"; }

private:

    //-----------------------------------------------------------------------
    /// UQ specific ECMWF reader routines
    //-----------------------------------------------------------------------

    double read_scalar(const std::string& Field) const;
    blitz::Array<double, 1> read_array(const std::string& Field) const;

    HdfFile h;
};
}
#endif
