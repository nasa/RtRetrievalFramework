#ifndef OCO_ECMWF_H
#define OCO_ECMWF_H
#include "meteorology.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the OCO specific ECMWF reading 
  functionality.
*******************************************************************/

class OcoEcmwf : public Meteorology {
public:
    OcoEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
    ~OcoEcmwf() {}

    // Define how to read various items
    using Meteorology::pressure_levels;
    blitz::Array<double, 1> pressure_levels() const
        { return read_array("ECMWF/vector_pressure_levels_ecmwf"); }

    using Meteorology::specific_humidity;
    blitz::Array<double, 1> specific_humidity() const
        { return read_array("ECMWF/specific_humidity_profile_ecmwf"); }

    using Meteorology::vmr;
    blitz::Array<double, 1> vmr(const std::string& Species) const;

    blitz::Array<double, 1> ozone_mmr() const
        { return read_array("/ECMWF/ozone_profile_ecmwf"); }

    virtual blitz::Array<double, 1> ozone_vmr() const;

    double surface_pressure() const
        { return read_scalar("ECMWF/surface_pressure_ecmwf"); }

    double windspeed_u() const
        { return read_scalar("ECMWF/windspeed_u_ecmwf"); }

    double windspeed_v() const
        { return read_scalar("ECMWF/windspeed_v_ecmwf"); }

    using Meteorology::temperature;
    blitz::Array<double, 1> temperature() const
        { return read_array("ECMWF/temperature_profile_ecmwf"); }

    void print(std::ostream& Os) const { Os << "OcoEcmwf"; }

private:

    //-----------------------------------------------------------------------
    /// OCO specific ECMWF reader routines
    //-----------------------------------------------------------------------

    double read_scalar(const std::string& Field) const;
    blitz::Array<double, 1> read_array(const std::string& Field) const;

    HdfFile h;
    boost::shared_ptr<HdfSoundingId> hsid;
    bool average_sounding_number;
};
}
#endif
