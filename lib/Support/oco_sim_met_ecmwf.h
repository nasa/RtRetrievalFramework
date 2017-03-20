#ifndef OCO_SIM_MET_ECMWF_H
#define OCO_SIM_MET_ECMWF_H
#include "meteorology.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the OCO specific ECMWF reading 
  functionality.

  This reads the simulator meteorology files. This are similar to 
  the OCO ECMWF, but the fields have different names.

  Note that the actual OCO simulator used "scene" files, which are
  somewhat lime the meteorology but in a different format, and with
  a different number of levels (not the normal 91 ECMWF). The
  meteorology files are this scene information resampled to the 91
  levels. 
*******************************************************************/

class OcoSimMetEcmwf : public Meteorology {
public:
    OcoSimMetEcmwf(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& Hdf_sounding_id);
    ~OcoSimMetEcmwf() {}

    // Define how to read various items
    using Meteorology::pressure_levels;
    blitz::Array<double, 1> pressure_levels() const
        { return read_array("ecmwf/specific_humidity_pressures"); }

    using Meteorology::specific_humidity;
    blitz::Array<double, 1> specific_humidity() const
        { return read_array("ecmwf/specific_humidity"); }

    double surface_pressure() const
        { return read_scalar("ecmwf/surface_pressure"); }

    double windspeed_u() const
        { return read_scalar("ecmwf/windspeed_u"); }

    double windspeed_v() const
        { return read_scalar("ecmwf/windspeed_v"); }

    using Meteorology::temperature;
    blitz::Array<double, 1> temperature() const
        { return read_array("ecmwf/temperature"); }

    void print(std::ostream& Os) const { Os << "OcoSimMetEcmwf"; }

private:

    //-----------------------------------------------------------------------
    /// OCO simulator meteorology specific ECMWF reader routines
    //-----------------------------------------------------------------------

    double read_scalar(const std::string& Field) const;
    blitz::Array<double, 1> read_array(const std::string& Field) const;

    HdfFile h;
    boost::shared_ptr<HdfSoundingId> hsid;
    bool average_sounding_number;
};
}
#endif
