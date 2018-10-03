#ifndef ACOS_MET_FILE_H
#define ACOS_MET_FILE_H
#include "meteorology.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements the ACOS specific ECMWF reading 
  functionality.
*******************************************************************/

class AcosMetFile : public Meteorology {
public:
    AcosMetFile(const std::string& Fname, const boost::shared_ptr<HdfSoundingId>& 
    Hdf_sounding_id, bool Avg_sounding_number);
    AcosMetFile(const std::string& Fname, const HeritageFile& Run_file);
    ~AcosMetFile() {}

    // Define how to read various items
    using Meteorology::pressure_levels;
    blitz::Array<double, 1> pressure_levels() const
        { return read_array("/ecmwf/specific_humidity_pressures"); }

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

    void print(std::ostream& Os) const { Os << "AcosMetFile"; }

private:

    //-----------------------------------------------------------------------
    /// ACOS specific ECMWF reader routines
    //-----------------------------------------------------------------------

    double read_scalar(const std::string& Field) const;
    blitz::Array<double, 1> read_array(const std::string& Field) const;

    HdfFile h;
    boost::shared_ptr<HdfSoundingId> hsid;
    bool average_sounding_number;
};
}
#endif
