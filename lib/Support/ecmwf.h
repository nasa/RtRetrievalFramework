#ifndef ECMWF_H
#define ECMWF_H
#include "meteorology.h"
#include "hdf_sounding_id.h"
#include "hdf_file.h"
#include "printable.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to read some of the fields from the ECMWF file,
  which can then be used for things such as the apriori.

  Since resampled ECMWF files can differ between instrument types,
  the read routines are pure virtual and need to be implemented 
  for the specifics of the instrument specific ECMWF files.
*******************************************************************/

class Ecmwf : public Meteorology {
public:
    virtual ~Ecmwf() {}
 
    //-----------------------------------------------------------------------
    /// Return the H20 VMR. This is the specific_humidity converted to a
    /// volume mixing ratio.
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> h2o_vmr() const;

    //-----------------------------------------------------------------------
    /// Ozone mass mixing ratio
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> ozone_mmr() const = 0;

    //-----------------------------------------------------------------------
    /// Return the Ozone VMR.
    //-----------------------------------------------------------------------
    virtual blitz::Array<double, 1> ozone_vmr() const;

    //-----------------------------------------------------------------------
    /// Return either the H2O or O3 (Ozone) vmr, error otherwise.
    //-----------------------------------------------------------------------
    using Meteorology::vmr;
    virtual blitz::Array<double, 1> vmr(const std::string& Species) const;

    void print(std::ostream& Os) const { Os << "Ecmwf"; }

};
}
#endif
