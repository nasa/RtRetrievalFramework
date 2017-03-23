#ifndef AEROSOL_CONSOLIDATED_OUTPUT_H
#define AEROSOL_CONSOLIDATED_OUTPUT_H
#include "register_output_base.h"
#include "aerosol_optical.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Aerosol class that should be
  written as output.

  Creates consolidated matrices with fill values for aerosol types
  not used.

  Creates:
  * aerosol_gaussian_log_param - num_all_aer * 3
  * aerosol_gaussian_log_param_apriori - num_all_aer * 3
  * aerosol_gaussian_log_param_uncert - num_all_aer * 3
  * aerosol_aod - num_all_aer * num_aod
  
  Where num_aod = 4 : [total, low, mid, high]


  See the discussion in RegisterOutputBase why this isn't just part of
  the Aerosol class.
*******************************************************************/
class AerosolConsolidatedOutput : public RegisterOutputBase {
public:

    //-----------------------------------------------------------------------
    /// Constructor. 
    ///
    /// All_aer_names defines the aerosol type index of the output matrix
    //-----------------------------------------------------------------------
    AerosolConsolidatedOutput(const boost::shared_ptr<Aerosol>& Aerosol, const std::vector<std::string>& All_aer_names);
    virtual ~AerosolConsolidatedOutput() {}
    virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
    static const double low_boundary;
    static const double high_boundary;

private:
    blitz::Array<double, 2>& aerosol_aod_matrix();
    boost::shared_ptr<AerosolOptical> aerosol;
    std::vector<std::string> all_aer_names;
};
}
#endif
