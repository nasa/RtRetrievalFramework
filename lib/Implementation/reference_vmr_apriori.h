#ifndef REFERENCE_VMR_APRIORI_H
#define REFERENCE_VMR_APRIORI_H
#include "level_1b.h"

namespace FullPhysics {

/****************************************************************//**
  Creates a VMR profile for a gas using a set of dated reference
  VMRs with a known latitude. These VMRs are then modified as so:
  1. Resampled to effective altitudes
  2. Latitude gradient applied
  3. Secular trends applied
  4. Season cycle applied

  This class is based on the TCCON 2014 release of gsetup. As
  per those techniques, values are interpolated based on altitudes.

  Inputs are expected to be inin increasing altitude
  decreasing pressure order.
*******************************************************************/

class ReferenceVmrApriori : public Printable<ReferenceVmrApriori> {
public:

    ReferenceVmrApriori(const blitz::Array<double, 1>& Model_altitude,
                        const blitz::Array<double, 1>& Model_temperature,
                        const double Model_latitude,
                        const blitz::Array<double, 1>& Ref_altitude,
                        const double Ref_latitude,
                        const double Ref_tropopause_altitude);

    double model_tropopause_altitude() const;
    const blitz::Array<double, 1> effective_altitude() const;
    const blitz::Array<double, 1> resample_to_model_grid(const blitz::Array<double, 1>& vmr) const;

    void print(std::ostream& Os) const { Os << "ReferenceVmrApriori"; }

private:
    
    void resample_at_effective_altitudes();
    void apply_latitude_gradients();
    void apply_secular_trends();
    void apply_seasonal_cycle();

    blitz::Array<double, 1> model_altitude;
    blitz::Array<double, 1> model_temperature;
    double model_latitude;

    blitz::Array<double, 1> ref_altitude;
    double ref_latitude;
    double ref_tropopause_altitude;

};
}

#endif

