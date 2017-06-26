#ifndef REFERENCE_VMR_APRIORI_H
#define REFERENCE_VMR_APRIORI_H

#include <blitz/array.h>
#include "fp_time.h"
#include "double_with_unit.h"

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

  NOTE: Inputs are expected to be in increasing altitude
  decreasing pressure order.

  Make sure gas names are capatilized.
*******************************************************************/

class ReferenceVmrApriori : public Printable<ReferenceVmrApriori> {
public:

    ReferenceVmrApriori(const blitz::Array<double, 1>& Model_pressure,
                        const blitz::Array<double, 1>& Model_altitude,
                        const blitz::Array<double, 1>& Model_temperature,
                        const blitz::Array<double, 1>& Ref_altitude,
                        const double Ref_latitude,
                        const Time& Ref_time,
                        const double Ref_tropopause_altitude,
                        const double Obs_latitude,
                        const Time& Obs_time);

    DoubleWithUnit model_tropopause_altitude() const;
    const blitz::Array<double, 1> effective_altitude() const;
    const double age_of_air(const double altitude) const;

    const blitz::Array<double, 1> resample_to_model_grid(const blitz::Array<double, 1>& vmr) const;
    const blitz::Array<double, 1> apply_latitude_gradient(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const;
    const blitz::Array<double, 1> apply_secular_trend(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const;
    const blitz::Array<double, 1> apply_seasonal_cycle(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const;

    const blitz::Array<double, 1> apriori_vmr(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const;

    void print(std::ostream& Os) const { Os << "ReferenceVmrApriori"; }

private:
    
    blitz::Array<double, 1> model_pressure;
    blitz::Array<double, 1> model_altitude;
    blitz::Array<double, 1> model_temperature;

    blitz::Array<double, 1> ref_altitude;
    double ref_latitude;
    Time ref_time;
    double ref_tropopause_altitude;

    double obs_latitude;
    Time obs_time;

};
}

#endif

