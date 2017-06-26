%include "common.i"
%{
#include "reference_vmr_apriori.h"
%}

%base_import(generic_object)
%import "fp_time.i"
%import "double_with_unit.i"
%fp_shared_ptr(FullPhysics::ReferenceVmrApriori);

namespace FullPhysics {
class ReferenceVmrApriori : public GenericObject {
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
    const blitz::Array<double, 1> apply_latitude_gradient(const blitz::Array<double, 1>& vmr, std::string& gas_name) const;
    const blitz::Array<double, 1> apply_secular_trend(const blitz::Array<double, 1>& vmr, std::string& gas_name) const;
    const blitz::Array<double, 1> apply_seasonal_cycle(const blitz::Array<double, 1>& vmr, std::string& gas_name) const;

    const blitz::Array<double, 1> apriori_vmr(const blitz::Array<double, 1>& vmr, std::string& gas_name) const;

    std::string print_to_string() const;
};
}
