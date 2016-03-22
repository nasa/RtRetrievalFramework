%include "common.i"
%{
#include "gas_vmr_apriori.h"
#include "temperature.h"
#include "altitude.h"
%}

%base_import(generic_object)
%import "ecmwf.i"
%import "level_1b.i"
%import "pressure.i"
%import "temperature.i"
%import "altitude.i"
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::GasVmrApriori);

namespace FullPhysics {

class GasVmrApriori : public GenericObject {
public:
    GasVmrApriori(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
                  const boost::shared_ptr<Level1b>& L1b_file,
                  const boost::shared_ptr<Altitude>& Alt,
                  const HdfFile& Hdf_static_input,
                  const std::string& Hdf_group,
                  const std::string& Gas_name);

    %python_attribute(apriori_vmr, blitz::Array<double, 1>)

    std::string print_to_string() const;
};
}
