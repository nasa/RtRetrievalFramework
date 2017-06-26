#include "oco_forward_model_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> oco_fm_out_create
(const boost::shared_ptr<ForwardModel>& fm)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new OcoForwardModelOutput
     (boost::dynamic_pointer_cast<OcoForwardModel>(fm)));
}
REGISTER_LUA_DERIVED_CLASS(OcoForwardModelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &oco_fm_out_create)
]
REGISTER_LUA_END()
#endif

// Adapter class to take things from OcoForwardModel into format needed
// for output
class OcoForwardModelCalc {
public:
  OcoForwardModelCalc(const boost::shared_ptr<OcoForwardModel>& Fm)
    : fm(Fm) {}

  Array<int, 1> samples() const 
  {
    int num_total_samples = fm->measured_radiance_all().spectral_domain().
      wavenumber().rows();
    Array<int, 1> samples(num_total_samples);
    for(int i = 0; i < fm->number_spectrometer(); ++i) {
      boost::optional<blitz::Range> pr = fm->pixel_range(i);
      if(pr)
	samples(*pr) = fm->spectral_grid()->low_resolution_grid(i).sample_index();
    }
    return samples;
  }

private:
    boost::shared_ptr<OcoForwardModel> fm;
};

// See base class for description

void OcoForwardModelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    boost::shared_ptr<OcoForwardModelCalc> fmc(new OcoForwardModelCalc(fm));
    out->register_data_source("/SpectralParameters/sample_indexes",
                              &OcoForwardModelCalc::samples, fmc);
}
