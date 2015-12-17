#include "forward_model_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ForwardModelOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<ForwardModel>& >())
REGISTER_LUA_END()
#endif

// Adapter class to take things from ForwardModel into format needed
// for output
class ForwardModelCalc {
public:
    ForwardModelCalc(const boost::shared_ptr<ForwardModel>& Fm)
        : fm(Fm) {}
    blitz::Array<double, 1> spectral_domain_values() const
    {
        return fm->measured_radiance_all().spectral_domain().data();
    }

    std::string spectral_domain_name() const
    {
        return fm->spectral_domain_type_preference() == SpectralDomain::PREFER_WAVELENGTH ?
               "wavelength" : "wavenumber" ;
    }

    blitz::Array<double, 1> measured_radiance() const
    {
        return fm->measured_radiance_all().spectral_range().data();
    }

    blitz::Array<double, 1> measured_radiance_uncert() const
    {
        return fm->measured_radiance_all().spectral_range().uncertainty();
    }

    int number_pixel() const
    {
        return fm->measured_radiance_all().spectral_domain().wavenumber().rows();
    }

    Array<int, 1> number_pixel_each_band() const
    {
        Array<int, 1> res(fm->number_spectrometer());

        for(int i = 0; i < res.rows(); ++i) {
            boost::optional<blitz::Range> pr = fm->pixel_range(i);

            if(pr) {
                res(i) = pr->length();
            } else {
                res(i) = 0;
            }
        }

        return res;
    }

private:
    boost::shared_ptr<ForwardModel> fm;
};

// See base class for description

void ForwardModelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    boost::shared_ptr<ForwardModelCalc> fmc(new ForwardModelCalc(fm));
    out->register_data_source("/SpectralParameters/num_colors",
                              &ForwardModelCalc::number_pixel, fmc);
    out->register_data_source("/SpectralParameters/num_colors_per_band",
                              &ForwardModelCalc::number_pixel_each_band, fmc);
    out->register_data_source("/SpectralParameters/" + fmc->spectral_domain_name(),
                              &ForwardModelCalc::spectral_domain_values, fmc);
    out->register_data_source("/SpectralParameters/measured_radiance",
                              &ForwardModelCalc::measured_radiance, fmc);
    out->register_data_source("/SpectralParameters/measured_radiance_uncert",
                              &ForwardModelCalc::measured_radiance_uncert, fmc);
}
