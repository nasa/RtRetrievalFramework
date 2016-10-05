#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "output_hdf.h"
#include "high_res_spectrum_output.h"

using namespace FullPhysics;
using namespace blitz;

class NamedSpectrumTestObservable : public Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> > {
  public:
    void send_spectrum(const Spectrum& spec, const std::string& name, int index) {
      notify_update_do(boost::shared_ptr<NamedSpectrum>(new NamedSpectrum(spec, name, index)));
    }

    virtual void add_observer(Observer<boost::shared_ptr<NamedSpectrum> > & Obs) 
    { add_observer_do(Obs); }
    virtual void remove_observer(Observer<boost::shared_ptr<NamedSpectrum> >& Obs) 
    { remove_observer_do(Obs); }

};

BOOST_FIXTURE_TEST_SUITE(high_res_spectrum_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<OutputHdf> out(new OutputHdf("high_res_spectrum_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("high_res_spectrum_output.h5");

  // Simple test, we just make sure that we can write output
  // Try writing some made up Spectrum classes with a made up name
  Array<double,1> sd1(20);
  sd1 = 1;
  Array<double,1> sd2(30);
  sd2 = 2;
  Array<double,1> sd3(40);
  sd3 = 3;
  Array<double,1> sr1(20);
  sr1 = 4;
  Array<double,1> sr2(30);
  sr2 = 5;
  Array<double,1> sr3(40);
  sr3 = 6;
  Spectrum spec1(SpectralDomain(sd1, units::inv_cm), SpectralRange(sr1, units::dimensionless));
  Spectrum spec2(SpectralDomain(sd2, units::inv_cm), SpectralRange(sr2, units::dimensionless));
  Spectrum spec3(SpectralDomain(sd3, units::inv_cm), SpectralRange(sr3, units::dimensionless));

  NamedSpectrumTestObservable spec_obs;
  HighResSpectrumOutput hi_res_out;

  spec_obs.add_observer(hi_res_out);
  hi_res_out.register_output(out);

  spec_obs.send_spectrum(spec1, "spectrum", 0);
  spec_obs.send_spectrum(spec2, "spectrum", 1);
  spec_obs.send_spectrum(spec3, "spectrum", 2);

  out->write();

  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


