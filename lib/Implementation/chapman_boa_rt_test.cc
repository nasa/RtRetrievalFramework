#include "chapman_boa_rt.h"
#include "unit_test_support.h"

#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;


SpectralDomain dummyd;
DoubleWithUnit dummyd2;
class ChapmanBoaRtFixture: public ConfigurationFtsFixture {
public:
  ChapmanBoaRtFixture() 
  {
    int spec_index = 0;
    spec_domain.reset(new SpectralDomain(config_spectrum_sampling->
 spectral_domain(spec_index, lowres_grid(spec_index), ils_half_width(spec_index))));
    // Spectral windows for better refraction calculation
    SpectralBound sb = config_spectral_window->spectral_bound();

    boost::shared_ptr<Level1b> l1b_fts = config_level_1b;

    // Create RT object to use
    Array<double, 1> sza(l1b_fts->number_spectrometer());
    for(int s = 0; s < l1b_fts->number_spectrometer(); s++)
      sza(s) = l1b_fts->solar_zenith(s).convert(units::deg).value;

    boost::shared_ptr<AtmosphereOco> atmosphere(boost::dynamic_pointer_cast<AtmosphereOco>(config_atmosphere));
    boa_rt.reset(new ChapmanBoaRT(atmosphere, sza, sb));
  }

  boost::shared_ptr<SpectralDomain> spec_domain;
  boost::shared_ptr<ChapmanBoaRT> boa_rt;
private:
};


BOOST_FIXTURE_TEST_SUITE(chapman_boa_rt, ChapmanBoaRtFixture)

BOOST_AUTO_TEST_CASE(check_jacobians)
{
  is_long_test();		// Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.

  // Load expected high res values
  Array<double, 1> wn_expt_high;
  Array<double, 1> refl_expt_high;
  Array<double, 2> jac_expt_high;
  IfstreamCs fd_high_res(test_data_dir() + "expected/chapman_boa_rt/high_res_fd");
  fd_high_res >> wn_expt_high
	      >> refl_expt_high
	      >> jac_expt_high;

  BOOST_CHECK_MATRIX_CLOSE(wn_expt_high, spec_domain->wavenumber());
  
  // Run reflectance only routine and due to how class is designed also the stokes routine
  Array<double, 1> refl_calc_high = 
    boa_rt->reflectance(*spec_domain, 0, true).spectral_range().data();

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt_high, refl_calc_high, 1e-5);

  // Run reflectance and jacobian routine and as a by product stokes_and_jacobian routine
  ArrayAd<double, 1> refl_jac_high = 
    boa_rt->reflectance(*spec_domain, 0).spectral_range().data_ad();

  BOOST_CHECK_MATRIX_CLOSE_TOL(refl_expt_high, refl_jac_high.value(), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_expt_high, refl_jac_high.jacobian(), 1e-6);

  if (false) {
    // Write out calculations for debugging
    std::cerr << std::scientific << std::setprecision(20)
	      << "# wavenumbers" << std::endl
	      << spec_domain->wavenumber() << std::endl
	      << "# reflectance" << std::endl
	      << refl_jac_high.value() << std::endl
	      << "# jacobian" << std::endl
	      << refl_jac_high.jacobian() << std::endl;
  }
  
}

BOOST_AUTO_TEST_CASE(generate_finite_diff_jac)
{
  // IMPORTANT: This case should be run to regenerate expected output should inputs ever change
  if(!getenv("L2_FP_RT_GEN")) {
    BOOST_WARN_MESSAGE(false, "Skipping generation case. To run, make rt_generate");
    return;
  }
  std::cerr << *config_state_vector << "\n";
  turn_on_logger();		// Have log output show up.
  
  // State vector object
  boost::shared_ptr<StateVector> sv_obj = config_state_vector;

  // Save initial state of statevector to set back at end of jacobian looping
  const Array<double, 1> initial_sv( sv_obj->state().copy() );

  // Retrieve unperturbed reflectance value
  Array<double, 1> refl_unpert = boa_rt->stokes(*spec_domain, 0)(Range::all(), 0);

  // Make sure statevector size is same as perturbation array
  if( initial_sv.extent(firstDim) != epsilon.extent(firstDim) ) {
    Exception err;
    err << "Statevector size does not match size of perturbation array, can not do FD jacobian calculations";
    throw(err);
  }

  // Loop over statevector perturbing each item in turn, saving value into jacobian array
  Array<double, 2> jac_fd(spec_domain->data().rows(), initial_sv.extent(firstDim));

  for( int sv_idx = 0; sv_idx < initial_sv.extent(firstDim); sv_idx++) {
    if( epsilon(sv_idx) > 0.0 ) {
      // Copy original and add perturbation for current element
      Array<double, 1> current_sv( initial_sv.copy() );
      current_sv(sv_idx) += epsilon(sv_idx);
    
      // Update statevector so calculations are done w/ perturbed values
      sv_obj->update_state(current_sv);
      
      // Calculate FD jacobian value
      jac_fd(Range::all(), sv_idx) = ( boa_rt->stokes(*spec_domain, 0)(Range::all(), 0) - refl_unpert(Range::all()) ) / epsilon(sv_idx);
    } else {
      // Set to all zeros if perturbation is 0.0 since 
      // the perturbation would have no effect
      jac_fd(Range::all(), sv_idx) = 0.0;
    }
  }

  // Set statevector object back to initial state
  sv_obj->update_state(initial_sv);
  
  // Will overwrite expected file for other test!!!
  std::ofstream fd_jac_file((test_data_dir() + "expected/chapman_boa_rt/high_res_fd").c_str());
  fd_jac_file << std::scientific << std::setprecision(20)
	      << "# wavenumbers" << std::endl
	      << spec_domain->wavenumber() << std::endl
	      << "# reflectance" << std::endl
	      << refl_unpert << std::endl
	      << "# jacobian" << std::endl
	      << jac_fd << std::endl;
  
}


BOOST_AUTO_TEST_SUITE_END()
