#include "unit_test_support.h"
#include "l_rad_driver.h"
#include "pressure_sigma.h"
#include "configuration_fixture.h"

#include "aerosol_property_hdf.h"

#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'
using namespace boost::assign; // bring 'operator+=()' into scope

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(l_rad_driver, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(lambertian_first_order)
{
    // Blitz tensor and range variables
    firstIndex i1; secondIndex i2; thirdIndex i3;
    Range ra(Range::all());

    int nstream = 8;
    int nstokes = 3;
    LRadDriver::PsMode ps_mode = LRadDriver::REGULAR;
     
    // Lambertian surface type
    int surface_type = 1;
   
    bool tms_correction = false;
    bool pure_nadir = false;

    // Load l_rad inputs from offline driver test data
    double sza; // theta0
    double zen; // theta
    double azm; // phi
    IfstreamCs geom_file(test_data_dir() + "in/l_rad_driver/geometry");
    geom_file >> sza >> zen >> azm;

    Array<double, 1> surface_parameters;
    IfstreamCs surf_params_file(test_data_dir() + "in/l_rad_driver/alb_01");
    surf_params_file >> surface_parameters;

    IfstreamCs wn_file(test_data_dir() + "in/l_rad_driver/wavenumbers");
    Array<double, 1> wavenumbers;
    wn_file >> wavenumbers;

    IfstreamCs od_ray_file(test_data_dir() + "in/l_rad_driver/od_ray_01");
    Array<double, 2> raydat;
    od_ray_file >> raydat;

    IfstreamCs od_aero_file(test_data_dir() + "in/l_rad_driver/od_aero_01");
    Array<double, 3> aerdat;
    od_aero_file >> aerdat;

    IfstreamCs od_gas_file(test_data_dir() + "in/l_rad_driver/od_gas_01");
    Array<double, 2> gasdat;
    od_gas_file >> gasdat;

    IfstreamCs ssa_file(test_data_dir() + "in/l_rad_driver/ssa_01");
    Array<double, 2> omega;
    ssa_file >> omega;

    IfstreamCs ssa_aero_file(test_data_dir() + "in/l_rad_driver/ssa_aero_01");
    Array<double, 3> omegaa;
    ssa_aero_file >> omegaa;

    std::vector<std::string> aerosol_types;
    IfstreamCs aero_types_file(test_data_dir() + "in/l_rad_driver/aerosol_types");
    std::string line;
    while(getline(aero_types_file, line)) {
        aerosol_types.push_back(line);
    }

    // Setup which wavenumber from inputs to use
    int wn_idx = 0;
    double wn = wavenumbers(wn_idx);

    // Sizes based on input data
    int n_layer = gasdat.cols();
    int n_aer = aerdat.depth();

    if ((int) aerosol_types.size() != n_aer) {
        Exception err;
        err << "Number of aerosol types from names file: " << aerosol_types.size() 
            << " does not match the number for data input: " << n_aer;
        throw err;
    }

    // Altitude grid sized for the number of levels
    Array<double, 1> altitude(n_layer+1);
    for (int i = 0; i < n_layer+1; i++) {
        altitude(i) = n_layer - i;
    }

    // Calculate total tau
    Array<double, 2> taudp(n_layer);
    for(int n = 0; n < n_layer; n++) {
        taudp(wn_idx, n) = sum(aerdat(wn_idx, n, Range::all())) + gasdat(wn_idx, n) + raydat(wn_idx, n);
    }

    // Calculate fractional rayleigh and aerosol for phase function computions
    Array<double, 1> frac_ray(n_layer);
    Array<double, 2> frac_aer(n_layer, n_aer);
    Array<double, 1> aer(n_aer);
    for(int n = 0; n < n_layer; n++) {
        double ray = raydat(wn_idx, n);
        for(int j = 0; j < n_aer; j++) {
            aer(j) = aerdat(wn_idx, n, j) * omegaa(wn_idx, n, j);
        }
        frac_ray(n) = ray / (ray + sum(aer));
        for(int j = 0; j < n_aer; j++) {
            frac_aer(n, j) = aer(j) / (ray + sum(aer));
        }
    }

    Array<double, 2> coefsr(RayleighGreekMoment::array());

    // Load aerosol properties from common aerosol input file
    HdfFile aerosol_prop_inp = HdfFile(input_dir() + "common/input/l2_aerosol_combined.h5", HdfFile::READ);
    Array<double, 1> a1(3), b(3);
    a1 = 0; b = 0.3, 0.6, 1.0;
    double psurf = 10;
    boost::shared_ptr<Pressure> p(new PressureSigma(a1,b, psurf, true));
    std::vector< boost::shared_ptr<AerosolProperty> > aer_properties;
    std::vector<Array<double, 2> > aer_pf;
    int s1 = 0;
    int s2 = 0;
    for(int aer_idx = 0; aer_idx < (int) aerosol_types.size(); aer_idx++) {
      aer_properties.push_back(boost::shared_ptr<AerosolProperty>(new AerosolPropertyHdf(aerosol_prop_inp, aerosol_types[aer_idx] + "/Properties", p)));
      aer_pf.push_back(aer_properties[aer_idx]->phase_function_moment_each_layer(wn).value()(ra, 0, ra));
      s1 = std::max(s1, aer_pf[aer_idx].rows());
      s2 = std::max(s2, aer_pf[aer_idx].cols());
    }

    // Phase function is sized according to the maximum number of moments and scattering coefficents for each aerosol type
    // Make sure to initialize to 0 since the computations that follow do not update very element
    ArrayAd<double, 3> phase_func(s1, n_layer, s2, 0);
    phase_func = 0;

    // Compute aerosol part of phase function
    // From aerosol.cc
    for(int j = 0; j < n_aer; ++j) {
        Range r1(0, aer_pf[j].rows() - 1);
        Range r2(0, aer_pf[j].cols() - 1);
        phase_func.value()(r1, ra, r2) += frac_aer(ra, j)(i2) * aer_pf[j](i1, i3);
    }

    // Rayleigh part of phase function
    // From atmosphere_oco.cc
    Range r3(0, coefsr.rows() - 1);
    phase_func.value()(r3, ra, ra) += frac_ray(i2) * coefsr(i1, i3);

    // Zeroth moment and scatterer is always 1
    phase_func(0, ra, 0) = 1;

    // Run l_rad with imputs
    boost::shared_ptr<LRadDriver> l_rad(new LRadDriver(nstream, nstokes, surface_type,
                                            tms_correction, pure_nadir, ps_mode));

    l_rad->setup_geometry(altitude, sza, zen, azm);
    l_rad->setup_surface_params(surface_parameters);

    // Must be called AFTER setup_geometry or else the correct values will not be used
    ArrayAd<double, 2> z_matrix(l_rad->z_matrix(phase_func));

    l_rad->setup_optical_inputs(taudp(wn_idx, ra), omega(wn_idx, ra), phase_func.value(), z_matrix.value());

    l_rad->clear_linear_inputs();
    l_rad->calculate_first_order();

    Array<double, 1> expt_stokes(3);
    expt_stokes =
        0.01471310439105843679, -0.00014581179209741254, 0.0;

    BOOST_CHECK_MATRIX_CLOSE_TOL(expt_stokes, l_rad->stokes(), 1e-7);
}

BOOST_AUTO_TEST_CASE(simple_brdf)
{
    int nstream = 4;
    int nmoms = 2*nstream;
    int nstokes = 3;
    LRadDriver::PsMode ps_mode = LRadDriver::REGULAR;
     
    // Lambertian surface type
    int surface_type = 3; // BREONVEG
    Array<double, 1> surface_params(5);
   
    bool tms_correction = false;
    bool pure_nadir = false;

    // Load l_rad inputs from offline driver test data
    double sza = 0.1; // theta0
    double zen = 0.0; // theta
    double azm = 0.0; // phi

    int nlayer = 1;
    Array<double, 1> heights(nlayer+1);
    ArrayAd<double, 1> od(nlayer, 1);
    ArrayAd<double, 1> ssa(nlayer, 1);
    double taug, taur;
    Range all = Range::all();
    
    // Simple height grid evenly spaced
    heights(0) = 100;
    for(int hidx = 1; hidx < nlayer+1; hidx++) {
      heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
    }
  
    // No aerosols, and depolization factor = 0 
    // so simplified phase function moments:
    ArrayAd<double, 3> phase_func(nmoms, nlayer, 6, 0);
    phase_func = 0.0;
    phase_func.value()(0, all, all) = 1.0;
    phase_func.value()(2, all, all) = 0.5;
  
    ////////////////
    // Surface only
    surface_params(0) = 1.0; // rahman kernel factor
    surface_params(1) = 0.1; // hotspot parameter
    surface_params(2) = 0.3; // asymmetry
    surface_params(3) = 1.5; // anisotropy_parameter
    surface_params(4) = 1.0; // breon kernel factor
 
    taur = 1.0e-6/nlayer;
    taug = 1.0e-6/nlayer;
  
    od = taur + taug;
    ssa = 0;
    ssa.value() = taur / od.value();

    // Run l_rad with imputs
    boost::shared_ptr<LRadDriver> l_rad(new LRadDriver(nstream, nstokes, surface_type,
                                                       tms_correction, pure_nadir, ps_mode));

    l_rad->setup_geometry(heights, sza, zen, azm);
    l_rad->setup_surface_params(surface_params);

    // Must be called AFTER setup_geometry or else the correct values will not be used
    ArrayAd<double, 2> z_matrix(l_rad->z_matrix(phase_func));
    z_matrix.resize_number_variable(1);

    l_rad->setup_optical_inputs(od.value(), ssa.value(), phase_func.value(), z_matrix.value());

    l_rad->setup_linear_inputs(od, ssa, phase_func, z_matrix);

    l_rad->calculate_first_order();

    double refl_calc = l_rad->stokes()(0);

    // Value from offline calculation, could also compare to LIDORT value
    double expt_val = 0.03540773173555763;
    BOOST_CHECK_CLOSE(expt_val, refl_calc, 1e-7);

    // Check surface jacobians against FD
  
    Array<double, 1> pert_values(surface_params.rows());
    pert_values = 1e-8;
  
    Array<double, 1> jac_surf( l_rad->surface_jacobian()(0, all) );
    Array<double, 1> jac_surf_fd( jac_surf.extent() );
    double refl_fd;
  
    jac_surf_fd = 0.0;
    refl_fd = 0.0;
  
    // First check PP mode against value just computed
    for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
        blitz::Array<double,1> surface_params_pert( surface_params.rows() );
        surface_params_pert = surface_params;
        surface_params_pert(p_idx) += pert_values(p_idx);
    
        l_rad->setup_surface_params(surface_params_pert);
        l_rad->calculate_first_order();
  
        refl_fd = l_rad->stokes()(0);
    
        jac_surf_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
    }

    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf, jac_surf_fd, 1e-6);

}

BOOST_AUTO_TEST_SUITE_END()
