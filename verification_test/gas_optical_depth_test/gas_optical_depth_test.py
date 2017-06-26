from __future__ import print_function
from nose.tools import *
from full_physics import *
from nose.plugins.skip import Skip, SkipTest
import os
import matplotlib.pyplot as plt

gosat_config = os.path.dirname(__file__) + "/../../input/gosat/config/config.lua"
tccon_small_set = os.path.dirname(__file__) + "/../../test/tccon_small_set/"
ecmwf_file = tccon_small_set + "acos_EcmB2900_tccon_5_good_qual.h5"
spectrum_file = tccon_small_set + "acos_L1bB2900_tccon_5_good_qual.h5"
cloud_file = tccon_small_set + "acos_CldB2900_tccon_5_good_qual.h5"
sounding_id = "20100223034944"

if(have_full_physics_swig):
    r = L2Run(gosat_config, sounding_id, ecmwf_file, spectrum_file,
              cloud_file)

def test_plot():
    '''Generate various profile plots.'''
    if(not have_full_physics_swig):
        raise SkipTest
    try:
        os.mkdir("gas_optical_depth")
    except OSError:
        pass
    r.atmosphere.temperature_profile_plot("gas_optical_depth/temperate.png")
    r.atmosphere.gravity_profile_plot("gas_optical_depth/gravity.png")
    r.atmosphere.volume_mixing_ratio_profile_plot("gas_optical_depth/vmr")
    r.atmosphere.absorber_integrand_independent_wn_profile_plot("gas_optical_depth/integrand_ind_wn")
    r.atmosphere.absorber_integrand_profile_plot(0, 12930.09, "gas_optical_depth/integrand_wn_0")
    r.atmosphere.absorber_integrand_profile_plot(1, 6146.14, "gas_optical_depth/integrand_wn_1")

def test_top_layer():
    '''Some detailed plots of particular trouble areas'''
    if(not have_full_physics_swig):
        raise SkipTest
    try:
        os.mkdir("gas_optical_depth")
    except OSError:
        pass
    wnmax = 13056.32
    func = lambda p : r.atmosphere.absorber_integrand(wnmax, p, 0, 1)
    pg = r.atmosphere.pressure_grid
    x = np.linspace(pg[0],pg[1], 100)
    y =[ func(v) for v in x]
    t = np.array([0,1,2,3,4,5,6,7,8,9])
    sublay_frac = (t + 0.5) / 10
    xsub = pg[0] * (1 - sublay_frac) + pg[1] * sublay_frac
    ysub = [ func(v) for v in xsub ]
    pecmwf = []
    for p in r.atmosphere.temperature.pressure_profile:
        if(p >= pg[0] and p <= pg[1]):
            pecmwf.append(p)
    yecmwf = [ func(v) for v in pecmwf ]
    plt.title("Integrand vs. Pressure Top Layer")
    plt.ylabel('Pressure (Pa)')
    plt.xlabel('Integrand')
    plt.plot(y, x)
    plt.scatter(ysub, xsub, marker='D', c='r', s=60)
    plt.scatter(yecmwf, pecmwf, marker='*', c='g', s=60)
    plt.legend(("Integrand", "Sublayer", "ECMWF Pressure"), loc='center right')
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.savefig("gas_optical_depth/top_layer.png", bbox_inches=0, dpi=150)
    plt.clf()

    x = np.linspace(pg[0],300, 100)
    y =[ func(v) for v in x]
    plt.title("Integrand vs. Pressure Top Layer (Zoomed in)")
    plt.ylabel('Pressure (Pa)')
    plt.xlabel('Integrand')
    plt.plot(y, x)
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ysub = [ func(v) for v in xsub ]
    plt.scatter(ysub, xsub, marker='D', c='r', s=60)
    plt.scatter(yecmwf, pecmwf, marker='*', c='g', s=60)
    plt.legend(("Integrand", "Sublayer", "ECMWF Pressure"), loc='center left')
    plt.savefig("gas_optical_depth/top_layer_zoom.png", bbox_inches=0, dpi=150)
    plt.clf()
    
def diff_integrate(spec_index, wn):
    '''Difference between what we get doing a full direct integration vs.
    our normal approximation'''
    odapprox = r.atmosphere.absorber.optical_depth_each_layer(wn, spec_index).value
    od = r.atmosphere.absorber.optical_depth_each_layer_direct_integrate(wn, spec_index, 0, 1e-2)
    return abs(odapprox - od) / np.where(odapprox == 0, 1, odapprox) * 100.0

def test_full_range():
    '''This tests the full spectral range, looking for the maximum difference
    between the full integration and our approximation. Note that this takes
    a fair chunk of time to run, because the full integration is slow.'''
    raise SkipTest
    for spec_index in range(3):
        mx = np.zeros((19, 3))
        mx_ind = np.zeros((19, 3))
        i = 0
        print("Starting spec_index", spec_index)
        for wn in r.forward_model.spectral_domain(spec_index).data:
            if(i % 100 == 0):
                print("Doing", i)
                print("Current max (species x layer):")
                print(mx)
                print("Index with current max (species x layer):")
                print(mx_ind)
            t = diff_integrate(spec_index, wn)
            mx_ind = np.where(mx < t, i, mx_ind)
            mx = np.where(mx < t, t, mx)
            i += 1
        print("Final results spec_index", spec_index, ":", mx, mx_ind)
    
def test_short_range():
    '''This is a small test that just makes sure we can still run the full
    range test, which we normally skip because it takes a long time to run.'''
    if(not have_full_physics_swig):
        raise SkipTest
    spec_index = 0
    print("Starting spec_index", spec_index)
    i = 3920
    mx = np.zeros((19, 3))
    mx_ind = np.zeros((19, 3))
    for wn in r.forward_model.spectral_domain(spec_index).data[3920:3940]:
        t = diff_integrate(spec_index, wn)
        mx_ind = np.where(mx < t, i, mx_ind)
        mx = np.where(mx < t, t, mx)
        i += 1
    print("Final results spec_index", spec_index, ":", mx, mx_ind)
    assert mx.max() < 2.5e-1
    
