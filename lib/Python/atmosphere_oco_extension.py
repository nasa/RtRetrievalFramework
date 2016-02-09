from __future__ import absolute_import
# The AtmosphereOco class is structured in a way that makes sense for the
# C++ code. However things tend to be deeply nested and fairly complicated
# structures (e.g., ArrayAdWithUnit). We provide a simpler interface here, 
# more suitable for use with python.

from .try_swig_load import *

def pressure_grid(self):
    '''Pressure grid, as an array of numbers in Pascals.'''
    return self.pressure.pressure_grid.convert("Pa").value.value

if(have_full_physics_swig):
    AtmosphereOco.pressure_grid = property(pressure_grid)

def __profile_plot(self, title, xlabel, x, y, fname = None):
    '''We have a number of similar profile plots. This captures this
    common behavior so we don't duplicate a lot of code. This either
    displays the plot, or if a file name is supplied writes to that file'''
    import matplotlib.pyplot as plt
    plt.title(title)
    plt.ylabel('Pressure (Pa) - Ticks are Pressure levels')
    plt.xlabel(xlabel)
    pg = self.pressure_grid
    plt.yticks(pg[::-1])
    plt.plot(y, x)
    # Reverse axis, so this show lower pressure higher up
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    if(fname is None):
        plt.show()
    else:
        # bbox_inches trims white space. We often want this for PNG.
        # We can play with this if we find out it is a problem
        plt.savefig(fname, bbox_inches=0, dpi=150)
        plt.clf()

if(have_full_physics_swig):
    AtmosphereOco.__profile_plot = __profile_plot

def temperature_func(self, p):
    '''Take a pressure value in Pascals, and return temperature in Kelvin.
    This is a more direct interface to the self.temperature object.'''
    x = AutoDerivativeWithUnitDouble(AutoDerivativeDouble(p), "Pa")
    return self.temperature.temperature(x).convert("K").value.value

if(have_full_physics_swig):
    AtmosphereOco.temperature_func = temperature_func

def temperature_profile_plot(self, fname = None):
    '''Quick method to plot out the temperature profile'''
    pg = self.pressure_grid
    x = np.linspace(pg[0], pg[-1], 100)
    y = [ self.temperature_func(v) for v in x ]
    self.__profile_plot('Temperature Profile', 'Temperature (K)', x, y, fname)
    
if(have_full_physics_swig):
    AtmosphereOco.temperature_profile_plot = temperature_profile_plot

def gravity_func(self, spec_index, p):
    '''Takes the pressure in Pascals and return the gravity in m/s^2.

    This is indexed by the spectral index, since different spectrometers
    may have different gravity (in practice they are all the same for GOSAT).

    This is a more direct interface to the self.altitude_obj object.'''
    x = AutoDerivativeWithUnitDouble(AutoDerivativeDouble(p), "Pa")
    return self.altitude_obj[spec_index].gravity(x).convert("m/s^2").value.value
    
if(have_full_physics_swig):
    AtmosphereOco.gravity_func = gravity_func

def gravity_profile_plot(self, fname = None):
    '''Quick method to plot out the gravity profile'''
    pg = self.pressure_grid
    x = np.linspace(pg[0], pg[-1], 100)
    y = [ self.gravity_func(0, v) for v in x ]
    self.__profile_plot('Gravity Profile (Spectrometer 0)', 
                        'Gravity (m/s^2)', x, y, fname)
    
if(have_full_physics_swig):
    AtmosphereOco.gravity_profile_plot = gravity_profile_plot

def volume_mixing_ratio_func(self, gas_name, p):
    '''Takes the pressure in Pascals and return the volumn mixing ratio.

    This is indexed by the gas name, usually "CO2", "O2", "H2O"

    This is a more direct interface to the self.absorber.absorber_vmr object.'''
    x = AutoDerivativeDouble(p)
    return self.absorber.absorber_vmr(gas_name).volume_mixing_ratio(x).value

if(have_full_physics_swig):
    AtmosphereOco.volume_mixing_ratio_func = volume_mixing_ratio_func

def volume_mixing_ratio_profile_plot(self, fname = None):
    '''Quick method to plot out the volume mixing ratio profile'''
    pg = self.pressure_grid
    x = np.linspace(pg[0], pg[-1], 100)
    for specie in ("CO2", "O2", "H2O"):
        y = [ self.volume_mixing_ratio_func(specie, v) for v in x ]
        if(fname is None):
            fname_specie = None
        else:
            fname_specie = "%s_%s.png" % (fname, specie)
        self.__profile_plot('%s Profile' % specie, 'Volume Mixing Ratio', x, y, fname_specie)

if(have_full_physics_swig):
    AtmosphereOco.volume_mixing_ratio_profile_plot = volume_mixing_ratio_profile_plot

def absorber_integrand_independent_wn(self, spec_index, species_index, p):
    '''Return the integrand that is independent of wave number in the
    absorber. This is in cm^-2 / Pa'''
    x = DoubleWithUnit(p, "Pa")
    return self.absorber.integrand_independent_wn(spec_index, species_index, x).convert("cm^-2 / Pa").value.value

if(have_full_physics_swig):
    AtmosphereOco.absorber_integrand_independent_wn = absorber_integrand_independent_wn

def absorber_integrand_independent_wn_profile_plot(self, fname = None):
    '''Quick method to plot out the absorber_integrand_independent_wn'''
    pg = self.pressure_grid
    x = np.linspace(pg[0], pg[-1], 100)
    for specie in ("CO2", "O2", "H2O"):
        i = self.absorber.gas_index(specie)
        y = [ self.absorber_integrand_independent_wn(0, i, v) for v in x ]
        if(fname is None):
            fname_specie = None
        else:
            fname_specie = "%s_%s.png" % (fname, specie)
        self.__profile_plot('%s Integrand ind WN Profile' % specie, 
                            'Integrand ind WN', x, y, fname_specie)

if(have_full_physics_swig):
    AtmosphereOco.absorber_integrand_independent_wn_profile_plot = absorber_integrand_independent_wn_profile_plot

def absorber_integrand(self, wn, p, spec_index, species_index):
    '''Return the integrand. Pressure should be in Pa, and wave number in cm^-1. Returned data is in Pa^-1'''
    return self.absorber.integrand(wn, p, spec_index, species_index)
    
if(have_full_physics_swig):
    AtmosphereOco.absorber_integrand = absorber_integrand

def absorber_integrand_profile_plot(self, spec_index, wn, fname = None):
    '''Quick method to plot out the absorber_integrand_independent_wn'''
    pg = self.pressure_grid
    x = np.linspace(pg[0], pg[-1], 100)
    for specie in ("CO2", "O2", "H2O"):
        i = self.absorber.gas_index(specie)
        y = [ self.absorber_integrand(wn, v, spec_index, i) for v in x ]
        if(fname is None):
            fname_specie = None
        else:
            fname_specie = "%s_%s.png" % (fname, specie)
        self.__profile_plot('%s Integrand WN=%d Profile' % (specie, wn), 
                            'Integrand', x, y, fname_specie)

if(have_full_physics_swig):
    AtmosphereOco.absorber_integrand_profile_plot = absorber_integrand_profile_plot
