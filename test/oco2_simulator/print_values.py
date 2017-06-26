from full_physics import *
import numpy as np
import os

bdir = os.path.dirname(__file__) + "/"
sdir = bdir + "simulator_result/"
sid_name = "2010090912004075"
met_name = "../input/oco2_sim_met.h5"
l1b_name = "../input/oco2_sim_l1b.h5"
scene_file = "../input/oco2_sim_scene.h5"
r = L2Run(bdir + "config/config_orbit_sim_match.lua", sid_name,met_name, 
          l1b_name, scene_file = scene_file)

hres_ref = np.genfromtxt(sdir + "hires_10774_1.txt", skip_header = 1)
sd = SpectralDomain(hres_ref[:,0], Unit("nanometer"))
sd = SpectralDomain(sd.wavelength(Unit("micron")), Unit("micron"))
refl = hres_ref[:,2]
rcalc = r.radiative_transfer.reflectance(sd, 0, True).value
stkcalc = r.radiative_transfer.stokes(sd, 0)
maxindex = np.argmax(abs(refl-rcalc))
wnmax = sd.wavenumber()[maxindex]
print("Maximum difference:\n")
print("Wavenumber:      ", wnmax)
print("CSU:             ", refl[maxindex])
print("JPL (LSI):       ", rcalc[maxindex])
print("Diff %:          ", (refl[maxindex] - rcalc[maxindex]) / max(refl) * 100.0)

# Calculate using high stream
lrad_rt = r.radiative_transfer.high_stream_radiative_transfer.rt
lidort_rt = lrad_rt.radiative_transfer
sd_max = SpectralDomain([wnmax])
rcalc_hs = lrad_rt.reflectance(sd_max, 0).value[0]
lrad_stoke = lrad_rt.stokes(sd_max, 0)[0,:]
print("Stokes coeff:    ", lrad_rt.stokes_coefficient.stokes_coefficient.value[0,0:3])
print("JPL High Stream (no LSI): ")
print("                 ", rcalc_hs)
print("JPL Stokes (LSI):", stkcalc[maxindex])
print("JPL Low Stream Stokes: ")
print("                 ", r.radiative_transfer.low_stream_radiative_transfer.stokes(sd_max,0)[0,:])
print("JPL High Stream Stokes (include LRAD):      ")
print("                 ", lrad_stoke)
print("LIDORT High Stream only (no LRAD, multiscatter only):")
print("                 ", lidort_rt.stokes(sd_max, 0)[0,0:3])
atm = lidort_rt.atmosphere
print("Altitude (KM): ")
print(atm.altitude(0).value.value)
print("Pressure: ")
print(atm.pressure.pressure_grid.value.value)
print("Optical depth: ")
print(atm.optical_depth_wrt_iv(wnmax, 0).value)

print("  Optical depth Rayleigh: ")
print(atm.rayleigh.optical_depth_each_layer(wnmax, 0).value)
for i in range(atm.absorber.number_species):
    print("  Optical depth ", atm.absorber.gas_name(i))
    print(atm.absorber.optical_depth_each_layer(wnmax,0).value[:, i])
print("Single Scattering Albedo: ")
print(atm.single_scattering_albedo_wrt_iv(wnmax, 0).value)
for i in range(19):
    print("Scattering moment layer", i, "(should be Rayleigh only)")
    print(atm.scattering_moment_wrt_iv(wnmax, 0).value[:,i,:])


print("O2 concentration")
print(atm.absorber.absorber_vmr("O2").coefficient.value[0])
print("Specific humidity layers:")
print(atm.absorber.specific_humidity_layer.value.value)

print("Dry air column thickness layer (molecules / m^2):")
print(atm.absorber.dry_air_column_thickness_layer.value.value)
print("Total: ", np.sum(atm.absorber.dry_air_column_thickness_layer.value.value))
print("O2 column thickness layer (molecules / m^2):")
print(atm.absorber.gas_column_thickness_layer("O2").value.value)
print("Total: ", atm.absorber.gas_total_column_thickness("O2").value.value)

