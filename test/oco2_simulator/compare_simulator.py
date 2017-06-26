from full_physics import *
import numpy as np
from matplotlib.pylab import *
from matplotlib.backends.backend_pdf import PdfPages
import os

pp = PdfPages('oco2_simulator/comparison_plot.pdf')

bdir = os.path.dirname(__file__) + "/"
sdir = bdir + "simulator_result/"
sid_name = "2010090912004075"
met_name = "../input/oco2_sim_met.h5"
l1b_name = "../input/oco2_sim_l1b.h5"
scene_file = "../input/oco2_sim_scene.h5"
r = L2Run(bdir + "config/config_orbit_sim_match.lua", sid_name,met_name, 
          l1b_name, scene_file = scene_file)

print (r.forward_model)
print (r.xco2)

# High resolution reflectance from CSU
refl = []
# High resolution radiance from CSU
rad = []
# Spectral domain for high resolution data
hres_sd = []
hres_ref = np.genfromtxt(sdir + "hires_10774_1.txt", skip_header = 1)
sd = SpectralDomain(hres_ref[:,0], Unit("nanometer"))
# Conversion currently needed by ILS. This is really a bug, which
# we will fix as Ticket #1260. But for now, work around
sd = SpectralDomain(sd.wavelength(Unit("micron")), Unit("micron"))
hres_sd.append(sd)
refl.append(hres_ref[:,2])
sr = SpectralRange(hres_ref[:,2] * hres_ref[:,3],
                   Unit("Ph / sec / m^2 / sr / um"))
rad.append(Spectrum(sd, sr))

hres_ref = np.genfromtxt(sdir + "hires_10774_2.txt", skip_header = 1)
sd = SpectralDomain(hres_ref[:,0], Unit("nanometer"))
sd = SpectralDomain(sd.wavelength(Unit("micron")), Unit("micron"))
hres_sd.append(sd)
refl.append(hres_ref[:,2])
sr = SpectralRange(hres_ref[:,2] * hres_ref[:,3],
                   Unit("Ph / sec / m^2 / sr / um"))
rad.append(Spectrum(sd, sr))

hres_ref = np.genfromtxt(sdir + "hires_10774_3.txt", skip_header = 1)
sd = SpectralDomain(hres_ref[:,0], Unit("nanometer"))
sd = SpectralDomain(sd.wavelength(Unit("micron")), Unit("micron"))
hres_sd.append(sd)
refl.append(hres_ref[:,2])
sr = SpectralRange(hres_ref[:,2] * hres_ref[:,3],
                   Unit("Ph / sec / m^2 / sr / um"))
rad.append(Spectrum(sd, sr))

# Solar absorption from Level 2
sabs = []
# Solar absorption from CSU
sabssim = []
# Full solar spectrum from Level 2
sspec = []
# Full solar spectrum from CSU
sspecsim = []

sm = r.solar_model[0]
ds = sm.doppler_shift
sabs.append(sm.absorption_spectrum.solar_absorption_spectrum(ds.doppler_stretch(hres_sd[0])))
sspec.append(sm.solar_spectrum(hres_sd[0]))
sabssim.append(np.genfromtxt(sdir + "solar_spectrum_10774_1.txt", skip_header = 5)[:,3])
sspecsim.append(np.genfromtxt(sdir + "solar_spectrum_10774_1.txt", skip_header = 5)[:,4])
sm = r.solar_model[1]
ds = sm.doppler_shift
sabs.append(sm.absorption_spectrum.solar_absorption_spectrum(ds.doppler_stretch(hres_sd[1])))
sspec.append(sm.solar_spectrum(hres_sd[1]))
sabssim.append(np.genfromtxt(sdir + "solar_spectrum_10774_2.txt", skip_header = 5)[:,3])
sspecsim.append(np.genfromtxt(sdir + "solar_spectrum_10774_2.txt", skip_header = 5)[:,4])
sm = r.solar_model[2]
ds = sm.doppler_shift
sabs.append(sm.absorption_spectrum.solar_absorption_spectrum(ds.doppler_stretch(hres_sd[2])))
sspec.append(sm.solar_spectrum(hres_sd[2]))
sabssim.append(np.genfromtxt(sdir + "solar_spectrum_10774_3.txt", skip_header = 5)[:,3])
sspecsim.append(np.genfromtxt(sdir + "solar_spectrum_10774_3.txt", skip_header = 5)[:,4])
figure(1)
subplot(311)
title("Difference solar transmittance band 1")
xlabel("Wavelength")
ylabel("Difference %")
plot(sabs[0].wavelength, (sabs[0].value -sabssim[0]) / max(sabssim[0]) * 100.0)

subplot(312)
title("Difference solar transmittance band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(sabs[1].wavelength, (sabs[1].value -sabssim[1]) / max(sabssim[1]) * 100.0)

subplot(313)
title("Difference solar transmittance band 3")
xlabel("Wavelength")
ylabel("Difference %")
plot(sabs[2].wavelength, (sabs[2].value -sabssim[2]) / max(sabssim[2]) * 100.0)
pp.savefig()

figure(2)
subplot(311)
title("Difference full solar spectrum band 1")
xlabel("Wavelength")
ylabel("Difference %")
plot(sspec[0].wavelength, (sspec[0].value -sspecsim[0]) / max(sspecsim[0]) * 100.0)

subplot(312)
title("Difference full solar spectrum band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(sspec[1].wavelength, (sspec[1].value -sspecsim[1]) / max(sspecsim[1]) * 100.0)

subplot(313)
title("Difference full solar spectrum band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(sspec[2].wavelength, (sspec[2].value -sspecsim[2]) / max(sspecsim[2]) * 100.0)

pp.savefig()

rcalc = []
rcalc.append(r.radiative_transfer.reflectance(hres_sd[0], 0, True).value)
rcalc.append(r.radiative_transfer.reflectance(hres_sd[1], 1, True).value)
rcalc.append(r.radiative_transfer.reflectance(hres_sd[2], 2, True).value)
figure(3)
subplot(311)
title("Difference High res spectra band 1")
xlabel("Wavelength")
ylabel("Difference %")
plot(hres_sd[0].wavelength(), (refl[0] - rcalc[0]) / max(refl[0]) * 100.0)
subplot(312)
title("Difference High res spectra band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(hres_sd[1].wavelength(), (refl[1] - rcalc[1]) / max(refl[1]) * 100.0)
subplot(313)
title("Difference High res spectra band 3")
xlabel("Wavelength")
ylabel("Difference %")
plot(hres_sd[2].wavelength(), (refl[2] - rcalc[2]) / max(refl[2]) * 100.0)

# Low resolution calculated
lcalc = []
r.forward_model.setup_grid()
lcalc.append(r.forward_model.radiance(0, True))
lcalc.append(r.forward_model.radiance(1, True))
lcalc.append(r.forward_model.radiance(2, True))
# Low resoluiton from CSU
mrad = []
mrad.append(r.forward_model.measured_radiance(0))
mrad.append(r.forward_model.measured_radiance(1))
mrad.append(r.forward_model.measured_radiance(2))

# Directly apply ILS to high resolution spectrum from CSU
lrad_ils = []
plist = r.spectral_window.grid_indexes(r.instrument.pixel_spectral_domain(0), 0)
lrad_ils.append(r.instrument.apply_instrument_model(rad[0], plist, 0))
plist = r.spectral_window.grid_indexes(r.instrument.pixel_spectral_domain(1), 1)
lrad_ils.append(r.instrument.apply_instrument_model(rad[1], plist, 1))
plist = r.spectral_window.grid_indexes(r.instrument.pixel_spectral_domain(2), 2)
lrad_ils.append(r.instrument.apply_instrument_model(rad[2], plist, 2))

pp.savefig()
figure(4)
subplot(311)
title("Difference Low res spectra band 1")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[0].wavelength, (mrad[0].value - lcalc[0].value) / max(mrad[0].value) * 100.0)
subplot(312)
title("Difference Low res spectra band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[1].wavelength, (mrad[1].value - lcalc[1].value) / max(mrad[1].value) * 100.0)
subplot(313)
title("Difference Low res spectra band 3")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[2].wavelength, (mrad[2].value - lcalc[2].value) / max(mrad[2].value) * 100.0)
pp.savefig()
figure(5)
subplot(311)
title("Difference direct ISL res spectra band 1")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[0].wavelength, (mrad[0].value - lrad_ils[0].value) / max(mrad[0].value) * 100.0)
subplot(312)
title("Difference direct ISL res spectra band 2")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[1].wavelength, (mrad[1].value - lrad_ils[1].value) / max(mrad[1].value) * 100.0)
subplot(313)
title("Difference direct ISL res spectra band 3")
xlabel("Wavelength")
ylabel("Difference %")
plot(mrad[2].wavelength, (mrad[2].value - lrad_ils[2].value) / max(mrad[2].value) * 100.0)
pp.savefig()
pp.close()
# Can display directly, but don't want to do this in the automated tests.
# show()







