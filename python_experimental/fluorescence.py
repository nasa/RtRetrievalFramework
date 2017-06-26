from full_physics import *
import numpy as np
import matplotlib.pyplot as pp
import sys
import math

from progressbar import ProgressBar, Percentage, Bar

class FluorescenceEffectPy(SpectrumEffectImpBase):
    def __init__(self, f_coeffs, ret_flag, atm, lza, spec_index, reference_wn, source_units):
        SpectrumEffectImpBase.__init__(self, f_coeffs, ret_flag)
        self.atm = atm
        self.lza = lza
        self.spec_index = spec_index
        self.reference_wn = reference_wn
        self.source_units = source_units

        # Save some values for plotting or writing to disk
        self.orig_radiance = None
        self.f_radiance = None
        self.f_contrib = None

    def clone(self):
        res = FluorescenceEffect(self.coefficient().value(), self.used_flag_value(), self.atm, self.lza)
        return res

    def apply_effect(self, spectrum):
        lza_sec = 1.0/np.cos(self.lza.convert(Unit("rad")).value)

        sys.stderr.write("Adding fluorescence to band: %d\n" % self.spec_index)

        wn_arr = spectrum.spectral_domain().wavenumber()
        rad_arr = spectrum.spectral_range().data_ad()

        self.orig_radiance = Spectrum(spectrum.spectral_domain(), spectrum.spectral_range())

        conv_factor = conversion(self.source_units, spectrum.spectral_range().units())

        # Always keep flouresence at reference point 0 or positive
        fs_ref = self.coefficient()[0] * conv_factor
        slope = self.coefficient()[1]

        f_contrib_ad = ArrayAd_double_1(rad_arr.rows(), rad_arr.number_variable())
        progress = ProgressBar(widgets=[Percentage(), Bar()], maxval=wn_arr.shape[0], term_width=51).start()
        for wn_idx in range(wn_arr.shape[0]):
            o2_col_abs = self.atm.column_optical_depth(wn_arr[wn_idx], self.spec_index, "O2")

            f_surf = fs_ref * (1.0 + slope * (wn_arr[wn_idx] - self.reference_wn))
            f_contrib_ad[wn_idx] = f_surf * exp(-1 * o2_col_abs * lza_sec)
            rad_arr[wn_idx] += f_contrib_ad[wn_idx]

            progress.update(wn_idx+1)
        progress.finish()

        self.f_contrib = Spectrum(spectrum.spectral_domain(), SpectralRange(f_contrib_ad, spectrum.spectral_range().units()))
        self.f_radiance = Spectrum(spectrum.spectral_domain(), spectrum.spectral_range())
                                 
    def state_vector_name_i(self, i):
        return "Fluorescence Surface Coefficient %d" % i

    def desc(self):
        return '''
FluorescenceEffect
   Coefficient:    %s
   Retrieval flag: %s
''' % ( self.coefficient().value().__str__(), self.used_flag_value().__str__())

class FluorescenceOutput(RegisterOutputBase):
    def __init__(self, f_effect, spec_win, instrument, conv_units):
        RegisterOutputBase.__init__(self)
        self.f_effect = f_effect
        self.spec_win = spec_win
        self.instrument = instrument
        self.conv_units = conv_units

    def convolved_fluorescence(self):
        spec_idx = 0
        plist = self.spec_win.grid_indexes(self.instrument.pixel_spectral_domain(spec_idx), spec_idx)
        
        print >>sys.stderr, "Convolving spectra"
        conv_f_contrib = self.instrument.apply_instrument_model(self.f_effect.f_contrib, plist, spec_idx)

        print >>sys.stderr, "Making new spectrum"
        conv_f_contrib = Spectrum(conv_f_contrib.spectral_domain(), conv_f_contrib.spectral_range().convert(self.conv_units))

        print "Returning contrib"
        return conv_f_contrib

    def register_output(self, out):
        def wavenumber_high_res(f_contrib):
            return f_contrib.spectral_domain().wavenumber()
        out.register_array_1d("/SpectralParameters/wavenumber_high_res",
                              lambda : wavenumber_high_res(self.f_effect.f_contrib))

        def fluorescence_high_res(f_contrib):
            return f_contrib.spectral_range().data()
        out.register_array_1d("/SpectralParameters/fluorescence_high_res",
                              lambda : fluorescence_high_res(self.f_effect.f_contrib))

        def wavenumber_conv_fluorescence(conv_f):
            return conv_f.spectral_domain().wavenumber()
        out.register_array_1d("/SpectralParameters/wavenumber_conv_fluorescence",
                              lambda : wavenumber_conv_fluorescence(self.convolved_fluorescence()))

        def fluorescence_convolved(conv_f):
            return conv_f.spectral_range().data()
        out.register_array_1d("/SpectralParameters/fluorescence_convolved",
                              lambda : fluorescence_convolved(self.convolved_fluorescence()))
        
    def register_output_apriori(self, out):
        pass
        
class FluorescencePlot(ObserverStateVector):
    def __init__(self, f_effect, solver=None):
        ObserverStateVector.__init__(self)
        self.f_effect = f_effect
        self.solver = solver

    def notify_add(self, sv):
        pass

    def notify_remove(self, sv):
        pass

    def notify_update(self, sv):
        self.plot()

    def plot(self):
        if self.f_effect.f_contrib == None or (abs(self.f_effect.coefficient().value()[0]) <= 0 and abs(self.f_effect.coefficient().value()[1]) <= 0):
            return

        pp.clf()
        pp.subplot(211)
        pp.plot(self.f_effect.orig_radiance.spectral_domain().wavenumber(), self.f_effect.orig_radiance.spectral_range().data())

        pp.plot(self.f_effect.f_radiance.spectral_domain().wavenumber(), self.f_effect.f_radiance.spectral_range().data())
        pp.legend(["no f", "with f"], 0)
        pp.title("Fluorescence")

        pp.subplot(212)
        pp.plot(self.f_effect.f_contrib.spectral_domain().wavenumber(), self.f_effect.f_contrib.spectral_range().data())
        pp.legend(["Ftoa"], 0)

        if self.solver != None:
            iteration = self.solver.number_iteration()
            pp.savefig("fluorescence_hr_i%02d.png" % (iteration-1))
        else:
            pp.savefig("fluorescence_hr.png")
