from full_physics import *

from matplotlib import pyplot
import numpy as np
import re

def get_sv_item_indexes(sv, sv_check_items):
    sv_names = sv.state_vector_name()

    item_sv_idx = []
    for nm_re in sv_check_items:
        found_idx = [y[0] for y in filter(lambda x: re.search(nm_re, x[1]), enumerate(sv_names))]
        item_sv_idx.append(found_idx)
    return item_sv_idx

def check_aer_coeff_fd_jac(fm, sv_check_items, item_perts, do_forward_fd=False):
    sv = fm.state_vector()

    sv_item_idxs = get_sv_item_indexes(sv, sv_check_items)

    spec_idx = 0
    print >>sys.stderr, "Calculating analytic jacobian"
    spectrum = fm.radiance(spec_idx, False)
    rad_orig = np.copy(spectrum.spectral_range().data_ad().value())
    jac_orig = np.copy(spectrum.spectral_range().data_ad().jacobian())

    # MUST copy original state vector here since the memory pointed
    # to by sv.state() will change with each call to sv.update_state(...)
    sv_orig = np.copy(sv.state())

    fd_jac = np.zeros(jac_orig.shape, dtype=float)
    for i_indexes, i_name, pert in zip(sv_item_idxs, sv_check_items, item_perts):

        sv_pert = np.zeros(sv_orig.shape, dtype=float)
        # Lots of copies here to make sure we are not referring to memory that
        # has changed on us in a subsequent step
        for count, c_index in enumerate(i_indexes):
            print >>sys.stderr, "SV value o = ", sv_orig[c_index]
            # Use more time consuming central difference method
            print >>sys.stderr,"\nCalculating fd +jacobian for %s #%d" % (i_name, count+1)
            print >>sys.stderr,"Coeff: %e, Pert: %e" % (sv_orig[c_index], pert[count])            
            sv_pert = np.copy(sv_orig)
            sv_pert[c_index] += pert[count]
            sv.update_state(sv_pert)
            spectrum_p = fm.radiance(spec_idx, False)
            spec_fd_p = np.copy(spectrum_p.spectral_range().data_ad().value())

            print >>sys.stderr,"\nCalculating fd -jacobian for %s #%d" % (i_name, count+1)
            print >>sys.stderr,"Coeff: %e, Pert: %e" % (sv_orig[c_index], -pert[count])
            sv_pert = np.copy(sv_orig)
            sv_pert[c_index] -= pert[count]
            sv.update_state(sv_pert)
            spectrum_m = fm.radiance(spec_idx, False)
            spec_fd_m = np.copy(spectrum_m.spectral_range().data_ad().value())

            fd_jac[:, c_index]  = (spec_fd_p - spec_fd_m) / (2*pert[count])

            if do_forward_fd:
                forward_jac = (spec_fd_p - rad_orig) / pert[count]

            jac_diff = jac_orig[:,c_index] - fd_jac[:, c_index]
            max_a = np.max(np.abs(jac_orig[:,c_index]))
            max_f = np.max(np.abs(fd_jac[:, c_index]))
            perc_diff = np.mean(np.abs(jac_diff) / max(max_a,max_f))
            print "Mean % Diff =", perc_diff
            rms_diff = np.sqrt(np.mean(np.square(jac_diff)))
            print "RMS Diff =", rms_diff

            # Plot stuff
            pyplot.figure(count+1)
            pyplot.clf()
            pyplot.subplot(2, 1, 1)
            pyplot.plot(jac_orig[:,c_index])
            pyplot.plot(fd_jac[:,c_index])
            legend_names = ["a","fd_c"]

            if do_forward_fd:
                pyplot.plot( forward_jac )
                legend_names.append("fd_f")

            pyplot.legend(legend_names, 0)
            pyplot.title("%s #%d Value: %.4e, Pert: %.2e\n%% Diff = %e, RMS Diff = %e" % (i_name, count+1, sv_orig[c_index], pert[count], perc_diff, rms_diff))

            pyplot.subplot(2, 1, 2)
            pyplot.plot(jac_orig[:,c_index] - fd_jac[:,c_index])
            legend_names = ["a - fd_c"]

            if do_forward_fd:
                pyplot.plot(jac_orig[:,c_index] - forward_jac)
                legend_names.append("a - fc_f")

            pyplot.legend(legend_names, 0)

            plot_file = "%s_%02d_%02d_fd.png" % (i_name.replace(" ", "_"), count+1, spec_idx)
            pyplot.savefig(plot_file)
