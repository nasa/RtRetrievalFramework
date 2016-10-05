# This script runs in the context of the analyze_l2.py environment and can not be run independently

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

pdf_filename = "aerosol_shape_treatment_1_analysis.pdf"
print "Creating plot file: %s" % pdf_filename
pdf_out_obj = PdfPages(pdf_filename)

def save_pdf_page(res, which_pdf=pdf_out_obj):
    if hasattr(res, "__iter__"):
        for item in res:
            save_pdf_page(item, which_pdf)
    else:
        which_pdf.savefig(res.figure)
        pyplot.close(res.figure)
    return res

save_pdf_page( plot_xco2_diff_corr() )
save_pdf_page( plot_psurf_diff_corr() )
save_pdf_page( plot_iter_histogram() )
save_pdf_page( plot_div_histogram() )
save_pdf_page( plot_abo2_chi2_diff_corr() )
save_pdf_page( plot_wco2_chi2_diff_corr() )
save_pdf_page( plot_sco2_chi2_diff_corr() )
save_pdf_page( plot_aerosol_total_diff_corr() )
save_pdf_page( plot_aerosol_kahn_2b_diff_corr() )
save_pdf_page( plot_aerosol_kahn_3b_diff_corr() )
save_pdf_page( plot_aerosol_ice_diff_corr() )
save_pdf_page( plot_aerosol_water_diff_corr() )
#save_pdf_page(  )

# Decide to plot TCCON data only if there are TCCON soundings available
plot_tccon_data = False
if len(analysis_env.addl_objs) > 0:
    addl_ids = analysis_env.addl_objs[0].get_sounding_ids()[:]
    l2_match_indexes = analysis_env.correlate_soundings([analysis_env.comparison_ids, addl_ids])['data_id_indexes']

    if l2_match_indexes != None:
        plot_tccon_data = True

# Add some TCCON bar charts
if plot_tccon_data:
    save_pdf_page( plot_tccon_xco2_bar_mean() )

    ax = plot_tccon_xco2_bar_lfit_slope(sigma_filter=2.5)
    ax.set_title(ax.get_title() + "\nFiltered for difference <= 2.5 * mean difference")
    save_pdf_page(ax)

    for ax in plot_tccon_xco2_diff_corr():
        ax.set_title("TCCON " + ax.get_title())
        save_pdf_page(ax)

# Close reference bar chart before making TCCON specific files below
pdf_out_obj.close()

print "Finished plotting"
