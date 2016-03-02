#!/usr/bin/env python

import os
from argparse import ArgumentParser
import logging
import itertools

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

from full_physics.log_util import init_logging
from full_physics import l2_analysis

def standalone_main():
    parser = ArgumentParser()

    parser.add_argument("l2_filename",
            help="L2 file from which to plot soundings")
    parser.add_argument("-i", "--sounding_id", metavar="SOUNDING_ID", action="append",
            help="Plot the single specified sounding id")
    parser.add_argument("-l", "--sounding_id_list", metavar="FILE", 
            help="A file with a list of sounding ids to plot")
    parser.add_argument("-o", "--output_file", metavar="FILE",
            help="Name of output file to be used instead of default")
    parser.add_argument("-v", "--verbose", action="store_true",
            help="Output additional verbose information")

    # Parse command line arguments
    args = parser.parse_args()

    if (args.verbose):
        init_logging(logging.DEBUG)
    else:
        init_logging(logging.INFO)

    if not os.path.exists(args.l2_filename):
        raise IOError("L2 filename specified does not exist: %s" % args.l2_filename)
    
    analysis = l2_analysis.AnalysisEnvironment([args.l2_filename])
    analysis.add_routines(l2_analysis.FrankenburgPlotRoutines(analysis_env=analysis))

    plot_sounding_ids = []
    if args.sounding_id or args.sounding_id_list:
        id_type = analysis.comparison_ids[0].dtype.type
        if args.sounding_id:
            plot_sounding_ids += [ id_type(snd_id) for snd_id in args.sounding_id ]
        if args.sounding_id_list:
            with open(args.sounding_id_list) as id_file:
                plot_sounding_ids += [ id_type(line.strip()) for line in id_file ]
    else:
        plot_sounding_ids = analysis.comparison_ids

    if args.output_file:
        output_file = args.output_file
    else:
        if len(plot_sounding_ids) > 1:
            output_file = "spectral_fit_%s-%s.pdf" % (plot_sounding_ids[0], plot_sounding_ids[-1])
        else:
            output_file = "spectral_fit_%s.pdf" % (plot_sounding_ids[0])

    logging.info("Writing plots to: %s" % output_file)
    pdf = PdfPages(output_file)

    num_soundings = len(plot_sounding_ids)

    for idx, (snd_id, figs) in enumerate(zip(plot_sounding_ids, analysis.call_analysis_routine("plot_spectral_fits", ids=(plot_sounding_ids,), all_soundings=True))):
        logging.info("[%d / %d] %s" % (idx+1, num_soundings, snd_id))
        for band_fig in figs:
            pdf.savefig(band_fig)
    pdf.close()
        
if __name__ == "__main__":
    standalone_main()
