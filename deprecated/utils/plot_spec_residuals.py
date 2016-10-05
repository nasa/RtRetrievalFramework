#!/usr/bin/env python

# Load standard modules
import os
import sys
import bisect
from optparse import OptionParser

# Load L2 modules
import L2_Version_Util
from L2_Plot_Util import plot_spec_residuals, combine_pdfs, VALID_GRID_NAMES, DEFAULT_DATA_COLUMN

def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [measured_spec_file] [simulated_spec_file]...")

    parser.add_option( "-d", "--detail", dest="detail",
                       metavar='NUM',
                       type="int",
                       default=0,
                       help="size in grid units of detail plots")

    parser.add_option( "-c", "--data_column", dest="data_column",
                       default=DEFAULT_DATA_COLUMN,
                       help="data column to be used from files being compared")

    parser.add_option( "-y", "--resid_ylims", dest="resid_ylims",
                       default=None,
                       help="vertical limits on residual plot (-percent, percent)")

    parser.add_option( "-g", "--grid", dest="grid",
                       default=VALID_GRID_NAMES[0], 
                       help="grid type to use for x axis. valid values: %s" % VALID_GRID_NAMES)

    parser.add_option( "--no_force", dest="no_force_grid",
                       default=True,
                       action='store_false',
                       help="do not force measured and simulation grids to be equal, ie do not compare by pixel")

    parser.add_option( "-o", "--out_file", dest="out_file",
                       metavar='FILE',
                       help="output filename")

    parser.add_option( "-t", "--title", dest="title",
                       metavar='TEXT',
                       help="additional title text")

    parser.add_option( "-n", "--notes", dest="notes",
                       metavar='TEXT',
                       action="append",
                       help="additional notes for bottom of first plot, option can be specified multiple times")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Need measured and simulated spectra files specified')

    meas_spec = args[0]
    sim_spec  = args[1]

    if options.out_file == None:
        meas_name = os.path.splitext(os.path.basename(meas_spec))[0]
        sim_name  = os.path.splitext(os.path.basename(sim_spec))[0]

        options.out_file = '%s-%s.pdf' % (meas_name, sim_name)

        print 'Writing output to: %s' % options.out_file

    if options.resid_ylims != None:
       options.resid_ylims = [float(f) for f in options.resid_ylims.split(",")] 

    bottom_notes = [os.path.realpath(meas_spec),
                    os.path.realpath(sim_spec),
                    '%s version: %s' % (os.path.basename(sys.argv[0]), L2_Version_Util.get_svn_version(sys.argv[0])),
                    ]

    if options.notes != None and len(options.notes) > 0:
        bottom_notes += options.notes
        
    plot_files = plot_spec_residuals( (meas_spec, sim_spec), 
                                      plot_title=options.title,
                                      bottom_notes=bottom_notes,
                                      grid_name=options.grid,
                                      detail_size=options.detail,
                                      data_column=options.data_column,
                                      force_grid_match=options.no_force_grid,
                                      resid_ylims = options.resid_ylims)

    if len(plot_files) != 0:
        combine_pdfs(plot_files, options.out_file, remove_source_files=True)

if __name__ == "__main__":
    standalone_main()


