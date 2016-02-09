from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
# General
import os
import sys
import copy
import bisect
import logging
import tempfile
import itertools
import subprocess

# Load additional modules
import numpy
import matplotlib.pyplot as pyplot

# L2 modules
from .oco_matrix import OcoMatrix

# General
OUTPUT_FORMAT = 'pdf'

# For combining pdfs
GHOST_SCRIPT_BIN = 'gs'

# For spectral residuals
VALID_GRID_NAMES = ['Wavenumber', 'Wavelength']
DEFAULT_DATA_COLUMN = 'Radiance'
NOISE_COLUMN = 'Error'
MEAS_SPEC_PLOT_OPT = 'black'
SIM_SPEC_PLOT_OPT  = 'red'
RES_SPEC_PLOT_OPT  = 'black'

def combine_pdfs(source_pdf_files, output_pdf_filename, remove_source_files=False):

    gs_cmd_arr = [GHOST_SCRIPT_BIN,
                  '-dBATCH', '-dNOPAUSE', '-q',
                  '-sDEVICE=pdfwrite',
                  '-sOutputFile=%s' % output_pdf_filename,
                  ] + source_pdf_files

    gs_process = None
    try:
        gs_process = subprocess.Popen(' '.join(gs_cmd_arr), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        gs_process.wait()
    except OSError:
        if gs_process == None:
            print('%s exited before subprocess could even do its work!' % GHOST_SCRIPT_BIN, file=sys.stderr)
        elif gs_process.returncode == None:
            # Premature exit due to error of inputs to gs
            # possibly just truncated file
            try:
                stdout_lines = gs_process.stdout.readlines()
            except IOError:
                stdout_lines = []
            
            if len(stdout_lines) > 0:
                print(stdout_lines, file=sys.stderr)
                raise IOError('Unexpected exit of %s' % GHOST_SCRIPT_BIN)
        elif gs_process.returncode != 0:
            raise OSError('Unexpected return of %s with return code: %d' % (GHOST_SCRIPT_BIN, gs_process.returncode))

    if remove_source_files:
        for tmp_file in source_pdf_files:
            os.remove(tmp_file)

def plot_jacobians(jacobian_files, group_name=None, label_names=None, figure_idx=1, page_indexes=None, **kwargs):
    """Routine for plotting groups of jacobian files.
    Returns created filenames.
    """

    logger = logging.getLogger()

    if group_name == None:
        group_name = '%s' % (os.path.commonprefix(jacobian_files),)

    if label_names == None:
        label_names = [ os.path.basename(curr_file) for curr_file in jacobian_files ]

    jacobian_objs = [OcoMatrix(curr_file) for curr_file in jacobian_files]

    common_pixels = jacobian_objs[0].pixels
    common_dims   = jacobian_objs[0].dims
    pixel_ranges = list(zip(common_pixels, common_pixels[1:] + [common_dims[0]]))

    logger.debug('common pixels: %s' % str(common_pixels))
    logger.debug('common dimensions: %s' % str(common_dims))

    for curr_obj in copy.copy(jacobian_objs):
        if list(curr_obj.pixels) != list(common_pixels) or list(curr_obj.dims) != list(common_dims):
            jacobian_objs.remove(curr_obj)
            print('Ignoring %s with inconsistent dimensions' % curr_obj.filename)

    if len(jacobian_objs) > 1:
        residual_combs = list(itertools.combinations(enumerate(jacobian_objs), 2))
        num_sep_subplots = 2
    else:
        residual_combs = None
        num_sep_subplots = 1

    if page_indexes == None:
        level_indexes = list(range(common_dims[1]))
    else:
        level_indexes = page_indexes

    output_filenames = []
    for level in level_indexes:
        logger.debug('Level %d, Figure Index: %d' % (level, figure_idx))
        fig_jac = pyplot.figure(figure_idx)
        fig_jac.clf()

        if num_sep_subplots > 1:
            fig_jac.set_size_inches(11.0*len(pixel_ranges), 8.5*num_sep_subplots)
            
            # Remove padding on left and right of plot
            fig_jac.subplots_adjust(left=0.04, right=0.98)
        else:
            fig_jac.set_size_inches(8.5*num_sep_subplots, 11.0)
            fig_jac.subplots_adjust()
            
        fig_jac.suptitle('%s, Level: %d' % (group_name, level), fontsize=14, fontweight='bold')

        plot_id = 1
        for range_idx, curr_range in enumerate(pixel_ranges):
            if num_sep_subplots > 1:
                curr_subplot = pyplot.subplot(num_sep_subplots, len(pixel_ranges), plot_id)
                curr_subplot.set_title('Over Plot, Window: %d' % range_idx)
            else:
                curr_subplot = pyplot.subplot(len(pixel_ranges), num_sep_subplots, plot_id)
                curr_subplot.set_title('Window: %d' % range_idx)
            
            plot_objs = []
            for file_idx, file_obj in enumerate(jacobian_objs):
                curr_plot = curr_subplot.plot(list(range(curr_range[0], curr_range[1])), file_obj.data[curr_range[0]:curr_range[1], level])
                plot_objs.append(curr_plot)

            # Build legend from first window
            if range_idx == len(pixel_ranges)-1 and num_sep_subplots > 1:
                legend = fig_jac.legend( plot_objs,
                                         label_names,
                                         loc='upper right',
                                         ncol=len(plot_objs),
                                         )
            
            plot_id += 1

        if residual_combs != None and len(residual_combs) > 0:
            for range_idx, curr_range in enumerate(pixel_ranges):
                curr_subplot = pyplot.subplot(num_sep_subplots, len(pixel_ranges), plot_id)
                curr_subplot.set_title('Residual. Window: %d' % range_idx)

                plot_objs = []
                res_labels = []
                for comb_idx, curr_comb in enumerate(residual_combs):
                    res_labels.append( '%s -\n %s' % (label_names[ curr_comb[0][0] ], label_names[ curr_comb[1][0] ]))
                    residual = curr_comb[0][1].data[curr_range[0]:curr_range[1], level] - curr_comb[1][1].data[curr_range[0]:curr_range[1], level]

                    curr_plot = curr_subplot.plot(list(range(curr_range[0], curr_range[1])), residual)
                    plot_objs.append(curr_plot)

                if range_idx == 1:
                    legend = fig_jac.legend( plot_objs,
                                             res_labels,
                                             loc='lower right',
                                             ncol=len(plot_objs),
                                             )

                plot_id += 1                                             

        with tempfile.NamedTemporaryFile(delete=False) as temp_fobj:
            pyplot.savefig(temp_fobj, format=OUTPUT_FORMAT)
            logger.debug('Saved plot to temporary file: %s' % temp_fobj.name)
            output_filenames.append(temp_fobj.name)
        figure_idx += 1

    return output_filenames


def plot_spec_figure(plot_title, meas_grid, sim_grid, grid_name, meas_data, meas_noise, sim_data, data_name, legend_labels, figure_index=1, bottom_notes=None, force_grid_match=True, rasterize_residual=False, resid_ylims=None):

    # If forcing a grid match make sure that both grids are same size
    if force_grid_match:
        if len(meas_data) != len(sim_data):
            raise Exception('Forcing grid match but data input does not have equal length')

        # Make sure grids are matching
        meas_grid       = sim_grid
        resid_grid      = sim_grid
        noise_for_resid = meas_noise
        meas_for_resid  = meas_data
        sim_for_resid   = sim_data
    else:
        # Find location where measurement grid and simulation grid overlap for residual
        meas_in_sim_idx = bisect.bisect_left(sim_grid, meas_grid[0])
        sim_in_meas_idx = bisect.bisect_left(meas_grid, sim_grid[0])

        if meas_in_sim_idx > sim_in_meas_idx:
            sim_beg  = meas_in_sim_idx
            sim_end  = bisect.bisect_left(sim_grid, meas_grid[-1])

            meas_beg = 0
            meas_end = sim_end-sim_beg

            resid_grid = sim_grid[sim_beg:sim_end]
        else:
            meas_beg = sim_in_meas_idx
            meas_end = bisect.bisect_left(meas_grid, sim_grid[-1])

            sim_beg  = 0
            sim_end  = meas_end-meas_beg

            resid_grid = meas_grid[meas_beg:meas_end]

        meas_for_resid  = meas_data[meas_beg:meas_end]
        sim_for_resid   = sim_data[sim_beg:sim_end]
        if meas_noise != None:
            noise_for_resid = meas_noise[meas_beg:meas_end]

    # Calculate spectral residual data
    if meas_noise != None:
        residual = old_div((meas_for_resid - sim_for_resid), noise_for_resid)
        chi2     = old_div(numpy.sum( numpy.power( old_div((meas_for_resid - sim_for_resid), noise_for_resid), 2 ) ), len(meas_for_resid))
    else:
        residual = 100.0 * (meas_for_resid - sim_for_resid) / max(meas_for_resid)
        chi2     = None

    # Determine xlimit based on greatest extent of both sets of data
    xextents = [ min(min(meas_grid), min(sim_grid)), max(max(meas_grid), max(sim_grid)) ]

    # Set where to plot
    fig = pyplot.figure(figure_index)
    fig.clf()
    
    fig.set_size_inches(8.5,11.0)
    fig.subplots_adjust()
    
    fig.suptitle(plot_title, fontsize=14, fontweight='bold')

    # Plot residual
    axis_top = pyplot.subplot(211)
    pos_top = axis_top.get_position()
    points_top = pos_top.get_points()
    points_top[0][1] += 0.2
    pos_top.set_points(points_top)
    axis_top.set_position(pos_top)
    
    # we provide the option of rasterizing the residual, preventing the problem where
    # thin features (like spikey residuals) disappear from the overview.  This
    # can increase computation time and file size.
    pyplot.plot(resid_grid, residual, RES_SPEC_PLOT_OPT, rasterized=rasterize_residual)
    pyplot.xlim( xextents ) # make sure axis in increasing order

    if resid_ylims != None:
        interval = (resid_ylims[0], resid_ylims[1])
        print("setting residual limits to[ [%f, %f]" % interval)
        pyplot.ylim(interval)
        axis_top.set_ylabel('Residual (%)')
    elif meas_noise != None:
        axis_top.set_ylabel('(Meas - Model) / Noise')
        residual_title = r'$\chi^2$ = %f' % chi2

        where_same = numpy.where(meas_noise == meas_noise[0])

        if len(where_same[0]) == len(meas_noise):
            residual_title += ', noise level = %e' % meas_noise[0]
            
        axis_top.set_title(residual_title)
    else:
        axis_top.set_ylabel('Residual (%)')
    
    # Plot overplot of measured and simulated data
    axis_bottom = pyplot.subplot(212)
    pos_bottom = axis_bottom.get_position()
    points_bottom = pos_bottom.get_points()

    points_bottom[1][1] += 0.2

    if bottom_notes != None:
        if len(bottom_notes) > 5:
            # we'll need some more room
            points_bottom[0][1] += 0.01*len(bottom_notes)
        else:
            points_bottom[0][1] += 0.002*len(bottom_notes)
    
    pos_bottom.set_points(points_bottom)
    axis_bottom.set_position(pos_bottom)

    plot_meas_spec = pyplot.plot(meas_grid, meas_data, MEAS_SPEC_PLOT_OPT)
    plot_sim_spec  = pyplot.plot(sim_grid, sim_data, SIM_SPEC_PLOT_OPT)
    pyplot.xlim( xextents ) # make axis consistent with residual
    pyplot.ylabel(data_name)
    pyplot.xlabel(grid_name)

    leg = axis_bottom.legend( (plot_meas_spec, plot_sim_spec),
                              legend_labels,
                              loc='lower right',
                              borderpad = 0.1,
                              labelspacing = 0.05,
                              borderaxespad = 0.2,
                              )
    leg.prop.set_size(8)

    if bottom_notes:
        text = pyplot.figtext(0.02, 0.02, '\n'.join(bottom_notes),
                              size='x-small')

    tmp_filename = None
    with tempfile.NamedTemporaryFile(delete=False) as temp_fobj:
        pyplot.savefig(temp_fobj, format=OUTPUT_FORMAT)
        tmp_filename = temp_fobj.name

    return tmp_filename

def plot_spec_residuals(spectra_files, plot_title=None, spectra_labels=None, bottom_notes=None, grid_name=VALID_GRID_NAMES[0], detail_size=0, force_grid_match=True, data_column=DEFAULT_DATA_COLUMN,  resid_ylims=None, **kwargs):

    if len(spectra_files) != 2:
        raise Exception('Can only handle plotting of two spectra files not %d files as specified: %s' % (len(spectra_files), spectra_files))

    meas_spec_file, sim_spec_file = spectra_files
    
    meas_obj = OcoMatrix(meas_spec_file)
    sim_obj  = OcoMatrix(sim_spec_file)

    # Check that Start_Pixels in two files are of the same length, ie same number of bands
    if len(meas_obj.pixels) != len(sim_obj.pixels):
        raise Exception('Measured spectra file %s and simlated spectra file %s have different number of bands represented by the start pixel index header keyword' % (meas_spec_file, sim_spec_file))

    # Check that both files have same measurement grid
    if grid_name.lower() in meas_obj.labels_lower:
        meas_grid_col_idx = meas_obj.labels_lower.index(grid_name.lower())
    else:
        meas_grid_col_idx = None
    if grid_name.lower() in sim_obj.labels_lower:
        sim_grid_col_idx = sim_obj.labels_lower.index(grid_name.lower())
    else:
        sim_grid_col_idx = None

    # Column indexes for data
    meas_data_col_idx = meas_obj.labels_lower.index(data_column.lower())
    sim_data_col_idx  = sim_obj.labels_lower.index(data_column.lower())

    if NOISE_COLUMN.lower() in meas_obj.labels_lower:
        noise_col_idx = meas_obj.labels_lower.index(NOISE_COLUMN.lower())
    else:
        noise_col_idx = None
        

    # Set up indexing into data of windows
    meas_start_indexes = meas_obj.pixels[:]
    meas_end_indexes   = [ pix for pix in meas_obj.pixels[1:] ]
    meas_end_indexes.append(meas_obj.dims[0])

    sim_start_indexes = sim_obj.pixels[:]
    sim_end_indexes   = [ pix for pix in sim_obj.pixels[1:] ]
    sim_end_indexes.append(sim_obj.dims[0])

    win_indexes = list(zip(meas_start_indexes, meas_end_indexes, sim_start_indexes, sim_end_indexes))

    if len(win_indexes) == 0:
        #raise Exception('Could not find any window indexes')
        win_indexes = ( (0, meas_obj.dims[0], 0, sim_obj.dims[0]), )

    tmp_file = None
    tmp_out_files = []
    figure_index = 1
    for window_idx, indexes in enumerate(win_indexes):
        meas_beg = indexes[0]
        meas_end = indexes[1]
        sim_beg  = indexes[2]
        sim_end  = indexes[3]
        
        # Save grid values for use in plotting
        if meas_grid_col_idx != None:
            meas_grid = meas_obj.data[meas_beg:meas_end, meas_grid_col_idx]
        else:
            meas_grid = numpy.arange(meas_beg,meas_end)

        if sim_grid_col_idx != None:
            sim_grid = sim_obj.data[sim_beg:sim_end, sim_grid_col_idx]
        else:
            sim_grid = numpy.arange(sim_beg,sim_end)

        # Get data values into local variables
        meas_data = meas_obj.data[meas_beg:meas_end, meas_data_col_idx]
        sim_data  = sim_obj.data[sim_beg:sim_end, sim_data_col_idx]

        if noise_col_idx != None:
            meas_noise = meas_obj.data[meas_beg:meas_end, noise_col_idx]
        else:
            meas_noise = None

        if plot_title != None:
            win_title = plot_title + '\n'
        else:
            win_title = ''

        if spectra_labels == None:
            spectra_labels = ('Measured', 'Simulated')
            
        win_title += 'Spectrum Residuals - Window %d' % (window_idx+1)

        # Make plot of desired data
        tmp_file = plot_spec_figure(win_title,
                                    meas_grid,
                                    sim_grid,
                                    grid_name,
                                    meas_data,
                                    meas_noise,
                                    sim_data,
                                    data_column,
                                    spectra_labels,
                                    bottom_notes=bottom_notes,
                                    figure_index=figure_index,
                                    force_grid_match=force_grid_match,
                                    rasterize_residual=True,
                                    resid_ylims=resid_ylims,
                                    )
        tmp_out_files.append(tmp_file)
        figure_index += 1

        if detail_size > 0:
            win_size =  max(meas_grid) - min(meas_grid)

            if detail_size > win_size:
                raise Exception('detail size: %f greater than window %d size: %f' % (detail_size, window_idx+1, win_size))

            num_detail = int(old_div(win_size, detail_size))
            num_win_points = int(old_div(len(meas_grid), num_detail))

            if num_detail < round(old_div(win_size, detail_size)):
                num_detail += 1

            beg_idx = 0
            end_idx = 0
            for detail_idx in range(num_detail):
                beg_idx = end_idx
                end_idx = min(beg_idx+num_win_points, len(meas_grid))

                # Subsumed all but single last point
                if (end_idx - beg_idx) == 0:
                    break

                if meas_noise != None:
                    detail_noise = meas_noise[beg_idx:end_idx]
                else:
                    detail_noise = None

                tmp_file = plot_spec_figure(win_title + ', Detail %d' % (detail_idx+1),
                                            meas_grid[beg_idx:end_idx],
                                            sim_grid[beg_idx:end_idx],
                                            grid_name,
                                            meas_data[beg_idx:end_idx],
                                            detail_noise,
                                            sim_data[beg_idx:end_idx],
                                            data_column,
                                            spectra_labels,
                                            figure_index=figure_index,
                                            force_grid_match=force_grid_match,
                                            rasterize_residual=False # keep file size low
                                            )
                tmp_out_files.append(tmp_file)
                figure_index += 1

    return tmp_out_files
