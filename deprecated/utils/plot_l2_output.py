#!/usr/bin/env python

import os
import pyx
import sys
import glob
from optparse import OptionParser

from OCO_Matrix import OCO_Matrix
from OCO_MathUtil import *
import OCO_Statistics
import OCO_TextUtils

def plot(matrix, abscissa, abscissa_string, out_filename=None, det_windows=True):
    
    if (matrix.dims[0] == 0 or matrix.dims[1] == 0): return
    
    letter = pyx.document.paperformat.Letter
    
    # Number of spectral points
    size = matrix.dims[0]
    
    # All pages will be stored in this object
    d = pyx.document.document(pages=[])
    
    # starting pixel indices for each window.  The last entry is
    # the total number of pixels.
    windows = matrix.pixels[:] # Make a copy instead of using a reference
    if (len(windows) == 0):
        if det_windows and matrix.filename.find("high_res") == -1:
            windows = [ wini * (matrix.dims[0] / 3) for wini in range(3) ]
        else:
            windows = [0]

    windows.append(matrix.dims[0])

    # Run through the columns
    for index in range(0, matrix.dims[1]):
        # Don't bother plotting wavelength
        try:
            if (matrix.labels_lower.index("wavelength") == index): continue
        except ValueError:
            pass
        
        try:
            xindex = matrix.labels_lower.index(abscissa)
            if (xindex == index): continue
        except ValueError:
            xindex = -1
            pass

        # Plot each window separately
        for w in range(0, len(windows) - 1):
            
            # Look for the abscissa_string.  If present, use that as
            # the x axis.
            if (xindex != -1):
                xaxis = matrix.data[windows[w]:windows[w+1], xindex]
                xtitle = abscissa_string
            else:
                xaxis = range(0, windows[w+1] - windows[w])
                xtitle = "Index"
                
            # Build the list object holding the data for this window.
            list = [(xaxis[0], matrix.data[windows[w], index])]
            for i in range(0, windows[w+1] - windows[w]):
                list.append((xaxis[i], matrix.data[i + windows[w], index]))
                
            # use the list to create the data object to plot
            title=OCO_TextUtils.escape_underscore(os.path.dirname(matrix.filename))
            data = pyx.graph.data.points(list, title=title, x=1, y=2)

            # Assign the axis ranges.
            xmin = xaxis[0]
            xmax = xaxis[windows[w+1] - windows[w] - 1]
            ymax = None
            ymin = None
            if (min(matrix.data[windows[w]:windows[w+1], index]) ==
                max(matrix.data[windows[w]:windows[w+1], index])):
                ymax = 1
                ymin = -1

            # create the graph object
            g = pyx.graph.graphxy(width=100,
                                  key=pyx.graph.key.key(pos="br",
                                                        symbolwidth=5*pyx.unit.t_cm),
                                  x=pyx.graph.axis.linear(min=xmin,
                                                          max=xmax,
                                                          title=xtitle),
                                  y=pyx.graph.axis.linear(min=ymin,
                                                          max=ymax,
                                                          title=matrix.ytitle[index]))
            # plot the data
            g.plot(data, [pyx.graph.style.line([pyx.color.rgb.red,
                                                pyx.style.linewidth.THICK])])
            # Add the title
            filename = OCO_TextUtils.escape_underscore(matrix.filename)
            title = "%s: Index %d, window %d" % (filename, index, w)
            g.text(g.width/2, g.height + 0.4, title,
                   [pyx.text.halign.center, pyx.text.valign.bottom,
                    pyx.text.size.Large])
            
            # Add this plot to the document
            d.append(pyx.document.page(g,
                                       paperformat=letter, 
                                       margin=3*pyx.unit.t_cm,
                                       fittosize=1,
                                       rotated=1, centered=1))

    # Write the postscript file
    if out_filename == None:
        out_filename = "%s.pdf" % matrix.filename
        
    print 'Writing: %s' % out_filename
    d.writetofile(out_filename)

            
def plotRMS(matrix, ref, abscissa, abscissa_string, out_filename=None, pressure_matrix=None, det_windows=True, ratio=False):
    
    if (matrix.dims[0] == 0 or matrix.dims[1] == 0): return
    if (ref.dims[0] == 0 or ref.dims[1] == 0): return
    
    # Number of spectral points
    size = matrix.dims[0]
    
    # All pages will be stored in this object
    d = pyx.document.document(pages=[])
    
    # starting pixel indices for each window.  The last entry is
    # the total number of pixels.
    windows = matrix.pixels[:] # Make a copy instead of using a reference
    if (len(windows) == 0):
        if det_windows and matrix.filename.find("high_res") == -1:
            windows = [ wini * (matrix.dims[0] / 3) for wini in range(3) ]
        else:
            windows = [0]

    windows.append(matrix.dims[0])

    # Get pressure levels to be used for partial derivate files to name the indexes
    if pressure_matrix != None:
        pressure_list = pressure_matrix.data[:, 0]
    else:
        pressure_list = []

    # Run through the columns
    for index in range(0, min(ref.dims[1], matrix.dims[1])):

        # Don't bother plotting wavelength or wavenumber
        if len(matrix.labels) > 0:
            if (matrix.labels[index].lower() == "wavelength"): continue
            if (matrix.labels[index].lower() == "wavenumber"): continue

        try:
            xindex = matrix.labels_lower.index(abscissa)
            if (xindex == index): pass # continue
        except ValueError:
            xindex = -1
            pass
        
        # This is the canvas that holds both graphs
        c = pyx.canvas.canvas()

        curr_ypos = 0

        # Plot each window separately
        num_windows = len(windows) - 1;
        for winEndIdx in range(num_windows, 0, -1):

            winBegIdx = winEndIdx-1

            # Look for the abscissa_string.  If present, use that as
            # the x axis.
            if (xindex != -1):
                xaxis = matrix.data[windows[winBegIdx]:windows[winEndIdx], xindex]
                xtitle = abscissa_string
            else:
                xaxis = range(0, windows[winEndIdx] - windows[winBegIdx])
                xtitle = "Index"


            selfdata = matrix.data[windows[winBegIdx]:windows[winEndIdx], index]
            refdata = ref.data[windows[winBegIdx]:windows[winEndIdx], index]

            if ratio:
                residual, rms = OCO_Statistics.Ratio_RMS(selfdata, refdata)
            else:
                residual, rms = OCO_Statistics.Residual_RMS(selfdata, refdata)

            if (rms < 0):
                print "residual has NaN's"
                return
            
            # Set up strings used to plot and for verbose output on screen
            if len(matrix.labels) > index:
                index_name = matrix.labels[index]
            elif len(pressure_list) > 0 and matrix.filename.find('_pd') > 0:
                #print 'pll:', len(pressure_list), 'index:', index
                # Use pressure level for label
                index_name = "Pressure %.2f" % pressure_list[index]
            else:
                index_name = "Index %d" % index

            if num_windows > 1:
                plot_description = "%s: %s, window %d"  % (matrix.filename, index_name, winBegIdx)
            else:
                plot_description = "%s: %s"  % (matrix.filename, index_name)
            
            if ratio:
                rms_description = "Ratio RMS = %f" % (rms)
            else:
                rms_description = "RMS = %f%%" % (rms * 100)
            
            print "%s, %s" % (plot_description, rms_description)

            # build the lists to plot
            list0 = [(xaxis[0], selfdata[0])]
            list1 = [(xaxis[0], refdata[0])]
            if ratio:
                diff = [(xaxis[0], residual[0])]
            else:
                diff = [(xaxis[0], residual[0]*100.)]
            for i in range(0, windows[winEndIdx] - windows[winBegIdx]):
                list0.append((xaxis[i], selfdata[i]))
                list1.append((xaxis[i], refdata[i]))
                if ratio:
                    diff.append((xaxis[i], residual[i]))
                else:
                    diff.append((xaxis[i], residual[i]*100.))
                
            # Create the PyX data structures
            title = OCO_TextUtils.escape_underscore(os.path.dirname(matrix.filename))
            data0 = pyx.graph.data.points(list0, title=title, x=1, y=2)
            title = OCO_TextUtils.escape_underscore(os.path.dirname(ref.filename))
            if (ref.filename != matrix.filename):
                title = "%s%s" % (title, OCO_TextUtils.escape_underscore(ref.filename))
            data1 = pyx.graph.data.points(list1, title=title, x=1, y=2)
            data2 = pyx.graph.data.points(diff, x=1, y=2)
            
            # Assign the axis ranges.
            xmin = xaxis[0]
            xmax = xaxis[windows[winEndIdx] - windows[winBegIdx] - 1]
            ymax = None
            ymin = None

            if (min(residual) == max(residual)):
                ymax = 1
                ymin = -1
                
            # This is the lower graph
            if ratio:
                y_title_txt = "ratio"
            else:
                y_title_txt = "difference"
            g1 = c.insert(pyx.graph.graphxy(width=100, ratio=1.5*num_windows,
                                            ypos=curr_ypos,
                                            x=pyx.graph.axis.linear(min=xmin, max=xmax,
                                                                    title=xtitle),
                                            y=pyx.graph.axis.linear(min=ymin, max=ymax,
                                                                    title=y_title_txt)))
            # Plot the differences
            g1.plot(data2,
                    [pyx.graph.style.line([pyx.color.rgb.red,
                                           pyx.style.linewidth.THICK])])
           
            # Add the title
            if ratio:
                title="Ratio RMS = %f" % (rms)
            else:
                title="Residual RMS = %f\%%" % (rms * 100)
            
            g1.text(g1.width/2, curr_ypos + g1.height + 0.4, title,
                    [pyx.text.halign.center, pyx.text.valign.bottom])

            curr_ypos += g1.height+6.0
            
            # Find the y axis range
            ymax = None
            ymin = None
            if (min(selfdata) == max(selfdata)
                and min(refdata) == max(refdata)
                and min(refdata) == min(selfdata)):
                ymax = 1
                ymin = -1
                
            # This is the upper graph
            g0 = c.insert(pyx.graph.graphxy(width=100, ratio=1.5*num_windows,
                                            ypos=curr_ypos,
                                            key=pyx.graph.key.key(pos="br", symbolwidth=5*pyx.unit.t_cm),
                                            x=pyx.graph.axis.linkedaxis(g1.axes["x"]),
                                            y=pyx.graph.axis.linear(min=ymin, max=ymax, title=matrix.ytitle[index])))           
            # Now plot the values
            g0.plot([data0, data1],
                    [pyx.graph.style.line([pyx.color.gradient.Rainbow,
                                           pyx.style.linewidth.THICK])])
            # Add the title
            title = OCO_TextUtils.escape_underscore(plot_description)
            g0.text(g0.width/2, curr_ypos + g0.height + 0.8, title,
                    [pyx.text.halign.center, pyx.text.valign.bottom])

            curr_ypos += g0.height+15.0
            
        # Add this plot to the document
        d.append(pyx.document.page(c,
                                   paperformat=pyx.document.paperformat.Letter,
                                   margin=pyx.unit.t_cm,
                                   fittosize=1,
                                   rotated=0, centered=1))
            
    # Write the postscript file
    if out_filename == None:
        out_filename = "%s_RMS.pdf" % matrix.filename

    print 'Writing: %s' % out_filename
    d.writetofile(out_filename)

if __name__ == "__main__":

    help_text = "To plot all files in the specified directory or compare two directories:\n \
\t %prog runDir [refDir]\n \
e.g.\t %prog test01/A/run [test01/A/std_output]\n\n \
To plot a particular file or compare two files:\n \
\t %prog runFile [refFile]\n \
e.g.\t %prog test01/A/run/out/press_pd.dat [test01/A/std_output/out/press_pd.dat]"

    parser = OptionParser(usage="usage: %prog [options] run_dir_or_file [ref_dir_or_file]\n" + help_text)

    parser.add_option( "-o", "--out_filename", dest="out_filename",
                       metavar="FILE",
                       help="name of file where plot is written, only used when plotting individual files, not directories",
                       )
    
    parser.add_option( "-i", "--iter", dest="iteration",
                       metavar="NUM",
                       type="int",
                       help="iteration number for appending to comparison files",
                       )

    parser.add_option( "--no_win", dest="det_windows",
                       action="store_false", default=True,
                       help="do not auto determine window ranges",
                       )

    parser.add_option( "--ratio", dest="ratio",
                       action="store_true", default=False,
                       help="use ratio comparison instead of residual",
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error('Need to specify at least one directory or filename')

    run_directory = args[0]

    if len(args) > 1:
        ref_directory = args[1]
    else:
        ref_directory = None

    if (len(args) < 1):
        parser.error("No run directory or filename specified")
    
    pyx.unit.set(xscale=8)
    pyx.text.set(mode="latex")
    pyx.text.preamble(r"\usepackage{times}")


    if (os.path.isdir(args[0])):
        comp_globs = ['rad_meas.dat',
                      'rad_conv.dat',
                      'control1/t_pd.dat',
                      'control1/aer_pd_species_?.dat',
                      'control1/alb_pd.dat',
                      'control1/alb_i_pd.dat',
                      'control1/cont_pd.dat', 
                      'control1/disp_pd.dat', 
                      'control1/mr_pd_species_?.dat',
                      'control1/mr_scaling_pd_species_?.dat', 
                      'control1/press_pd.dat',
                      'high_res??.dat',
                      'high_res??_??.dat',
                      'solar_all??.dat',
                      ]

        for curr_glob in comp_globs:
            run_file_list = glob.glob("%s/out/%s" % (run_directory, curr_glob))
            ref_file_list = glob.glob("%s/out/%s" % (ref_directory, curr_glob))

            run_file_list.sort()
            ref_file_list.sort()

            if len(run_file_list) != len(ref_file_list):
                print >>sys.stderr, 'Did not find same number of files in %s and %s for glob %s' % (run_directory, ref_directory, curr_glob)
                    
            for run_filename, ref_filename in zip(run_file_list, ref_file_list):
                if not (os.path.basename(run_filename) == os.path.basename(ref_filename)):
                    raise IOError('Mismatched comparison of %s and %s' % (run_filename, ref_filename))
            
                if options.iteration != None:
                    iter_ext = ".iter%02d" % options.iteration

                    if os.path.exists(run_filename + iter_ext) and os.path.exists(ref_filename + iter_ext):
                        run_filename += iter_ext
                        ref_filename += iter_ext

                try:
                    print "reading %s" % os.path.basename(run_filename)
                    run = OCO_Matrix(run_filename)

                    # Create a difference plot if two directory names
                    # or create a single plot for the files in the directory
                    if (ref_directory != None):
                        ref = OCO_Matrix(ref_filename)
                        plotRMS(run, ref, "wavenumber", "Wavenumber (cm$^{-1}$)", det_windows=options.det_windows, ratio=options.ratio)
                    else:
                        plot(run, "wavenumber", "Wavenumber (cm$^{-1}$)", det_windows=options.det_windows)

                except IOError:
                    continue
                except:
                    print "Error plotting %s" % os.path.basename(run_filename)
    else:
        try:
            print "reading %s" % run_directory
            run = OCO_Matrix("%s" % run_directory)

            if (ref_directory != None):
                print "reading %s" % ref_directory
                ref = OCO_Matrix("%s" % ref_directory)
                plotRMS(run, ref, "wavenumber", "Wavenumber (cm$^{-1}$)", out_filename=options.out_filename, det_windows=options.det_windows, ratio=options.ratio)
            else:
                plot(run, "wavenumber", "Wavenumber (cm$^{-1}$)", out_filename=options.out_filename, det_windows=options.det_windows)
        except IOError:
            pass
