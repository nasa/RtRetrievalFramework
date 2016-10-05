#!/usr/bin/env python

import os
import sys
from regression_testing import *
from optparse import OptionParser

if os.getenv('L2_SUPPORT_PATH') == None:
    print 'L2_SUPPORT_PATH must be defined.'
    sys.exit(1)

l2_support_path = os.getenv('L2_SUPPORT_PATH')

class cr_config:
    reg_compare_bin_name        = l2_support_path + "/reg_compare/regression_compare.pl"
    reg_compare_cmd_prefix      = ""
    reg_compare_verbose         = False
    reg_compare_matrix_map_name = 'reg_comp_matrix_match_map.txt'

    plot_diff_bin_name          = l2_support_path + "/utils/plot_l2_output.py"

    do_diff_plotting            = False

    pressure_file_name          = 'pressure.dat'

    # List of files that should be plotted if there are differences
    # above tolerance level detected
    plot_comparison_files = ['out/rad_meas.dat',
                             'out/rad_conv.dat',
                             'out/control1/t_pd.dat',
                             'out/control1/aer_pd_species_1.dat',
                             'out/control1/alb_pd.dat',
                             'out/control1/alb_i_pd.dat',
                             'out/control1/cont_pd.dat', 
                             'out/control1/disp_pd.dat', 
                             'out/control1/mr_pd_species_1.dat',
                             'out/control1/mr_pd_species_3.dat',
                             'out/control1/mr_scaling_pd_species_1.dat', 
                             'out/control1/mr_scaling_pd_species_3.dat', 
                             'out/control1/press_pd.dat',
                             'out/high_res01.dat',
                             'out/high_res02.dat', 
                             'out/high_res03.dat',
                             'out/solar_all01.dat',
                             'out/solar_all02.dat', 
                             'out/solar_all03.dat' ]
    
class compare_results(regression_testing):

    def __init__ (self):
        pass

    def locate_file(self, base_name, location_file):
        
        file_location = ''
        if os.path.exists(location_file):
            file_loc_fobj = open(location_file, 'r')
            file_location = file_loc_fobj.readline()
            file_loc_fobj.close()
        else:
            print 'Looking for needed file: ' + base_name
            findCmd = 'find ~/ -name ' + base_name + ' -maxdepth 5'

            ff = os.popen(findCmd, 'r')
            file_location = ff.readline()
            ff.close()

        # Remove any newlines on locations string
        file_location = file_location.strip()

        if not os.path.exists(file_location):
            return None
        elif not os.path.exists(location_file):
            file_loc_fobj = open(location_file, 'w')
            file_loc_fobj.write(file_location)
            file_loc_fobj.close()

        return file_location.strip()

    def get_pressure_filename(self, std_dir):
        findCmd = 'find ' + std_dir + ' -name ' + cr_config.pressure_file_name + ' -maxdepth 5'
        
        ff = os.popen(findCmd, 'r')
        pressure_file = ff.readline()
        ff.close()

        return pressure_file.strip()

    def plot_test_diffs(self, plot_module, test_dir, std_dir, run_dir, plot_filenames):

        # Default to plotting all files if none were specified to function
        if plot_filenames == None:
            plot_filenames = cr_config.plot_comparison_files

        pressure_file = self.get_pressure_filename(std_dir)

        if os.path.exists(pressure_file):
            pressure_data = plot_module.OCO_Matrix(pressure_file)
        else:
            pressure_data = None

        for comp_file in plot_filenames:
            # Only plot files that were found to be difference
            if not comp_file in cr_config.plot_comparison_files:
                continue
        
            try:
                run_data = plot_module.OCO_Matrix("%s/%s" % (run_dir, comp_file))
                std_data = plot_module.OCO_Matrix("%s/%s" % (std_dir, comp_file))
            
                plot_filename = "%s/%s_%s-%s.pdf" % (test_dir, os.path.basename(run_dir), os.path.basename(std_dir), os.path.basename(comp_file))
                plot_module.plotRMS(run_data, std_data, "wavenumber", "Wavenumber (cm$^{-1}$)", out_filename=plot_filename, pressure_matrix=pressure_data)
            except IOError:
                continue
            except:
                print "Error plotting %s" % comp_file

    def compare_test_dir(self, comparison_dir_name, baseline_dir_name, plot_module, test_dir):
        baseline_dir = test_dir + '/' + baseline_dir_name       
        comparison_dir = test_dir + '/' + comparison_dir_name

        report_file = test_dir + "/" + comparison_dir_name + '_' + baseline_dir_name + "-file_diffs.log"

        if not os.path.exists(baseline_dir):
            print 'baseline directory does not exist for test: "%s" as "%s"' % (test_dir, baseline_dir)
            return

        if not os.path.exists(comparison_dir):
            print 'comparison directory does not exist for test: "%s" as "%s"' % (test_dir, comparison_dir)
            return
            
        print 'Comparing %s against %s for test:\n%s' % (os.path.basename(comparison_dir), os.path.basename(baseline_dir), test_dir)

        print '- Creating file comparisons'
        sys.stdout.flush()

        compare_map_file = '%s/%s' % (os.path.realpath(test_dir), cr_config.reg_compare_matrix_map_name)
        
        compare_cmd = cr_config.reg_compare_cmd_prefix + \
                      os.path.realpath(cr_config.reg_compare_bin_name) + \
                      " -baseline_dir " + os.path.realpath(baseline_dir) + "/out " + \
                      "-comparison_dir " + os.path.realpath(comparison_dir) + "/out " + \
                      "-all " + \
                      "-find_ignore 'pid' " + \
                      "-find_ignore '\.svn' " + \
                      "-find_ignore 'sv_names\.dat' " + \
                      "-report_file " + os.path.realpath(report_file) + \
                      " -module_options \"comparison_map_file=" + compare_map_file + \
                      ",delete_comparison_map=0\"" + \
                      " 2>&1"
        
        ff = os.popen(compare_cmd, 'r')

        for compare_line in ff.readlines():
            if cr_config.reg_compare_verbose:
                print compare_line.rstrip()

        ff.close()
        print 'Created log ' + report_file
        sys.stdout.flush()

        # Shell to use tail because it is fast...
        if os.path.exists(report_file):
            ff = os.popen('tail -n 1 ' + report_file, 'r')
            diff_count = ff.readline()
            ff.close()

            print diff_count.replace('There are ', '').strip()

        if os.path.exists(compare_map_file):
            map_fobj = open(compare_map_file, 'r')
            map_out_contents = map_fobj.readlines()
            map_fobj.close()
            os.remove(compare_map_file)
           
            plot_filenames = []
            for map_line in map_out_contents:
                map_parts = map_line.split()
                if map_parts[2].isdigit() and int(map_parts[2]) > 0:
                    base_filename = map_parts[0].replace(os.path.realpath(baseline_dir), '').lstrip('/')
                    plot_filenames.append(base_filename)

        else:
            # If the map file does not exist then do not present an empty list
            # meaning that all files are in agreement, instead make plot_filenames
            # list none meaning that the list could not be constructed
            plot_filenames = None
        
        if cr_config.do_diff_plotting:
            print '- Creating plot comparisons'
            sys.stdout.flush()
            self.plot_test_diffs(plot_module, test_dir, baseline_dir, comparison_dir, plot_filenames)

        print ''
        sys.stdout.flush()

    def main(self):

        # Set up command line arguments
        parser = OptionParser(usage = '%prog [options] <testcase_pattern> [<testcase_pattern> ...]')

        parser.add_option( "-c", "--comparison_directory", dest="comparison_dir_name",
                           metavar="DIR_NAME",
                           help="basename of directory to compare to baseline results.",
                           )

        parser.add_option( "-b", "--baseline_directory", dest="baseline_dir_name",
                           metavar="DIR_NAME",
                           help="basename of directory considered to be the correct results. Default ['%s']" % rt_config.std_output_dir_name,
                           default=rt_config.std_output_dir_name
                           )

        parser.add_option( "-t", "--test_base_dir", dest="test_base_dir",
                           metavar="DIR_PATH",
                           default=rt_config.default_test_base_dir,
                           help="base directory for finding testcases. Default ['%s']" % rt_config.default_test_base_dir
                           )

        parser.add_option( "-p", "--plot_differences", dest="diff_plotting",
                           action="store_true",
                           help="turn on plotting of differences of files",
                           default=cr_config.do_diff_plotting
                           )

        (options, args) = parser.parse_args()

        #####
        # Parse options related to the binary name used for testing
        if options.comparison_dir_name == None or len(options.comparison_dir_name) == 0:
            parser.error('A comparison directory must be specified witht he -c argument')

        print 'Comparing run directories named %s with baseline directories named %s' % (options.comparison_dir_name, options.baseline_dir_name)

        ####
        # Parse options related to the names of the testcases to be run
        run_test_dirs = self.parse_test_names(options.test_base_dir, args)

        if run_test_dirs == None:
            sys.exit(1)
        elif len(run_test_dirs) > 0:
            print ''
            print 'Comparing test directories:'
            for test_dir in run_test_dirs:
                print test_dir
        else:
            parser.error('No test directories found')

        plot_diff_path = os.path.dirname(cr_config.plot_diff_bin_name)
        if not plot_diff_path in sys.path:
            sys.path.append(plot_diff_path)

        # Import plotting module, get basename and remove .py extension
        # so that python can properly find the module
        if options.diff_plotting:
            plot_module = __import__(os.path.basename(cr_config.plot_diff_bin_name).replace('.py', ''), globals(), {}, [])

            plot_module.pyx.unit.set(xscale=8)
            plot_module.pyx.text.set(mode="latex")
            plot_module.pyx.text.preamble(r"\usepackage{times}")
        else:
            plot_module = None        

        # Store option controlling difference plotting into global config
        cr_config.do_diff_plotting = options.diff_plotting

        print ''
        for curr_dir in run_test_dirs:
            self.compare_test_dir(options.comparison_dir_name, options.baseline_dir_name, plot_module, curr_dir)

# Main
g_ct_obj = compare_results()
g_ct_obj.main()
