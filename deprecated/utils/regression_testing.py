import os
import re
import sys

bin_dir = os.path.dirname(sys.argv[0])

class rt_config:
    # Name of directory to search for when finding test directories, this dir is copied to the
    # run dir
    std_input_dir_name    = 'std_input'
    std_output_dir_name   = 'std_output'

    # Default location binaries to run with
    binary_dir            = bin_dir + '/bin'

    # Prefix name of binary that is replaced by the run prefix when creating test run dir
    binary_prefix         = "oco_l2"
    run_prefix            = "run"

    # Name of a directory to not be included when doing all tests
    test_devel_dir_name   = 'test_devel/'

    # Symbols to look for using 'nm' to determine if valid binaries
    bin_check_symbol      = 'oco_input_file'
    parallel_check_symbol = 'mpi_comm'

    # Where to start looking for test directories
    default_test_base_dir = os.getcwd()

class regression_testing:

    def find_test_dirs(self, search_base):
        'Return list of standard input directories representing testcases'
        
        findCmd = 'find ' + search_base + ' -name "' + rt_config.std_input_dir_name + '" -maxdepth 4 | sort'
        
        ff = os.popen(findCmd, 'r')
        test_dirs = [x.strip().replace(rt_config.std_input_dir_name, '') for x in ff.readlines()]
        ff.close()
        
        return test_dirs

    def get_binary_names(self, use_parallel):

        bin_names = []
        for bin_dir_file in os.listdir( rt_config.binary_dir ):
            bin_full_path = rt_config.binary_dir + '/' + bin_dir_file
            if self.is_valid_binary(bin_full_path, use_parallel):
                bin_names.append(bin_full_path)

        def bin_name_sort(x, y):
            # Make sure we can split based on a hypen
            if x.find('-') > 0 and x.find('-') < len(x) - 1 and \
               y.find('-') > 0 and y.find('-') < len(y) - 1:
                return cmp(x.split('-')[1], y.split('-')[1])
            else:
                return cmp(x, y)

        bin_names.sort(bin_name_sort)
        bin_names.reverse()
        
        return bin_names

    def is_valid_binary(self, binary_filename, use_parallel):
        'Ensures that a binary is a valid L2 binary'
    
        if not os.path.exists(binary_filename) or os.path.isdir(binary_filename):
            return False

        is_valid_l2 = False
        has_mpi = False

        nm_cmd = 'nm ' + binary_filename
        
        ff = os.popen(nm_cmd, 'r')

        for out_line in ff:
            line_parts = out_line.split()

            if len(line_parts) >= 3:
                symbol = line_parts[2]

                if symbol.lower().find( rt_config.bin_check_symbol.lower() ) >= 0:
                    is_valid_l2 = True
                    
                if symbol.lower().find( rt_config.parallel_check_symbol.lower() ) >= 0:
                    has_mpi = True
                    
        ff.close()

        if is_valid_l2:
            if use_parallel and has_mpi:
                return True
            elif not use_parallel and not has_mpi:
                return True

        return False

    def parse_binary_option(self, binary_name, use_parallel):
        #####
        # Parse options related to the binary name used for testing
        if binary_name == None:
            if not os.path.exists(rt_config.binary_dir):
                print 'Binary path does not exist "%s"\nPlease make a symlink to the L2_EXE binaries directory' % rt_config.binary_dir
                return None

            bin_list = self.get_binary_names(use_parallel)

            if len(bin_list) == 0:
                print 'No binaries available in standard location: "%s"' % rt_config.binary_dir
                return None

            print 'Possible binaries matching arguments:'
            for bin_idx in range(len(bin_list)):
                print '[%d] %s' % (bin_idx, bin_list[bin_idx])
            return None
        elif binary_name.isdigit():
            bin_list = self.get_binary_names(use_parallel)
            bin_int = int(binary_name)
            if bin_int >= len(bin_list):
                print 'Invalid index into list of possible binaries'
                return None
                
            binary_name = bin_list[bin_int]
            
        if not self.is_valid_binary(binary_name, use_parallel):
            print '%s is not a valid binary for running with the supplied arguments' % binary_name
            return None

        # Make sure pass on full path especially when launching on nodes
        binary_name = os.path.realpath( binary_name )

        return binary_name

    def parse_test_names(self, test_base_dir, test_names):
        all_test_dirs = self.find_test_dirs(test_base_dir)
        run_test_dirs = []

        if len(test_names) <= 0:
            print 'Possible test directories:'
            for dir_name in all_test_dirs:
                print dir_name.replace(test_base_dir, '').lstrip('/')
            return None
        elif len(test_names) == 1 and test_names[0].lower() == "all":
            # For all add all but the development testcases
            for dir_name in all_test_dirs:
                if test_base_dir != rt_config.default_test_base_dir or dir_name.find(rt_config.test_devel_dir_name) < 0:
                    run_test_dirs.append(dir_name)
        else:
            for curr_name in test_names:
                for curr_test_dir in all_test_dirs:
                    if re.search(curr_name, curr_test_dir) and ((test_base_dir != rt_config.default_test_base_dir or curr_test_dir.find(rt_config.test_devel_dir_name) < 0) or curr_name.find(rt_config.test_devel_dir_name) >= 0):
                        run_test_dirs.append(curr_test_dir)

        return run_test_dirs
