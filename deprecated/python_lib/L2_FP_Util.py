import os
import sys

try:
    import platform

    py_minor_ver = int(platform.python_version_tuple()[1])

    if py_minor_ver <= 3:
        import popen2
    else:
        import subprocess
except:
    py_minor_ver = 2
    import popen2

# Symbols to look for using 'nm' to determine if valid binaries
BIN_CHECK_SYMBOL      = 'oco_input_file'
PARALLEL_CHECK_SYMBOL = 'mpi_comm'

# In case there is an error just running an external command
BOGUS_ERROR_CODE = 999

def is_valid_binary(binary_filename, use_parallel):
    'Ensures that a binary is a valid L2 binary'

    if not os.path.exists(binary_filename) or os.path.isdir(binary_filename):
        return False

    is_valid_l2 = False
    has_mpi     = False

    nm_cmd = 'nm ' + binary_filename

    ff = os.popen(nm_cmd, 'r')

    for out_line in ff:
        line_parts = out_line.split()

        if len(line_parts) >= 3:
            symbol = line_parts[2]

            if symbol.lower().find( BIN_CHECK_SYMBOL.lower() ) >= 0:
                is_valid_l2 = True

            if symbol.lower().find( PARALLEL_CHECK_SYMBOL.lower() ) >= 0:
                has_mpi = True

    ff.close()

    if is_valid_l2:
        if use_parallel and has_mpi:
            return True
        elif not use_parallel and not has_mpi:
            return True

    return False

def run_wait_command(command, verbose=True, show_output=True):

    if verbose:
        print 'Running and waiting on command: %s' % command

    if py_minor_ver <= 3:
        run_instance = popen2.Popen3(command, True)
    else:
        if show_output:
            run_instance = subprocess.Popen(command, stdin=subprocess.PIPE, shell=True, close_fds=True)
        else:
            run_instance = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, close_fds=True)

    if py_minor_ver <= 3:
        if verbose or show_output:
            sys.stdout.writelines(run_instance.fromchild.readlines())
            sys.stderr.writelines(run_instance.childerr.readlines())

    status_code = run_instance.wait()

    if os.WIFEXITED(status_code):
        error_code = os.WEXITSTATUS(status_code)
    else:
        error_code = BOGUS_ERROR_CODE

    return error_code
