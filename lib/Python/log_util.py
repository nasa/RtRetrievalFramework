import logging
import subprocess

# For traceability tracking
SVN_COMMAND = 'svn info -R'

def init_logging(log_level=logging.INFO, format='%(levelname)s - %(message)s'):
    # Initialize global logging so all messages are reported
    # even if another handler does not want as verbose an output
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Start up console handler
    console = logging.StreamHandler()
    console.setLevel(log_level)
    console.setFormatter(logging.Formatter(fmt=format))

    logger.addHandler(console)

def log_revision_info(check_path='./'):

    logger = logging.getLogger('subversion')

    logger.debug('Loading subversion revision information')
    logger.debug('Running subversion command: %s' % SVN_COMMAND)
    info_process = subprocess.Popen('%s %s' % (SVN_COMMAND, check_path), shell=True, close_fds=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

    for line in info_process.stdout.readlines():
        logger.debug(line.strip())
    info_process.wait()

    # If the subinfo_process is interrupted then exit the python program
    if info_process.returncode < 0:
        raise OSError('%s subprocess exited abnormally' % SVN_COMMAND)

def open_log_file(log_filename, append=False):
    logger = logging.getLogger()

    if append:
        mode = 'a'
    else:
        mode = 'w'
        
    fileout = logging.FileHandler(log_filename, mode)
    fileout.setLevel(logging.DEBUG)
   
    logger.addHandler(fileout)

    return fileout

def close_log_file(handler):
    logger = logging.getLogger()
    logger.removeHandler(handler)
    handler.close()
