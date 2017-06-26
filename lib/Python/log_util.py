import sys
import logging
import subprocess

def init_logging(log_level=logging.INFO, name=None, format='%(levelname)s - %(message)s', stream=sys.stdout):
    # Initialize global logging so all messages are reported
    # even if another handler does not want as verbose an output
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Start up console handler
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(log_level)
    console.setFormatter(logging.Formatter(fmt=format))

    logger.addHandler(console)

    return logger

def open_log_file(log_filename, logger=logging.getLogger(), append=False):
    if append:
        mode = 'a'
    else:
        mode = 'w'
        
    fileout = logging.FileHandler(log_filename, mode)
    fileout.setLevel(logging.DEBUG)
   
    logger.addHandler(fileout)

    return fileout
