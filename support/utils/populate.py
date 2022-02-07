#!/usr/bin/env python

import os
import re
import sys
import logging
import traceback

from optparse import OptionParser
from full_physics import PopulatorBase, parse_keyval_str_list, docopt_simple
import full_physics.log_util as log_util

version = "September 18, 2019"
usage = '''Usage:
  populate.py [options] <config_file> [<config_file>...]
  populate.py -h | --help
  populate.py -v | --version

Create the run scripts needed to run Level 2 Full Physics retrieval code
on the given set of configuration files.

Options:
  -h, --help
      Print this message

  -b, --binary=FILE
      Location of L2 FP binary to use for jobs

  -c, --l2_config=FILE
      Name of config file to use

  -a, --aggregate       
      After generating level 2 output files, aggregate them into one
      common file.

  --use-subdirectory
      Break log and output files into subdirectories by the grouping given
      with -g. We also aggregate the subdirectories separately, which can
      significantly reduce the time it takes to do aggregation.

  --skip-check
      Skip check to see if data is processable. This can take a while, and
      is unnecessary if we already know the data is processable (e.g., we
      generated the QTS data by looking at real data file).

  --email-address=d
      Email address to mail when aggregation starts and ends.

  -g, --group_size=NUM
      Number of L2 jobs to group together in a single instance
      on the cluster.
      [default: 1]

  -p, --parallel_size=NUM
      Number of L2 jobs from a group to run in parallel
      [default: 1]
    
  --absco-version=V
      Version of ABSCO to use. This is only used by Pleiades runs.
      [default: v5.2_final]

  -t, --target_cluster=NAME
      Name of cluster management system to be run on.
      choices: [pbs_pro, torque, pleiades]
      [default: pbs_pro]

  -q, --quiet           
      Do not output details of program operation to console

  -l, --log_file=FILE
      Log file where populator operations are stored instead of printed 
      to stdout
'''

args = docopt_simple(usage, version=version)

# Logger for file operations
logger = logging.getLogger(os.path.basename(__file__))

if args.quiet:
    log_util.init_logging(logging.ERROR)
else:
    log_util.init_logging(logging.INFO)
        
# Initialize logging
if args.log_file:
    log_obj = log_util.open_log_file(args.log_file)
    log_obj.setFormatter( logging.Formatter("%(asctime)s: %(name)25s - %(levelname)7s - %(message)s") )
else:
    log_obj = None
populate_options = {}
if(args.binary):
    populate_options["l2_binary_filename"] = args.binary
if(args.l2_config):
    populate_options["l2_config_filename"] = args.l2_config
populate_options["skip_check"] = args.skip_check
populate_options["email_address"] = ""
if(args.email_address is not None):
    populate_options["email_address"] = args.email_address
populate_options["aggregate"] = args.aggregate
populate_options["use_subdirectory"] = args.use_subdirectory
populate_options["abscoversion"] = args.absco_version
populate_options["target_cluster"] = args.target_cluster
populate_options["group_size"] = int(args.group_size)
populate_options["parallel_size"] = int(args.parallel_size)

for config_file in args.config_file:
    logger.info('Populating using configuration file: %s' % config_file)

    # Log any exceptions to disk then re-raise
    success = False
    try:
        p = PopulatorBase.create_populator_from_config_file(config_file, **populate_options)
        if(p is None):
            raise LookupError('Could not determine config type from config contents or filename')
        success = p.populate(config_file)
    except:
        # Will cause exception to print to screen
        logger.error('Error when populating using file: %s' % config_file)
        for tb_line in traceback.format_exception(*sys.exc_info()):
            logger.error(tb_line.strip())

if success:
    logger.info("Population was successful")
    sys.exit(0)
else:
    logger.error("Population failed")
    sys.exit(1)
