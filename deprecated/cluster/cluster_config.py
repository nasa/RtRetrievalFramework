import os
import sys
import copy

COSMOS_COMMON = '-o "%J.out" -e "%J.err"'

SCF_USERNAME = os.environ.get('SCF_USERNAME', None) and os.environ.get('SCF_USERNAME') or os.environ.get('USER')

# Interface through buzz into fullerene
RUNS_REMOTE_MACHINE = { 'nebula':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'galaxy':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'cosmos':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'fullerene':     '',
                        'ococluster01':  'ococluster01.cira.colostate.edu',
                        }

WORK_REMOTE_MACHINE = { 'nebula':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'galaxy':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'cosmos':        '%s@buzz.jpl.nasa.gov' % SCF_USERNAME,
                        'fullerene':     '',
                        'ococluster01':  '%s@nephthys.jpl.nasa.gov' % SCF_USERNAME,
                        }

SCRATCH_BASE = { 'galaxy': '/lscratch/' + os.environ['USER'],
                 'nebula': '/lscratch/' + os.environ['USER'],
                 'ococluster01': '/data4/groups/algorithm/%s/scratch' % os.environ['USER'],
                 'cosmos':    { 'single':   '/lscratch/' + os.environ['USER'],
                                'parallel': '/work03/bfk/%s/scratch' % os.environ['USER'] },
                 'fullerene': { 'single':   '/state/partition1/' + os.environ['USER'],
                                'parallel': None,
                                }               
                 }
                 
SUBMIT_EXEC = { 'nebula':
                { 'long':  'qsub -V -q longq',
                  'short': 'qsub -V -q shortq',
                  },
                'galaxy': 'qsub -V',
                'fullerene':
                { 'long':  'qsub -V -q long',
                  'short': 'qsub -V -q short',
                  'ops':   'qsub -V -q ops',
                  'amd':   'qsub -V -q amd',
                  },
                'cosmos':
                { 'long':  'bsub -W 720 -q long ' + COSMOS_COMMON, 
                  'short': 'bsub -W 180 -q short ' + COSMOS_COMMON,
                  'debug': 'bsub -W 30 -q debug ' + COSMOS_COMMON,
                  },
                'ococluster01': 'qsub -V',
                }

JOB_ARRAY_OPT = { 'nebula':    '-J %d-%d',
                  'galaxy':    '-J %d-%d',
                  'fullerene': '-t %d-%d',
                  'cosmos':    '-n 1 -J "l2_fp[%d-%d]"',
                  'ococluster01': '-t %d-%d',
                  }

SINGLE_MODE_OPT = { 'nebula':    '',
                    'galaxy':    '',
                    'fullerene': '',
                    'cosmos':    '-n 1',
                    'ococluster01': '',
                    }

PARALLEL_SUBMIT_OPT    = { 'fullerene': '-lnodes=%d:ppn=%d',
                           'cosmos':    '-n %d -a mpich_gm',
                           'ococluster01':'-lnodes=%d:ppn=%d',
                           }
PARALLEL_RUN_PREFIX    = { 'fullerene': 'mpiexec -n %d',
                           'cosmos':    'mpirun.lsf',
                           'ococluster01':'mpirun_now -np %d',
                           }
PARALLEL_PROC_PER_NODE = { 'fullerene': 4,
                           'cosmos': 1,
                           'ococluster01':4,
                           }
PARALLEL_LAUNCH_MPD = { 'fullerene': True, 'cosmos': False, 'ococluster01': False}

ABSCO_DIRS = {'nebula':     '/work00/bfk/mcduffie/absorber/v3.0.0_alpha/',
              'galaxy':     '/work01/bfk/mcduffie/absorber/v3.0.0_alpha/',
              'fullerene': '/groups/algorithm/l2_fp/absco/v3.0.0_alpha/',
              'cosmos':    '/work03/bfk/mcduffie/absorber/v3.0.0_alpha/',
              'ococluster01': '/data4/OCO/l2/absco/v3.0.0_alpha/',
              }

QUERY_SCRIPT_CONFIG = { 'nebula':    { 'query_command': 'qstat', 'index_tmpl': '%s[%d]', },
                        'galaxy':    { 'query_command': 'qstat', 'index_tmpl': '%s[%d]', },
                        'fullerene': { 'query_command': 'qstat', 'index_tmpl': '%s-%d',  },
                        'cosmos':    { 'query_command': 'bstat', 'index_tmpl': '%s[%d]', },
                        'ococluster01': { 'query_command': 'qstat', 'index_tmpl': '%s[%d]', },
                      }

RUN_SCRIPT = os.path.dirname(sys.argv[0]) + '/run_job.py'
JOB_ARRAY_SCRIPT_TEMPLATE = os.path.dirname(sys.argv[0]) + '/job_array_template.py'
JOB_ARRAY_SCRIPT_OUTPUT = 'l2_run_%s.py'

QUERY_SCRIPT_TEMPLATE = os.path.dirname(sys.argv[0]) + '/query_script_template.py'
QUERY_SCRIPT_OUTPUT = 'query_job_%s.py'

GLOBAL_RUN_TMP_DIR = '%s/%s' % (os.environ['HOME'].rstrip('/'), '.l2_fp_run_tmp')
RSYNC_TIME_QUEUE_DIR = { 'nebula': GLOBAL_RUN_TMP_DIR,
                         'galaxy': GLOBAL_RUN_TMP_DIR,
                         'cosmos': GLOBAL_RUN_TMP_DIR,
                         'fullerene': GLOBAL_RUN_TMP_DIR,
                         'ococluster01': GLOBAL_RUN_TMP_DIR,
                         }

# Rush ssh quietly and do not do any x11 forwarding which slows down the connection
SSH_EXEC = 'ssh -q -o ForwardX11=no'

TIME_OUTPUT_FILE = 'runtime.txt'
TIME_COMMAND = '/usr/bin/time -v -o "%s"' % TIME_OUTPUT_FILE

######

VERBOSE_RUN_JOB = True

RSYNC_BIN = "rsync"
RSYNC_GENERAL_OPTS = '-Cauzb --delete-after --partial --timeout=120'
RSYNC_INPUT_DFLT_OPTS = RSYNC_GENERAL_OPTS + ' --copy-unsafe-links'
RSYNC_OUTPUT_DFLT_OPTS = RSYNC_GENERAL_OPTS

RSYNC_INPUT_OPTS = { 'nebula': RSYNC_INPUT_DFLT_OPTS,
                     'galaxy': RSYNC_INPUT_DFLT_OPTS,
                     'cosmos': RSYNC_INPUT_DFLT_OPTS,
                     'fullerene': RSYNC_INPUT_DFLT_OPTS,
                     'ococluster01': RSYNC_INPUT_DFLT_OPTS + ' --no-g',
                     }
RSYNC_OUTPUT_OPTS = { 'nebula': RSYNC_OUTPUT_DFLT_OPTS,
                      'galaxy': RSYNC_OUTPUT_DFLT_OPTS,
                      'cosmos': RSYNC_OUTPUT_DFLT_OPTS,
                      'fullerene': RSYNC_OUTPUT_DFLT_OPTS,
                      'ococluster01': RSYNC_OUTPUT_DFLT_OPTS + ' --no-g',
                      }

RSYNC_NUM_SIMUL = 40

RSYNC_REMOTE_OPTS = '--rsh "%s"' % SSH_EXEC

RSYNC_RETRY_EXIT_CODES = [12, 23, 30, 255, 999]
IN_RETRY_COUNT = 20
OUT_RETRY_COUNT = 40
RSYNC_RETRY_WAIT = 60

QUEUE_RETRY_COUNT = 180*2
QUEUE_RETRY_WAIT = 30

STDOUT_FILENAME = "stdout"
STDERR_FILENAME = "stderr"
RUN_FILENAME = "oco_l2.run"
RUN_MODE_KEYWORD = "run_mode"

SHARED_SYNC_LIST = [STDOUT_FILENAME, STDERR_FILENAME, TIME_OUTPUT_FILE, "*.log", "out/", "fort.*", "core*"]

SYNC_OUT_LIST = { 'FORWARD_MODEL' : SHARED_SYNC_LIST + ["in/l1b/spec/"],
                  'JACOBIAN_ONLY' : SHARED_SYNC_LIST + ["in/l1b/spec/"],
                  'RETRIEVAL'     : SHARED_SYNC_LIST,
                  'SOLAR'         : SHARED_SYNC_LIST,
                  }

CLEAN_FILE_LIST = SYNC_OUT_LIST

ERROR_CHECK_STRINGS = [ "Error in forward model",
                        "RADIANT ACTION:",
                        "FATAL:",
                        "glibc detected",
                        "Program terminated",
                        "Arithmetic exception:",
                        "FORTRAN Runtime Error",
                        ]
LOG_VERBOSE_PREFIX = '--->'
