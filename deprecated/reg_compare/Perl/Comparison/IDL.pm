#  -------------------------------------------------------------------------
#   Name:        IDL Comparison Base Class
#   File:        Comparison/IDL.pm
#   Created:     04/06/2004
#   Comment:
#
#   Contains commmon functionality for running an external IDL comparison
#   routine. A subclass only needs to implement the can_compare_file and
#   the _init routine where it will specify the location of the IDL save 
#   file.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package Comparison::IDL;

use strict;
use Carp;
use File::Basename;
use Cwd;

# Inherit from our base class
use Comparison::Base;
our @ISA = ("Comparison::Base");

# ===============================================================
# Public Class Methods
# ===============================================================

# ---------------------------------------------------------------
# new:
# Creates a new instance of this class
# ---------------------------------------------------------------

sub new {
  my $classname = shift;
  my $self      = $classname->SUPER::new(@_);
  return $self;
}

# ===============================================================
# Public Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# compare_pair:
# The compare method uses the filenames passed as arguments
# and performs a comparison on them. It should output verbose
# information into the REPORT attribute.
# Return:
#    True if the two files have differences false otherwise.
# ---------------------------------------------------------------

sub compare_pair {
  my ($self, $base_filename, $comp_filename) = @_;

  # Set up enviromental variables that pass info to IDL process
  $self->_setup_enviroment_variables($base_filename, $comp_filename);

  # Run idl process to do comparison
  return $self->_run_idl_process();

}

# ===============================================================
# Private Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# _init:
# Initializes attributes of the class
# ---------------------------------------------------------------

sub _init {
  my $self = shift;

  # Call super class's init routine
  $self->SUPER::_init(@_);

  # Initialize attributes
  $self->{IDL_BIN_PATH}  = $ENV{"IDL_DIR"} . "/bin";
  $self->{IDL_COMP_CODE_DIR} = $self->{regression_base_dir} . "/IDL";

  # These attributes need to be initialzed in the sub class
  $self->{IDL_SAVE_FILE} = "";
  $self->{IDL_MODULE_DIR} = "";
  $self->{IDL_PROCEDURE_NAME} = "";

  # Used for inspection mode
  $self->{TMP_STARTUP_FILE} = "tmp_inspection_startup.pro";

  # Initialize test case directory
  my $dir = getcwd(  );
  $self->{TESTCASE_DIRECTORY} = $dir;

  # Make sure there is not some mistake and that the IDL comp code dir
  # exists!
  die "IDL Comparison Code directory does not exist"
    unless -e $self->{IDL_COMP_CODE_DIR};

}

# ---------------------------------------------------------------
# _setup_enviroment_variables:
# ---------------------------------------------------------------

sub _setup_enviroment_variables {
  my ($self, $base_filename, $comp_filename) = @_;

  $ENV{"REG_BASELINE_FILENAME"}   = $base_filename;
  $ENV{"REG_COMPARISON_FILENAME"} = $comp_filename;
  $ENV{"REG_CHECK_ALL"}           = $self->{check_all};

  # This environment variable is needed so that IDL knows where
  # to locate the test case directory.
  $ENV{"TESTCASE_WORKING_DIRECTORY"} = $self->{TESTCASE_DIRECTORY};

}

# ---------------------------------------------------------------
# _run_idl_process:
# Runs IDL process that does comparison
# ---------------------------------------------------------------

sub _run_idl_process {
  my $self = shift;

  unless(-e $self->{IDL_SAVE_FILE}) {
    croak "IDL save file '" . $self->{IDL_SAVE_FILE} . "' does not exist.";
  }

  my $idl_bin = $self->{IDL_BIN_PATH} . "/idl";

  unless(-e $idl_bin) {
    $idl_bin = `which idl`;
    chomp($idl_bin);
  }

  unless(-e $idl_bin) {
    croak "IDL binary '" . $idl_bin . "' does not exist.";
  }

  my $idl_exit_code = 0;

  unless( $self->{debug} ) {
    my $idl_cmd = "$idl_bin -32 -quiet -rt='" .
      $self->{IDL_SAVE_FILE} . "'" ;

    # Run IDL process redirecting STDERR to STDOUT during the run
    open(IDL_PROCESS, "$idl_cmd 2>&1 |");

    my $execution_halted = 0;
    # Process remaining lines
    while (<IDL_PROCESS>) {

      # Gets rid of the IDL lines that need not be printed
      if($_ =~ /^\%\sRestored file/) {
	next;
      }
      $self->_process_idl_output_line($_);
      $execution_halted = 1 if /Execution halted at/;
    }

    close(IDL_PROCESS);

    if ( $execution_halted ) {
      # If the IDL program's execution was halted unexpectedly there
      # was most likely an error and so there must be problems in the
      # comparison
      $idl_exit_code = 1;
    } else {
      # Find out whether differences were found
      $idl_exit_code = ($? >> 8);
    }
  }

  # If we run in inspection mode then IDL is run interactively for failures
  if ( $self->{debug} || ($self->{inspect} && $idl_exit_code) ) {
    # Instructs IDL comparison routine to stop if there are differences
    $ENV{"REG_INSPECT"} = $self->{inspect};

    $self->_report("Differences found. Now running comparison module in inspection mode\n");

    $self->_create_inspect_startup_file();
    my $idl_cmd = $idl_bin . " " . $self->{TMP_STARTUP_FILE};

    system($idl_cmd);

    # Remove the temporary startup file that we created before running IDL
    unlink  $self->{TMP_STARTUP_FILE};
  }

  return $idl_exit_code;
}

# ---------------------------------------------------------------
# _process_idl_output_line:
# Processes an output line from an IDL process in any way the 
# class wants to
# ---------------------------------------------------------------

sub _process_idl_output_line {
  my $self = shift;
  my $output_line = shift;

  if($output_line =~ /^\-/ || $output_line =~ /^Comparison/) {
    $self->_report($output_line, 1);
  } else {
    $self->_report($output_line);
  }
}

# ---------------------------------------------------------------
# _create_inspect_startup_file
# Creates a temporary startup file which will be used to load IDL
# in an interactive like mode.
# ---------------------------------------------------------------

sub _create_inspect_startup_file {
  my $self = shift;

  die "No module code directory specified" unless $self->{IDL_MODULE_DIR};

  unless ( -e $self->{IDL_MODULE_DIR} ) {
    die "Module code dir " .$self->{IDL_MODULE_DIR} . " does not exist";
  }

  die "No procedure name specified" unless $self->{IDL_PROCEDURE_NAME};

  my $startup_contents = <<STARTUP_CONTENTS;

L1BRootDir = '/vobs/L1B/idl/'
L1BSharedDir   = '/vobs/Shared/L1B/idl'
RegCommon   = '$self->{IDL_COMP_CODE_DIR}/Common'
CompModule  = '$self->{IDL_MODULE_DIR}'

dirsToExpand = [                \$
  L1BRootDir,                   \$
  L1BSharedDir,                 \$
  RegCommon,                    \$
  CompModule                    \$
]

; Expand and added to !PATH all directories with .pro files.
;
CASE !D.NAME OF         & \$
  'VMS' : delim = ','   & \$
  'WIN' : delim = ';'   & \$
ELSE    : delim = ':'   & \$
ENDCASE

FOR iPath=0, N_Elements(dirsToExpand)-1L DO \$
  !PATH = Expand_Path('+' + dirsToExpand[iPath]) + delim + !PATH

$self->{IDL_PROCEDURE_NAME}

STARTUP_CONTENTS

  open(STARTUP_FILE, ">" . $self->{TMP_STARTUP_FILE});
  print STARTUP_FILE $startup_contents;
  close(STARTUP_FILE);
}

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
