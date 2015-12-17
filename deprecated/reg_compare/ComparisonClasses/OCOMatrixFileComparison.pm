#  -------------------------------------------------------------------------
#   Name:        OCO Matrix File Comparison
#   File:        OCOMatrixFileComparison.pm
#   Created:     03/27/2006
#   Comment:
#
#   Compares files that were written out by OCO L2 PGE
#
#  -------------------------------------------------------------------------

package OCOMatrixFileComparison;

use strict;
use Carp;
use Cwd;

use File::Basename;

# Inherit from our base class
use Comparison::IDL;
our @ISA = ("Comparison::IDL");

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

  unless ( defined($self->{comparison_map_file}) ) {
    # File to pass comparison map file contents, use pid to prevent overlap from
    # simultatenous runs in the same dir
    $self->{comparison_map_file} = "reg_comp_match_map_" . getppid() . ".txt";
  }

  unless ( defined($self->{delete_comparison_map}) ) {
    $self->{delete_comparison_map} = 1;
  }

  return $self;
}

# ---------------------------------------------------------------
# can_compare_files:
# Returns true if the class determines it can do a comparison of
# the type of file specified by the paramater filename.
# Parameter:
#    The baseline and comparison filenames that would be compared
# Return:
#    True if a comparison can be performed false otherwise
# ---------------------------------------------------------------

sub can_compare_files {
  my $classname = shift;
  my($baseline_filename, $comparison_filename) = @_;

  if ( $baseline_filename =~ /.dat/ ) {
    my $grep_cmd = "grep '^[ \t]*begin' $baseline_filename";
    open(GREP, "$grep_cmd |");
    my $has_header = <GREP>;
    chomp($has_header);
    close(GREP);

    if( $has_header ne "" ) {
      return 1;
    } else {
      return undef;
    }
  } else {
    return undef;
  }
}

# ===============================================================
# Public Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# compare:
# Method that is called to perform comparison actions on file
# pairs listed in class.
#
# An filename is created with the list of file pairs to compare.
# This filename is passed to an IDL module which loops over the
# files to perform the comparisons.
#
# Return:
#    Array: (number diff files, number files compared)
# ---------------------------------------------------------------

sub compare {
  my $self = shift;

  $self->_setup_enviroment_variables();
  $self->_create_pair_list_file();
  $self->_run_idl_process();
  my($num_diff, $num_comp) = $self->_retrieve_counts();
  $self->_remove_pair_list_file();

  return ($num_diff, $num_comp);
}

# ---------------------------------------------------------------
# compare_pair:
# For fulfillment of the interface, this method creates a file
# pair file for the IDL module that contains only one pair.
#
# Return:
#    True if the two files have differences false otherwise.
# ---------------------------------------------------------------

sub compare_pair {
  my ($self, $base_filename, $comp_filename) = @_;

  my $single_pair = {};
  $single_pair->{$base_filename} = $comp_filename;

  $self->_setup_enviroment_variables();
  $self->_create_pair_list_file($single_pair);
  $self->_run_idl_process();
  my($num_diff, $num_comp) = $self->_retrieve_counts();
  $self->_remove_pair_list_file();

  return $num_diff > 0 ? 1 : undef;
}

# ===============================================================
# Private Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# init:
# Initializes attributes of the class
# ---------------------------------------------------------------

sub _init {
  my $self = shift;

  # Call super class's init routine
  $self->SUPER::_init(@_);

  # Initialize attributes
  $self->{IDL_SAVE_FILE} =
    $self->{IDL_COMP_CODE_DIR} . "/" . "oco_l2_matrix_compare.sav";

  $self->{IDL_MODULE_DIR} = 
    $self->{IDL_COMP_CODE_DIR} . "/" . "OCOL2MatrixCompare";

  $self->{IDL_PROCEDURE_NAME} = "oco_l2_matrix_compare";

  # Initialize testcase directory
  my $dir = getcwd();
  $self->{TESTCASE_DIRECTORY} = $dir;

}

# ---------------------------------------------------------------
# _setup_enviroment_variables:
# Sets up enviromental variables necessary for performing
# comparison. Overrides definition of that in IDL.pm since this
# module passes a filename with a list of pairs.
# ---------------------------------------------------------------

sub _setup_enviroment_variables {
  my ($self) = shift;

  # This environmental variable is needed so that IDL knows where
  # to locate the test case directory.
  $ENV{"TESTCASE_WORKING_DIRECTORY"} = $self->{TESTCASE_DIRECTORY};

  $ENV{"REG_PAIR_LIST_FILENAME"}  = $self->{comparison_map_file};
  $ENV{"REG_CHECK_ALL"}           = $self->{check_all};

}

# ---------------------------------------------------------------
# _create_pair_list_file:
# Creates a pair list file which will be used by the IDL module
# to know which files to loop over. If a pair map is supplied it 
# is used instead of the class attribute.
#
# Return:
# Filename of pair list file
# ---------------------------------------------------------------

sub _create_pair_list_file {
  my($self, $pair_map) = @_;

  # Use supplied map instead of object attribute
  unless ( defined($pair_map) ) {
    $pair_map = $self->{compare_pair_map};
  }

  open(PAIR_FILE, ">" . $self->{comparison_map_file}) or
    die "Count not write to pair map file: " . $self->{comparison_map_file};

  foreach my $base_filename (keys(%{$pair_map})) {
    my $comp_filename = $pair_map->{$base_filename};

    if(-f $base_filename && -f $comp_filename) {
	print PAIR_FILE $base_filename . "\t" . $comp_filename . "\n";
      }
  }

  close(PAIR_FILE);

}

# ---------------------------------------------------------------
# _retrieve_counts:
# Retrieves the counts for how many files are different and
# were compared from the pairs file which is overwritten with the
# values by the IDL module.
# ---------------------------------------------------------------

sub _retrieve_counts {
  my ($self) = shift;

  open(PAIR_FILE, $self->{comparison_map_file}) or
    die "Count not read from pair map file: " . $self->{comparison_map_file};

  my @lines = <PAIR_FILE>;
  close(PAIR_FILE);

  my $num_files_compared = scalar(@lines);

  my $num_diff_files = 0;
  foreach my $pair_count_line (@lines) {
    chomp($pair_count_line);
    my($base_name, $comp_name, $num_pair_diffs) = split(/\s+/, $pair_count_line);
    $num_diff_files++ if $num_pair_diffs > 0;
  }

  return ($num_diff_files, $num_files_compared);
}

# ---------------------------------------------------------------
# _remove_pair_list_file:
# Deletes the pair map file passed to IDL
# ---------------------------------------------------------------

sub _remove_pair_list_file {
  my($self) = @_;

  if ( $self->{delete_comparison_map} ) {
    unlink($self->{comparison_map_file});
  }
}

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
