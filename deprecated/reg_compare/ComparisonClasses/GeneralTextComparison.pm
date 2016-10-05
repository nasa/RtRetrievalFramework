#  -------------------------------------------------------------------------
#   Name:        General Text Comparison Module
#   File:        GeneralTextComparison.pm
#   Created:     04/05/04
#   Comment:
#
#   Does a simple textual comparison using diff.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package GeneralTextComparison;

use strict;
use Carp;

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

  # Files both appear to be text files
  return (-T $baseline_filename && -T $comparison_filename);
}

# ---------------------------------------------------------------
# get_priority:
# Return the lowest priority. We should only be run as a last
# resort.
# ---------------------------------------------------------------

sub get_priority {
  my $classname = shift;

  return 0;
}

# ===============================================================
# Public Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# compare_pair:
# The compare method uses the filenames given as arguments
# and performs a comparison on them.
# Return:
#    True if the two files have differences false otherwise.
# ---------------------------------------------------------------

sub compare_pair {
  my ($self, $base_filename, $comp_filename) = @_;

  my $diff_command =
    "diff " . $base_filename . " " . $comp_filename;

  open(DIFF, "$diff_command |");
  my @diff_results = <DIFF>;
  close(DIFF);

  # Parse counts of lines from output
  my $base_diff_count = 0;
  my $comp_diff_count = 0;
  foreach my $count (@diff_results) {

    if($count =~ /^(\d+),?(\d+)?\w(\d+),?(\d+)?$/) {

      if(defined($2)) {
	$base_diff_count += $2 - $1 + 1;
      } else {
	$base_diff_count++;
      }

      if(defined($4)) {
	$comp_diff_count += $4 - $3 + 1;
      } else {
	$comp_diff_count++;
      }

    } # end if matches count line

    # Add each line to the report
    $self->_report($count);

  } # end for each diff results

  # Take the greater num of differences
  my $num_differences = $base_diff_count > $comp_diff_count ?
    $base_diff_count : $comp_diff_count;

  if($num_differences > 0) {
    if($self->{check_all}) {
      # Diff always checks the whole file so we always report the
      # number of diffs
      $self->_report("Files have $num_differences differences\n");
    }

    return 1;
  } else {
    return undef;
  }
}

# ===============================================================
# Private Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
