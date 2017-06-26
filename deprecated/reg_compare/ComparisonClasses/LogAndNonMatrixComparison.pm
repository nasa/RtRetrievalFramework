#  -------------------------------------------------------------------------
#   Name:        Log File and Non Matrix Text Comparison Module
#   File:        LogAndNonMatrixComparison.pm
#   Created:     10/11/08
#   Comment:
#
#   Compares log files and non matrix .dat files but does not count 
#   date stamps, svn revisons, etc.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package LogAndNonMatrixComparison;

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

  my @tolerable_pats =
    qw( OCO\sL2\sVERSION\s.*
        Date\s+\d{14}\.\d{3}
	Built\s+.*\s+on\s+
	SVN\srevision\s.*
        INFO\:\s+Total.*time.*((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
	INFO\:\s+GEN_SPECTRA\sRUNTIME\:\s+((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
	INFO\:\s+RUNTIME\sCOMPUTE_SPEC_PD\s+((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
	INFO\:\s+RUNTIME\sRETRIEVAL\s+((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
	INFO\:\s+RUNTIME\sITERATION\s+\d+\s+((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
	INFO\:\s+Time\sfor.*per\swn\s+((\d*\.\d+)([Ee][\+-]?\d+)?|([1-9]\d*[Ee][\+-]?\d+)|(\d+\.))
      );

  $self->{tolerable_patterns} = \@tolerable_pats;

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

  if ( $baseline_filename =~ /.dat/ ||
       $baseline_filename =~ /.log/ ||
       $baseline_filename =~ /stdout/ ||
       $baseline_filename =~ /stderr/ ) {

    return 1;
  } else {
    return undef;
  }
}

# ---------------------------------------------------------------
# get_priority:
# Return the lowest priority. We should only be run as a last
# resort.
# ---------------------------------------------------------------

sub get_priority {
  my $classname = shift;

  return 5;
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
  my $diff_results = join("", <DIFF>);
  close(DIFF);

  $diff_results =~ s/^[\s\n]+//;
  my @diff_sections = split(/^(\d+,?(?:\d+)?\w\d+,?(?:\d+)?)$/m, $diff_results);

  # Take out empty section at beginning
  shift(@diff_sections);

  # Parse counts of lines from output
  my $base_diff_count = 0;
  my $comp_diff_count = 0;

  for (my $sect_index = 0; $sect_index < scalar(@diff_sections); $sect_index += 2) {

    my $diff_count    = $diff_sections[$sect_index];
    my $diff_contents = $diff_sections[$sect_index+1];

    my $diff_check = $diff_contents;

    # Check each of the tolerable diff patterns for removal
    # from the difference section
    foreach my $tolerable_pat (@{$self->{tolerable_patterns}}) {
      $diff_check =~ s/$tolerable_pat//g;
    }

    # Extract base and comparison parts of the modified diff section
    my ($base_part, $comp_part) = split(/^---$/m, $diff_check, 2);

    # Remove direction indicators
    $base_part =~ s/^\<//mg;
    $comp_part =~ s/^\>//mg;

    # Check once more that the difference is not tolerable
    if ($base_part ne $comp_part) {

      $diff_count =~ /^(\d+),?(\d+)?\w(\d+),?(\d+)?.*$/;

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

      # Add each line to the report
#      $self->_report($diff_check);
      $self->_report($diff_count . $diff_contents . "\n");
    }

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
