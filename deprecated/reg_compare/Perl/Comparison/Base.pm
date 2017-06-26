#  -------------------------------------------------------------------------
#   Name:        Comparison Base Defnition Class
#   File:        Comparison/Base.pm
#   Created:     04/02/2004
#   Comment:
#
#   This class serves as mostly an abstract definition of the contract that
#   comparison classes need to override.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package Comparison::Base;

use strict;

use POSIX;
use Text::Wrap;
use Carp;

# Wraps the reported text at a certain column
my $WRAP_COLUMN = 80;

# ===============================================================
# Public Class Methods
# ===============================================================

# ---------------------------------------------------------------
# new:
# Creates a new instance of this class
# ---------------------------------------------------------------

sub new {
  my $classname = shift;
  my $self      = { };
  bless($self, $classname);
  $self->_init(@_);
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

  croak "Comparison::Base->can_compare_file has either not been overridden " .
    "or is being called directly."
}

# ---------------------------------------------------------------
# get_priority:
# Returns the priority for consideration by the Finder class.
# A higher number indicates a higher priority. Numbers should be
# greater than or equal to zero.
# ---------------------------------------------------------------

sub get_priority {
  my $classname = shift;

  # Override this method to produce comparors with lower or
  # higher priorities. This serves as a sort of default
  return 100;
}

# ===============================================================
# Private Class Methods
# ===============================================================

# ---------------------------------------------------------------
# _compare_scalar_values:
# Compares two scalar values and returns true if the values match
# and false otherwise
# ---------------------------------------------------------------

sub _compare_scalar_values {
  my ($baseline_val, $comparison_val) = @_;

  my $values_match = undef;
  my $frac_diff = undef;
  if($baseline_val =~ /^-?[0-9]+\.[0-9+e]+/) {
    my $eps = POSIX::DBL_EPSILON;
    my $xmin = POSIX::DBL_MIN;

    # Start off assuming that value match and have tests prove that wrong
    my $mBaseVal = $baseline_val + 2*$xmin;
    my $mCompVal = $comparison_val + 2*$xmin;

    my $absDiff = abs($mCompVal - $mBaseVal);
    $frac_diff = $absDiff / abs($mBaseVal);

    if (lc($frac_diff) eq "nan") {
      # Really big or really small numbers, compare as strings
      $values_match = $baseline_val eq $comparison_val;
    } else {
      $values_match = ($frac_diff < $eps);
    }

  } else {
    $values_match = $baseline_val eq $comparison_val;
  } # end if baseline value digit or not

  # Return fractional difference when called in list mode
  if (wantarray()) {
    return ($values_match, $frac_diff);
  } elsif (defined wantarray()) {
    return $values_match;
  } else {
    return;
  }
}

# ===============================================================
# Public Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# compare:
# Method that is called to perform comparison actions on file
# pairs listed in class. By default compare_all will be called.
#
# Return:
#    Array: (number diff files, number files compared)
# ---------------------------------------------------------------

sub compare {
  my $self = shift;

  return $self->compare_all();
}

# ---------------------------------------------------------------
# compare_all:
# Iterates through the hash of matched filename pairs using the
# compare_one method for each pair.
# Return:
#    Array: (number diff files, number files compared)
# ---------------------------------------------------------------

sub compare_all {
  my $self = shift;

  my $number_diff_files = 0;
  my $number_compared_files = 0;
  foreach my $base_filename ( keys( %{ $self->{compare_pair_map} } ) ) {
    my $comp_filename = $self->{compare_pair_map}->{$base_filename};

    # Make sure we are looking at two files
    if(-f $base_filename && -f $comp_filename) {

      # Increment number of files we have compared for output reporting
      $number_compared_files++;

      $self->_report_buffer_enable();

      $self->_report( "-" x 80 . "\n", 1 );
      $self->_report( "Comparison of $base_filename with $comp_filename\n\n", 1 );

      # Run comparison object on files
      my $has_differences = $self->compare_pair($base_filename, $comp_filename);

      if($has_differences) {
	$number_diff_files++;

	$self->_report_buffer_disable();
	$self->_report_buffer_output();
	
	unless($self->{check_all}) {
	  $self->_report( "-" x 80 . "\n", 1 );
	  $self->_report( "Regression test case has differences\n" );

	  return($number_diff_files, $number_compared_files);
	}
      } else {
	$self->_report( "Files are identical within tolerances\n" );

	$self->_report_buffer_disable();
	$self->_report_buffer_clear();
      }

    } elsif (-f $base_filename && -d $comp_filename) {
      $self->_report( "WARNING: Baseline filename: $base_filename is a file, " .
		      "however comparison filename: $comp_filename is a " .
		      "directory.\n" );
    } elsif (-d $base_filename && -f $comp_filename) {
      $self->_report( "WARNING: Baseline filename: $base_filename is a " .
		      "directory, however comparison filename: " .
		      "$comp_filename is a file.\n" );
    } elsif (-d $base_filename && -d $comp_filename) {
      # Ignore this case since two directories having the same name is ok
    } else {
      $self->_report( "WARNING: baseline filename: $base_filename and " .
		      "comparison filename: $comp_filename are neither " .
		      "directories or files\n!" );
    } # End type of file comparison

  } # end of foreach file pair

  return($number_diff_files, $number_compared_files);

}


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

  croak "Comparison::Base->compare_pair has either not been overridden or is " .
    " being called directly."
}

# ---------------------------------------------------------------
# add_file_pair:
# Adds a pair filenames to the compare file map that need to be
# compared.
# ---------------------------------------------------------------

sub add_file_pair {
  my($self, $base_filename, $comp_filename) = @_;

  $self->{compare_pair_map}->{$base_filename} = $comp_filename;

}

# ---------------------------------------------------------------
# get_file_pairs:
# Returns the map of file pairs
# ---------------------------------------------------------------

sub get_file_pairs {
  my($self) = @_;

  return $self->{compare_pair_map};
}

# ---------------------------------------------------------------
# num_file_pairs:
# Returns the number of file pairs
# ---------------------------------------------------------------

sub num_file_pairs {
  my($self) = @_;

  return scalar(keys(%{$self->{compare_pair_map}}));
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

  # Default values for parameters
  $self->{check_all}           = undef;
  $self->{regression_base_dir} = "";
  $self->{report_function}     = \&print;

  # Buffer the reported strings instead of outputting directly, allowing for
  # strings to be disregarded if there are no differences
  $self->{report_buffer_enable} = 0;
  $self->{report_buffer_str}    = "";

  # This is where we store the matching of files this module needs
  # to compare
  $self->{compare_pair_map} = {};

  # Process constructor/init arguments
  my %params = @_;
  @$self{keys %params} = values %params;

}

# ---------------------------------------------------------------
# _report:
# Adds string to the report attribute
# ---------------------------------------------------------------

sub _report {
  my ($self, $string_to_report, $nowrap) = @_;

  unless( $nowrap ) {
    # Find extra newlines
    my($lone_newline) = ($string_to_report =~ /^\n+$/);
    my($extra_newlines) = ($string_to_report =~ /\n(\n*)$/);

    # Wrap text at a certain column
    unless($lone_newline) {
      $Text::Wrap::columns = $WRAP_COLUMN;
      $string_to_report = wrap("", "", $string_to_report) . $extra_newlines;
    }
  }

  if ( $self->{report_buffer_enable} ) {
    $self->{report_buffer_str} .= $string_to_report;
  } else {
    $self->{report_function}($string_to_report);
  }
}

# ---------------------------------------------------------------
# _report_buffer_enable:
# Enable the report buffering
# ---------------------------------------------------------------

sub _report_buffer_enable {
  my ($self) = @_;

  $self->{report_buffer_enable} = 1;
}

# ---------------------------------------------------------------
# _report_buffer_disable:
# Disable report buffering
# ---------------------------------------------------------------

sub _report_buffer_disable {
  my ($self) = @_;

  $self->{report_buffer_enable} = 0;
}

# ---------------------------------------------------------------
# _report_buffer_clear:
# Clears contents of report buffer
# ---------------------------------------------------------------

sub _report_buffer_clear {
  my ($self) = @_;

  $self->{report_buffer_str} = "";
}


# ---------------------------------------------------------------
# _report_buffer_output:
# Outputs report buffer using reporting function
# ---------------------------------------------------------------

sub _report_buffer_output {
  my ($self) = @_;

  $self->{report_function}($self->{report_buffer_str});
  $self->{report_buffer_str} = "";
}

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
