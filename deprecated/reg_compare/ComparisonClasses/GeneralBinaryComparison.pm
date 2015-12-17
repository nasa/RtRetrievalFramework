#  -------------------------------------------------------------------------
#   Name:        General Binary Comparison
#   File:        GeneralBinaryComparison.pm
#   Created:     04/02/04
#   Comment:
#
#   Performs a general binary comparison that takes into consideration a
#   possible ASCII header.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package GeneralBinaryComparison;

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

  # If the filename has .bin in the filename then most likely a
  # binary file!
  return 1 if($baseline_filename =~ /.bin$/o);

  # Or if Perl tells us the file is binary file then it probably is
  return (-B $baseline_filename && -B $comparison_filename);
}

# ---------------------------------------------------------------
# get_priority:
# This class should have higher priority than the general ascii
# comparison routine but not more than any other specialized
# comparison routine.
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
  # If -B returns binary for the files then that probably means they
  # have no ASCII header
  unless(-B $base_filename && -B $comp_filename) {
    my $num_header_diff =
      $self->_compare_headers($base_filename, $comp_filename);
  }

  # Compare binary portion of file
  my $has_differences =
    $self->_compare_binary_contents($base_filename, $comp_filename);

  # Report on the number of differences if we are indeed checking everything
  if($self->{check_all}) {
      $self->_report("Files have " . $self->{NUM_BIN_DIFFERENCES} . " differences\n");
  }

  return $has_differences;
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
  $self->{BASELINE_HEADER_SIZE}   = 0;
  $self->{COMPARISON_HEADER_SIZE} = 0;
  $self->{NUM_ASCII_DIFFERENCES}  = 0;
  $self->{NUM_BIN_DIFFERENCES}    = 0;

  # Default values for arguments
  $self->{frac_difference_threshold} = 1.0e-8;
  $self->{byte_increment} = 8;
  $self->{pack_template} = "d";
  $self->{header_strip_list} = undef;

}

# ---------------------------------------------------------------
# _compare_binary_contents:
# Compares and reports on the binary portions of the files.
# Returns the number of differences.
# ---------------------------------------------------------------

sub _compare_binary_contents {
  my ($self, $base_filename, $comp_filename) = @_;

  my ($base_contents, $base_size) =
    $self->_read_file_bytes( $base_filename,
			     $self->{BASELINE_HEADER_SIZE} );

  my ($comp_contents, $comp_size) = 
    $self->_read_file_bytes( $comp_filename,
			     $self->{COMPARISON_HEADER_SIZE} );

  # If the files have different sizes then they are necessarily different!
  if($base_size != $comp_size) {
    $self->_report("Binary portion of files differ in size!\n");
    return 1;
  }

  # Report how many items were looked at
  $self->_report("Comparing $base_size elements.\n");

  # Compare all items
  for(my $count = 0; $count < $base_size; $count++) {

    my $base_value = $self->_convert_value($base_contents->[$count]);
    my $comp_value = $self->_convert_value($comp_contents->[$count]);

    # Just incase convert returns undef values causing a divide by zero
    unless(defined($base_value) && defined($comp_value)) {
      $self->_report("WARNING: Convert in GeneralBinaryComparison returned undefined values!\n");
      next;
    }

    # Compute the fractional difference between the two values
    my $abs_diff = abs($comp_value - $base_value);

    # Add a small amount to the denominators to keep from dividing by 0
    my $frac_diff = 0;
    if($base_value > $comp_value) {
      $frac_diff = $abs_diff / ($base_value + 1.0e-20);
    } else {
      $frac_diff = $abs_diff / ($comp_value + 1.0e-20);
    }

    if($frac_diff > $self->{frac_difference_threshold}) {
      # These differences are above our threshold so output them
      $self->_report("  baseline_value[$count] = $base_value\n");
      $self->_report("comparison_value[$count] = $comp_value\n");
      $self->_report("Fractional Difference = " . $frac_diff . "\n");
      $self->_report("\n");

      $self->{NUM_BIN_DIFFERENCES}++;

      unless($self->{check_all}) {
	return 1;
      }
    }

  }

  if($self->{NUM_BIN_DIFFERENCES} > 0) {
    return 1;
  } else {
    return undef;
  }
}

# ---------------------------------------------------------------
# _compare_headers:
# Compares the headers of files, store the header sizes as
# attributes and turn return the number of differences.
# ---------------------------------------------------------------

sub _compare_headers {
  my ($self, $base_filename, $comp_filename) = @_;

  my($base_header, $comp_header);

  ($base_header, $self->{BASELINE_HEADER_SIZE}) =
    $self->_read_header($base_filename);

  ($comp_header, $self->{COMPARISON_HEADER_SIZE}) =
    $self->_read_header($comp_filename);

  my $base_header_lines = scalar(@{$base_header});
  my $comp_header_lines = scalar(@{$comp_header});

  if($base_header_lines != $comp_header_lines) {
    $self->_report("Number of header lines unequal:\n");
    $self->_report("Baseline: $base_header_lines Comparison: $comp_header_lines\n");
  }

  my $max_lines;
  if($base_header_lines > $comp_header_lines) {
    $max_lines = $base_header_lines;
  } else {
    $max_lines = $comp_header_lines;
  }

  for(my $header_line = 0; $header_line < $max_lines; $header_line++) {
    my $base_line = $base_header->[$header_line];
    my $comp_line = $comp_header->[$header_line];

    chomp($base_line);
    chomp($comp_line);

    # Ignore things that we dont care if they are different in the header
    my $skip_current_line = 0;
    my $output_base_line = 1;
    my $output_comp_line = 1;
    foreach my $strip_string (@{ $self->{header_strip_list} }) {
      if($base_line =~ /^$strip_string/ && $comp_line =~ /^$strip_string/) {
	$skip_current_line = 1;
	last;
      } elsif($base_line =~ /^$strip_string/) {
	$output_base_line = 0;
	last;
      } elsif($comp_line =~ /^$strip_string/) {
	$output_comp_line = 0;
	last;
      }
    }

    if($skip_current_line) {
      next;
    }

    if($base_line ne $comp_line) {
      $self->_report("Line: " . ($header_line + 1) . "\n");
      if($output_base_line) {
	$self->_report("< $base_line\n");
      }
      if($output_comp_line) {
	$self->_report("> $comp_line\n");
      }

      $self->{NUM_ASCII_DIFFERENCES}++;
    }
  }

  if($self->{NUM_ASCII_DIFFERENCES} > 0) {
    return 1;
  } else {
    return undef;
  }
}

# ---------------------------------------------------------------
# _read_file_bytes:
# Reads bytes from a files and converts those bytes into an array
# of numbers. The number of bytes to read every increment
# determines how the number is read.
# ---------------------------------------------------------------

sub _read_file_bytes {
  my $self = shift;
  my($filename, $offset) = @_;

  open(FILE, "$filename");
  binmode(FILE);
  seek(FILE, $offset, "0");

  my $num_read;
  my $array_index = 0;

  my @contents_array;
  do {
    my $temp_val;

    $num_read = read(FILE, $temp_val, $self->{byte_increment});

    if($num_read > 0) {
      $contents_array[$array_index] = $temp_val;
      $array_index++;
    }
  } while($num_read > 0);

  close(FILE);

  return (\@contents_array, $array_index);
}

# ---------------------------------------------------------------
# _convert_value:
# Converts binary data into a variable understandable to perl
# ---------------------------------------------------------------

sub _convert_value {
  my $self = shift;
  my($value) = @_;

  return unpack $self->{pack_template}, $value;
}


# ---------------------------------------------------------------
# _read_header:
# Reads the ASCII header of a binary file looking for the start
# of the binary data.
# Parameter:
#   filename of binary file to read header from
# Return:
#   reference to array containing header
#   size of header in bytes
# ---------------------------------------------------------------

sub _read_header {
  my $self = shift;
  my($filename) = @_;

  open(BIN_FILE, $filename);

  # Location where we will place the contents of the header
  my @header_array;

  my $current_line;
  my $header_size;
  do {
    $current_line = <BIN_FILE>;
    $header_size += length($current_line);
    push(@header_array, $current_line);
  } while (!eof(BIN_FILE) && !($current_line =~ /^End_of_Header/));

  # We saw an end of header mark, but there might be 0 to 2 more lines
  # of information that we want to discard this is usually column
  # header information. So search for a possibly two lines beginning with
  # alphabetic info

  for(my $ah_count; $ah_count < 2; $ah_count++) {
    $current_line = <BIN_FILE>;

    # If this line starts with a alphanumeric character then
    # account for its length in the header size. If not then we quit
    # looking
    my $chomped_line = $current_line;
    chomp($chomped_line);
    $chomped_line =~ s/^\s+//;

    if($chomped_line =~ /^([\w-\/]+\s*)+$/ && length($chomped_line) > 3) {
      $header_size += length($current_line);
      push(@header_array, $current_line);
    } else {
      last;
    }
  }

  close(BIN_FILE);

  return (\@header_array, $header_size);
}

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
