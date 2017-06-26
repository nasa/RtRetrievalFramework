#!/usr/bin/env perl
#  -------------------------------------------------------------------------
#   Name:        PGE Regession Output Comparison Program
#   File:        regression_compare.pl
#   Author:      James McDuffie
#
#   Description:
#
#   Compares the file in a baseline and comparison directory. Figures out
#   how to compare each file based upon the type of the file.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

use strict;

# Initialize location of modules
use File::Basename;
my($run_name, $run_path, $run_suffix);

BEGIN {
  ($run_name, $run_path, $run_suffix) = fileparse($0, '\..*');
  $run_path =~ s/\/$//; # Remove trailing slash
  push(@INC, $run_path . "/Perl");
}

# Modules from Perl distribution
use Getopt::Long;
use File::Find;
use Carp;
use IO::Handle;

# For debugging status information
use POSIX;
use Term::Cap;

# Custom modules
use Comparison::Finder;

# Turn on autoflushing for output streams
STDOUT->autoflush(1);
STDERR->autoflush(1);

# ---------------------------------------------------------------
# Setup variables used to store command line parameters.
# ---------------------------------------------------------------

# Use hash to store flags.

my %flags =
  ( help                 => undef,
    output_report        => undef,
    check_all            => undef,
    verbose              => undef,
    inspect              => undef,
    file_append          => undef,
    debug                => undef,
    list_mapping         => undef,
  );

# Use hash to store parameters with values.

my %params =
  ( baseline_directory           => "",
    comparison_directory         => "",
    comparison_classes_directory => $run_path . "/ComparisonClasses",
    report_file                  => "",
    module_options               => "",
    find_ignore                  => [] );

# ---------------------------------------------------------------
# Process command line arguments with Getopt module.
# ---------------------------------------------------------------

my $get_opt_results =
  GetOptions( 'help'             => \$flags{help},
	      'output_report'    => \$flags{output_report},
	      'all'              => \$flags{check_all},
	      'verbose'          => \$flags{verbose},
	      'inspect'          => \$flags{inspect},
	      'file_append'      => \$flags{file_append},
	      'debug'            => \$flags{debug},
	      'list_mapping'     => \$flags{list_mapping},
	      'baseline_dir=s'   => \$params{baseline_directory},
	      'comparison_dir=s' => \$params{comparison_directory},
	      'dir_for_cc=s'     => \$params{comparison_classes_directory},
	      'report_file=s'    => \$params{report_file},
	      'module_options=s' => \$params{module_options},
	      'find_ignore=s'    => \@{ $params{find_ignore} } );

if ( !$get_opt_results ) {

    print_help() and exit(1);
}

# Supply help if requested.

if( $flags{help} ) {
    print_help() and exit(1);
}

if( $flags{inspect} && $flags{check_all} ) {
  print STDERR "The -inspect and -all options can not be specified together.";
  print_help() and exit(1);
}

# Enforce required parameters and flags.

if( !$params{baseline_directory}  or
    !$params{comparison_directory} ) {

    print STDERR "Required Parameter(s) missing.\n";
    print_help() and exit(1);
}

# ---------------------------------------------------------------
# End options setup code
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Variables used for finding files in subdirectories
# ---------------------------------------------------------------

my %FIND_OPTIONS = ( wanted => \&wanted, no_chdir => 1, follow => 1 );
my $current_files_array;

# ---------------------------------------------------------------
# File handle for report file if one is specified
# ---------------------------------------------------------------

*REPORT_FILE;

# ---------------------------------------------------------------
# Set error code so that if the program dies it will return that
# comparison failed
# ---------------------------------------------------------------

$! = 1;

# ---------------------------------------------------------------
# Initialize the TermCap to control screen output
# ---------------------------------------------------------------

my $termcap = init_termcap();

# Maximum screen sizes for tcap
my ($max_screen_row, $max_screen_col) =
  ($termcap->{_li} - 1, $termcap->{_co} - 1);


# ---------------------------------------------------------------
# Finally call main method to run program
# ---------------------------------------------------------------

main();

# ===============================================================
# Subroutines
# ===============================================================

# ---------------------------------------------------------------
# print_help:
# Prints information describing command line options
# ---------------------------------------------------------------

sub print_help {

  print <<HELP;
Compares Baseline and Comparison Directories Using Comparison Modules

Matches filenames in the baseline and comparison directories and then
performs comparisons on those files using comparison modules that 
understand the file format of various files.

Required arguments:
   -baseline_dir "baseline_dir/"
      Directory of baseline output
   -comparison_dir "comparison_dir/"
      Directory of comparison output

Optional arguments:
   -output_report
      Output verbose info reported by modules to the screen. Without this
      option only a return code is returned. Use the -report_file option
      to save output to a file without outputting to the screen.
   -all
      Check all values of files instead of quiting after the first difference
   -list_mapping
      Instead of comparing files, simply output a text filename mapping
   -verbose
      Turns on verbose output which may produce additional information from
      various modules
   -inspect_idl
      This option will cause IDL comparison routines to be run in a way such
      that inspection is possible when an error is found. Not compatible with
      the -all option.
   -debug
      Runs the comparison module in debug mode that runs interactive IDL for
      each set of files compared.

   -report_file "report.log"
      If the option is specified verbose information is stored in this file
   -file_append
      Appends to the report file instead of overwriting
   -dir_for_cc "comp_class_dir/"
      Directory where the comparison classes other than those distributed with
      the program are located.
   -module_options "option1_name=option1_value option2_name=option_2=value"
      A list of strings separated by spaces that specify options passed to
      comparison modules. Not all comparison modules use these options.
      Ex: -module_options "do_db_check=0"
   -find_ignore "regular_expression"
      Regular expression that allows files to be ignored when searching
      for files to compare.
      Ex: -find_ignore '.*\.txt'
HELP

} # end of sub print_help

# ---------------------------------------------------------------
# display:
# Sends string to the appropriate place, either to a file or to
# be printed on standard output. Display must be declared before
# all other subroutines since it has a prototype that is used to
# allow calling it like this:
# display "Some string\n";
# ---------------------------------------------------------------

sub display (@) {
  if($flags{output_report}) {
    print @_;
  }
  if($params{report_file}) {
    print REPORT_FILE @_;
  }
}

# ---------------------------------------------------------------
# init_termcap:
# Initializes the term cap screen control library
# ---------------------------------------------------------------

sub init_termcap {
  my $delay = (shift() || 0) * 0.005;
  my $termios = POSIX::Termios->new();
  $termios->getattr;
  my $ospeed = $termios->getospeed;

  unless ( defined($ENV{TERM}) ) {
    $ENV{TERM} = "vt100";
  }

  my $tcap = Term::Cap->Tgetent ({ TERM => undef, OSPEED => $ospeed });
  $tcap->Trequire(qw(cl cm cd));

  return $tcap;
}

# ---------------------------------------------------------------
# main:
# Starts off he whole big shebang
# ---------------------------------------------------------------

sub main {
  my $base_dir = $params{'baseline_directory'};
  my $comp_dir = $params{'comparison_directory'};

  if ( $flags{verbose} ) {
    $termcap->Tputs('cl', 1, *STDOUT);
    $termcap->Tgoto('cm', 0, $max_screen_row, *STDOUT);
  }

  display "Comparison program beginning.\n";

  # Strip leading and trailing unnecessaries from directory filenames
  $base_dir =~ s/^\.\///;
  $comp_dir =~ s/^\.\///;

  $base_dir =~ s/\/$//;
  $comp_dir =~ s/\/$//;

  unless(-e $base_dir) {
    print STDERR "The directory $base_dir does not exist!\n\n";
    print_help() and exit(1);
  }

  unless(-e $comp_dir) {
    print STDERR "The directory $comp_dir does not exist!\n\n";
    print_help() and exit(1);
  }

  # Open report file if this option is specified
  if ($params{report_file}) {
    if ($flags{file_append}) {
      open(REPORT_FILE, ">>" . $params{report_file});
    } else {
      open(REPORT_FILE, ">" . $params{report_file});
    }
  }

  # Retrieve list of files under both the baseline and comparison directories
  display "Finding files to compare.\n";
  my ($base_files, $comp_files) = find_files($base_dir, $comp_dir);

  # Find the comparison objects to be used for comparisons, in the process a map
  # is created that matches the files under the baseline directory to those
  # in the comparison directory. Also retrieve a list of files only found in
  # either directory
  display "Finding comparison objects for found files.\n";
  my ($comp_objects, $match_map, $originals_list) =
    find_comparison_objects( $base_files, $comp_files,
			     $base_dir, $comp_dir,
			     \%params, \%flags );

  if( $flags{list_mapping} ) {
    foreach my $key (keys(%{$match_map})) {
      my $base_filename = $key;
      my $comp_filename = $match_map->{$key};
      if(-f $base_filename && -f $comp_filename) {
	display $base_filename . "\t" . $comp_filename . "\n";
      }
    }
    exit;
  }

  # Log list of original files in the way this script is configured
  log_originals($originals_list);

  display "Beginning file comparison...\n";
  my $result = compare_files( $comp_objects );

  # Close report file if this option is specified
  if($params{report_file}) {
    close(REPORT_FILE);
  }

  exit($result);

}

# ---------------------------------------------------------------
# find_files:
# Starts at the baseline and comparison directories and searches
# for files matched by the wanted function. Returns array 
# references to lists of files in the two directories
# ---------------------------------------------------------------

sub find_files {
  my($base_dir, $comp_dir) = @_;

  print "Finding in baseline directory:\n\n" if $flags{verbose};

  my @base_files;
  $current_files_array = \@base_files;
  find(\%FIND_OPTIONS, $base_dir);
  @base_files = sort(@base_files);

  print "\nFinding in comparison directory:\n\n" if $flags{verbose};

  my @comp_files;
  $current_files_array = \@comp_files;
  find(\%FIND_OPTIONS, $comp_dir);
  @comp_files = sort(@comp_files);

  print "\n" if $flags{verbose};

  return (\@base_files, \@comp_files);
}

# ---------------------------------------------------------------
# wanted:
# Called by File::Find for each file it looks at
# ---------------------------------------------------------------

sub wanted {

  # Ignore files that match the find_ignore regexp
  for my $ignore_name ( @{ $params{find_ignore} } ) {
    if($params{find_ignore} && $File::Find::name =~ /$ignore_name/) {
      return;
    }
  }

  push(@{ $current_files_array }, $File::Find::name);

  if ($flags{verbose}) {
      $termcap->Tgoto('cm', 0, $max_screen_row, *STDOUT);
      print scalar(@{ $current_files_array });
    }
}

# ---------------------------------------------------------------
# find_comparison_objects:
# Looks through two arrays containing filenames one for the
# baseline directory and one for the comparison directory. Tries
# to match files from one list with the other.
# Uses the finder to obtain objects to be used for the comparisons.
# ---------------------------------------------------------------

sub find_comparison_objects {
  my($base_files, $comp_files, $base_dir, $comp_dir, $params, $flags) = @_;

  # Create object that will let us get references to comparison classes
  my $comp_classes_dir = $params->{comparison_classes_directory};

  # Parse command line options and pass as regular options
  my @command_option_parts = split(/[ \,]/, $params->{module_options});
  my %command_line_options =
    map { my($k, $v) = split(/[=]/);
	  $v = 1 unless defined($v);
	  $k => $v } @command_option_parts;

  # Create finder object to locate modules to perform the comparisons
  my $comp_finder =
    Comparison::Finder->new( comparison_class_dirs => [ $comp_classes_dir ],
			     module_options =>
			     { regression_base_dir  => $run_path,
			       baseline_directory   => $base_dir,
			       comparison_directory => $comp_dir,
			       check_all            => $flags->{check_all},
			       verbose              => $flags->{verbose},
			       inspect              => $flags->{inspect},
			       debug                => $flags->{debug},
			       report_function      => \&display,
			       %command_line_options,
			     }
			   );

  # Filenames will be sorted into these two variables
  my %match_map;
  my @originals_list;

  # Hash comparison files basenames for easy comparison to original files
  my %comp_basename_hash;
  foreach my $comp_filename (@{ $comp_files }) {
    # Get a unique name that may include some path information
    my $comp_rel_filename = get_relative_unique_name($comp_filename, $comp_dir);
    $comp_basename_hash{$comp_rel_filename} = $comp_filename;
  }

  # New line so that termcap stuff does not overwrite last line
  print "\n" if $flags{verbose};

  # Now compare each orig filename with the comp filenames
  my $file_count = 0;
  foreach my $base_filename (@{ $base_files }) {
    my $base_rel_filename = get_relative_unique_name($base_filename, $base_dir);

    if( exists($comp_basename_hash{$base_rel_filename}) ) {
      my $comp_filename = $comp_basename_hash{$base_rel_filename};

      # Save filename match and mark as found in comparison dir hash for 
      # later originals check
      $match_map{$base_filename} = $comp_filename;
      $comp_basename_hash{$base_rel_filename} = "found";

      # Print a debug message for times when number of files to
      # find comparison objects for is huge
      if ($flags{verbose}) {
	$termcap->Tgoto('cm', 0, $max_screen_row, *STDOUT);
	print ++$file_count;
      }

      # Store the pair into the appropriate finder object
      if(-f $base_filename && -f $comp_filename) {
	$comp_finder->find_comparison_pair($base_filename, $comp_filename);
      }

    } elsif ($base_filename ne $base_dir) {
      # The file in the list of base files is original to the baseline dir
      push(@originals_list, $base_filename);
    }
  }

  print "\n" if $flags{verbose};

  # Now look bach through the comp_basename_hash to find the files
  # that only exist in the comparison directory
  foreach my $comp_basename ( keys(%comp_basename_hash) ) {
    unless( $comp_basename_hash{$comp_basename} eq "found" ||
	    $comp_basename_hash{$comp_basename} eq $comp_dir ) {
      push(@originals_list, $comp_basename_hash{$comp_basename});
    }
  }

  return ( $comp_finder->get_comparison_objects(), \%match_map, \@originals_list );
}

# ---------------------------------------------------------------
# get_relative_unique_name:
# Retrieves the relative unique part of a filename
# ---------------------------------------------------------------

sub get_relative_unique_name {
  my($input_filename, $input_dir) = @_;

  my $output_filename = $input_filename;
  $output_filename =~ s/$input_dir\/+//;

  return $output_filename;
}

# ---------------------------------------------------------------
# log_originals:
#
# Outputs to the log the names of the files that are not located 
# in both directories
# ---------------------------------------------------------------

sub log_originals {
  my($originals_list) = @_;

  # Log list of original files in the way this script is configured
  if(scalar(@{$originals_list}) > 0) {
    display "-" x 80, "\n";

    # For now just print out the list of original files
    display "Files not in both directories:\n";
    foreach my $originals_filename ( @{ $originals_list } ) {
      display $originals_filename . "\n";
    }

    # Fail if there are original files in either directory
    unless($flags{check_all} || $flags{inspect}) {
      # Close report file if this option is specified
      if($params{report_file}) {
	close(REPORT_FILE);
      }

      exit(1);
    }
  } # end if have originals

}

# ---------------------------------------------------------------
# compare_files:
# Iterates through an array of comparison module objects and
# instructs each to compare its files. Keeps track of the total
# number of different files and total compared.
# ---------------------------------------------------------------

sub compare_files {
  my ($comparison_objects) = @_;

  my $number_diff_files = 0;
  my $number_compared_files = 0;
  foreach my $c_obj ( @{ $comparison_objects } ) {

    my ($diff, $compared) = $c_obj->compare();

    # If none compared then there has been an error so all the files 
    # are considered bogus
    if($compared == 0) {
      $number_diff_files = $c_obj->num_file_pairs();
    } else {
      $number_diff_files += $diff;
      $number_compared_files += $compared;
    }
  } # end of foreach match_map key

  display "-" x 80, "\n";
  display "There are " . $number_diff_files .
    " different files of " . $number_compared_files . " compared\n";

  if($number_diff_files > 0) {
    return 1;
  } else {
    return 0;
  }
}
