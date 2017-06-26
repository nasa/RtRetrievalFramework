#  -------------------------------------------------------------------------
#   Name:        Comparison Class Finder
#   File:        Comparison/Finder.pm
#   Created:     04/02/2004
#   Comment:
#
#   Locates the appropriate class that should be used for comparing two
#   files. This class is supplied with a list of directories to search where
#   comparison classes are derived from Comparison::Base. The class builds
#   up a list of comparison objects that contain the filenames they will
#   need to compare once instructed to do so.
#
#   Version: See ClearCase History
#
#  -------------------------------------------------------------------------

package Comparison::Finder;

use Carp;
use File::Basename;

# ===============================================================
# Public Class Methods
# ===============================================================

# ---------------------------------------------------------------
# new:
# Creates a new instance of this class
# Parameters:
# comparison_class_dirs - A references to a list of directories
# containing comparison classes.
# ---------------------------------------------------------------

sub new {
  my $classname = shift;
  my $self      = { };
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

# ===============================================================
# Public Instance Methods
# ===============================================================

# ---------------------------------------------------------------
# find_comparison_pair:
# Given the baseline and comparison filenames will add the
# filenames to an object that will compare the files.
# ---------------------------------------------------------------

sub find_comparison_pair {
  my ($self, $baseline_filename, $comparison_filename) = @_;

  # Keep track of the highest priority comparison class we have seen so far
  my $best_seen_cc_obj = undef;
  my $best_seen_cc_class = undef;
  my $best_seen_cc_pri_num = -1;

  foreach my $cc_class (@{ $self->{comparison_classes} }) {

    # Return a new object of this class if it can handle the comparison
    if($cc_class->can_compare_files($baseline_filename, $comparison_filename)) {

      my $cc_obj = undef;
      if( exists( $self->{comparison_objects}->{$cc_class} ) ) {
	$cc_obj = $self->{comparison_objects}->{$cc_class};
      } else {
	$cc_obj =
	  $cc_class->new( %{ $self->{module_options} });
      }

      # If we found only a secondary type of comparison class then
      # store that object reference away and keep looking otherwise
      # we have found a primary handler and so we return the handle
      # right away
      if($cc_class->get_priority() > $best_seen_cc_pri_num) {
	$best_seen_cc_obj = $cc_obj;
	$best_seen_cc_class = $cc_class;
	$best_seen_cc_pri_num = $cc_class->get_priority();
      }

    } # end if can compare
  } # end foreach cc file

  if (defined($best_seen_cc_obj)) {
    $best_seen_cc_obj->add_file_pair($baseline_filename, $comparison_filename);

    # Make assignment here of comparison objects so do not store objects
    # with no filenames
    $self->{comparison_objects}->{$best_seen_cc_class} = $best_seen_cc_obj;
  } else {
    die "No comparison object found for " .
      "$baseline_filename and $comparison_filename."
  }

}

# ---------------------------------------------------------------
# get_comparison_objects:
# Returns an array reference of the objects stored in the class
# that have had filenames assigned to them.
# ---------------------------------------------------------------

sub get_comparison_objects {
  my ($self) = @_;

  # Assign to an array so ref will always be to an array even if only one
  # object found
  my @comparison_objects_array = values( %{ $self->{comparison_objects} } );
  return \@comparison_objects_array;
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

  # List of classes to search for handling of comparisons
  $self->{comparison_classes} = [];

  # Map of comparison objects that have had the proper filenames
  # assigned to them
  $self->{comparison_objects} = {};

  # Process constructor/init arguments
  my %params = @_;
  @$self{keys %params} = values %params;
  if (exists($self->{comparison_class_dirs})) {

    foreach my $curr_dir (@{ $self->{comparison_class_dirs} }) {
      croak "Comparison dir does not exist: " . $curr_dir
	unless -e $curr_dir;

      push(@INC, $curr_dir);
    }

  } else {
    croak "Comparison::Finder needs to have at least the parameter " .
      "comparison_class_dirs specified.";
  }

  # Load list of modules to use for comparisons
  foreach my $classes_dir (@{ $self->{comparison_class_dirs} }) {
    opendir(DIR, $classes_dir);
    my @comp_class_files = grep { /.pm$/ } readdir(DIR);
    closedir DIR;

    foreach my $curr_cc_file (@comp_class_files) {
      my($name, $path, $suffix) = fileparse($curr_cc_file, '\..*');

      # Load the comparison class into perl space
      require $name . $suffix;

      # Store the name the finder will use now that the module has been
      # loaded
      push(@{ $self->{comparison_classes} }, $name);

    } # end for each cc file
  } # end for each cc dir

}

# ---------------------------------------------------------------
# Necessary return for modules
# ---------------------------------------------------------------

return 1;
