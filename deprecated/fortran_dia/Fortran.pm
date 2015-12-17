################################################################
# AutoDIA - Automatic Dia XML.   (C)Copyright 2001 A Trevena   #
#                                                              #
# AutoDIA comes with ABSOLUTELY NO WARRANTY; see COPYING file  #
# This is free software, and you are welcome to redistribute   #
# it under certain conditions; see COPYING file for details    #
################################################################
package Fortran;

require Exporter;

use strict;

use Data::Dumper;

use vars qw($VERSION @ISA @EXPORT);
use Autodia::Handler;

@ISA = qw(Autodia::Handler Exporter);

use Autodia::Diagram;

use File::Basename;

#---------------------------------------------------------------

#####################
# Constructor / Class Methods

# new inherited from Autodia::Handler

sub find_files_by_packagename {
    my $config = shift;
    my $args = $config->{args};
    my @filenames = ();
    die "not implemented yet, sorry\n";
}

#------------------------------------------------------------------------
# Access Methods

# parse_file inherited from Autodia::Handler

#-----------------------------------------------------------------------------
# Internal Methods

# _initialise inherited from Autodia::Handler

sub _parse {
    my $self     = shift;
    my $fh       = shift;
    my $filename = shift;
    my $Diagram  = $self->{Diagram};
    my $Class = undef;

    my $curr_operation = undef;
    my $in_module_head = undef;
    my $type_class     = undef;
    my $type_attrs     = undef;

    my $inc_attr_text = undef;

    my $in_interface  = undef;
    my $scope_text    = undef;

    my $file_basename = basename($filename);
    $file_basename =~ s/\.\w+$//;

    my $line_no = 0;
    foreach my $line (<$fh>) {
      $line_no++;
      chomp $line;
      if ($self->_discard_line($line)) {
	next;
      }

      if (defined($in_interface) && $line =~ /^\s*end\s*interface/i) {
	$in_interface = undef;
	next;
      } elsif (!defined($in_interface) && $line =~ /^\s*interface/i) {
	$in_interface = 1;
	next;
      } elsif (defined($in_interface)) {
	print "\tIgnoring interface contents: $line\n";
	next;
      }

      if (!defined($in_module_head) && $line =~ /^\s*program\s+(\w+)/i) {
	my $class_name = $1;
	$Class = $self->_add_class($Diagram, $class_name);
      } elsif (!defined($in_module_head) && $line =~ /^\s*module\s+(\w+)/i) {
	my $class_name = $1;
	$Class = $self->_add_class($Diagram, $class_name);

	$in_module_head = 1;
	next;
      } elsif ($line =~ /^\s*contains\s*/i) {
	print "\tBegin module contents\n";
	$in_module_head = undef;
	next;

      } elsif (defined($in_module_head) && !defined($inc_attr_text) && ($line =~ /^\s*public\s+/i || $line =~ /^\s*private\s+/i)) {
	print "\tIgnoring scope begin\n";
	$line =~ s/\!.*//;
	$line =~ s/&//;
	$scope_text .= $line;
	# Ignore for now
	next;

      } elsif (defined($in_module_head) && defined($scope_text) && $line =~ /\&/) {
	$line =~ s/\!.*//;
	$line =~ s/&//;
	$scope_text .= $line;

	print "\tIgnoring scope contents\n";
	next;
      } elsif (defined($in_module_head) && defined($scope_text)) {
	$scope_text = undef;
	print "\tIgnoring scope end\n";
	next;
      } elsif (defined($in_module_head) && $line =~ /^\s*type\s*\:\:\s*(\w+)/i) {
	my $type_name = $1;
	$type_class = Autodia::Diagram::Class->new(lc($type_name));
	print "\tFound type: " . $type_name . "\n";
	$type_class = $Diagram->add_class($type_class);
	$type_attrs = [];

      } elsif (defined($in_module_head) && defined($type_class) && ($line =~ /^\s*public\s*$/i || $line =~ /^\s*private\s*$/i)) {
	print "\tIgnoring class scope: $line\n";
      } elsif (defined($in_module_head) && defined($type_class) && $line =~ /^\s*end\s*type/i) {

	foreach my $curr_attr (@{$type_attrs}) {
	  print "\t\tAttr: " . $curr_attr->{'name'} . " <- " . $curr_attr->{'type'} . "\n";
	  $type_class->add_attribute($curr_attr);
	}
	print "\tEnd type\n";

	my $dependancy = Autodia::Diagram::Dependancy->new($Class, $type_class);
	
	# add dependancy to diagram
	$Diagram->add_dependancy($dependancy);
	
	# add dependancy to class
	$Class->add_dependancy($dependancy);
	
	# add dependancy to component
	$type_class->add_dependancy($dependancy);

	$type_class = undef;
	$type_attrs = undef;

      } elsif (defined($in_module_head) && $line =~ /\&/) {
	$line =~ s/\&//;
	$line =~ s/\!.*//;

	if (defined($type_class)) {
	  print "\t";
	}

	print "\tFound incomplete attribute info\n";

	if(defined($inc_attr_text)) {
	  $inc_attr_text .= $line;
	} else {
	  $inc_attr_text = $line;
	}
	next;

      } elsif (defined($in_module_head) && ( defined($inc_attr_text) || 
					     (!defined($inc_attr_text) && $line =~ /\:\:/) ||
					     ($line =~ /\s*integer\s+/) ||
					     ($line =~ /\s*double\s+precision\s+/) ||
					     ($line =~ /\s*logical\s+/)
					   )) {
	$line =~ s/\!.*//;

	if (defined($type_class)) {
	  print "\t";
	}

	print "\tFound attribute completion\n";

	my $attr_text = undef;
	if(defined($inc_attr_text)) {
	  $attr_text = $inc_attr_text . $line;
	} else {
	  $attr_text = $line;
	}

	$attr_text =~ s/\s+/ /g;
	my @line_parts;
	@line_parts = split(/\s*\:\:\s*/, $attr_text);
	if (scalar(@line_parts) != 2) {
	  @line_parts = split(/\s+/, $attr_text, 2);
	}

	my $type_def   = $line_parts[0];
	my @type_names = split(/\s*,\s*/, $line_parts[1]);


	$type_def =~ s/^\s+//;
	$type_def =~ s/\s+$//;

	for my $curr_name (@type_names) {
	  $curr_name =~ s/^\s+//;
	  $curr_name =~ s/\s+$//;

	  my $curr_attr = { name => $curr_name,
			    type => $type_def,
			  };

	  if (defined($type_class)) {
	    push(@{$type_attrs}, $curr_attr);
	  } else {
	    print "\tAttr: " . $curr_attr->{'name'} . " <- " . $curr_attr->{'type'} . "\n";
	    $Class->add_attribute($curr_attr);
	  }
	}

	$inc_attr_text = undef;
	
      } elsif (defined($in_module_head) && defined($type_class)) {
	die "Unhandled line in type declaration: $line at line $line_no\n";
      }

      if (!defined($curr_operation) && ($line =~ /^(?:\s|\w(?!end))*subroutine\s+(\w+)/i || $line =~ /^(?:\s|\w(?!end))*function\s+(\w+)/i)) {
	my $routine = $1;
	print "\tRoutine: " . $routine . "\n";

	$curr_operation = {};
	$curr_operation->{'name'} = $routine;
	$curr_operation->{'all_params'} = 0;
	$curr_operation->{'Param'} = [];
	
	$self->_get_param_names($line, $curr_operation);

	next;
      } elsif (defined($curr_operation) && !$curr_operation->{'all_params'}) {
	# Handle once we are inside of an operation
	$self->_get_param_names($line, $curr_operation);

	next;
      } elsif (defined($curr_operation) && $curr_operation->{'all_params'} && $line =~ /::\s*(\w+)/) {
	my $param_name = $1;

	my $type_str = $line;
	$type_str =~ s/::.*$//;
	my @type_parts = split(/,\s*/, $type_str);
	my $type_name = $type_parts[0];
	$type_name =~ s/^\s+//;
	$type_name =~ s/\s+$//;
	$type_name = lc($type_name);
	
	foreach my $saved_param (@{$curr_operation->{'Param'}}) {
	  if ($saved_param->{'Name'} =~ /$param_name/i) {
	    print "\t\t\t" . $saved_param->{'Name'} . " has type: " . $type_name . "\n";
	    $saved_param->{'Type'} = $type_name;
	    last;
	  }
	}

	next;
      } elsif ($line =~ /^\s*end\s*subroutine\s+(\w+)/i || $line =~ /^\s*end\s*function\s+(\w+)/i) {
	my $end_name = $1;

	unless(ref $Class) {
	  $Class = $self->_add_class($Diagram, $file_basename);
	}

	if (!defined($curr_operation->{'name'})) {
	  die "End of subroutine without properly parsed beginning: $line at line: $line_no\n";
	}

	if(lc($curr_operation->{'name'}) eq lc($end_name)) {
	  print "\tAdding: " . $curr_operation->{'name'} . "\n";
	  $Class->add_operation($curr_operation);
	  $curr_operation = undef;
	} else {
	  print STDERR "Routine found end name: " . $end_name . " does not match beginning name: " . $curr_operation->{'name'} . " for file: $filename, line: $line_no\n";
	}

	next;
      }

      # Handle dependancies to other classes and components
      if ($line =~ /^\s*use\s+(\w+)/i) {
	unless(ref $Class) {
	  $Class = $self->_add_class($Diagram, $file_basename);
	}

	my $component_name = $1;

	if (defined($curr_operation)) {
	  print "\t";
	}
	print "\tUses: " . $component_name . "\n";
	
	
	#my $Component = Autodia::Diagram::Component->new($component_name);

	# add component to diagram
	#my $exists = $Diagram->add_component($Component);
	
	# replace component if redundant
	#if (ref $exists) {
	#  $Component = $exists;
	#}

	my $Component = $self->_add_class($Diagram, $component_name);

	# create new dependancy
	my $Dependancy = Autodia::Diagram::Dependancy->new($Class, $Component);
	
	# add dependancy to diagram
	$Diagram->add_dependancy($Dependancy);
	
	# add dependancy to class
	$Class->add_dependancy($Dependancy);
	
	# add dependancy to component
	$Component->add_dependancy($Dependancy);
	next;
      }


      # create new class with name
      #$Class = Autodia::Diagram::Class->new($className);
      # add class to diagram
      # $Class = $Diagram->add_class($Class);

      # create new inheritance
      #$my $Inheritance = Autodia::Diagram::Inheritance->new($Class, $Superclass);
      # add inheritance to superclass
      #$Superclass->add_inheritance($Inheritance);
      # add inheritance to class
      #$Class->add_inheritance($Inheritance);
      # add inheritance to diagram
      #$Diagram->add_inheritance($Inheritance);

      # create new class with name
      #$Class = Autodia::Diagram::Class->new($filename);
      # add class to diagram
      #$Class = $Diagram->add_class($Class);

      #$Class->add_attribute({
      #			     name => $field,
      #			     visibility => $attribute_visibility,
      #			     Id => $Diagram->_object_count,
      #			    }) unless ($field =~ /^\$/);

      # create component
      #my $Component = Autodia::Diagram::Component->new($componentName);
      # add component to diagram
      #my $exists = $Diagram->add_component($Component);

      # create new dependancy
      #my $Dependancy = Autodia::Diagram::Dependancy->new($Class, $Component);
      # add dependancy to diagram
      #$Diagram->add_dependancy($Dependancy);
      # add dependancy to class
      #$Class->add_dependancy($Dependancy);
      # add dependancy to component
      #$Component->add_dependancy($Dependancy);

      #my $Inheritance = Autodia::Diagram::Inheritance->new($Class, $Superclass);
      # add inheritance to superclass
      #$Superclass->add_inheritance($Inheritance);
      # add inheritance to class
      #$Class->add_inheritance($Inheritance);
      # add inheritance to diagram
      #$Diagram->add_inheritance($Inheritance);

      # add accessor
      #$Class->add_operation({ name => $col, visibility => $visibility, Id => $Diagram->_object_count() } );
      #$Class->add_operation({ name => $col, visibility => $visibility, Id => $Diagram->_object_count() } );
    }

    $self->{Diagram} = $Diagram;
    close $fh;
    return;
}

sub _add_class {
  my $self       = shift;
  my $diagram    = shift;
  my $class_name = shift;

  my $class_obj = Autodia::Diagram::Class->new(lc($class_name));
  print "Found class: " . $class_name . "\n";

  $class_obj = $diagram->add_class($class_obj);

  return $class_obj
}

sub _discard_line
{
  my $self    = shift;
  my $line    = shift;
  my $discard = 0;

  SWITCH:
    {
	if ($line =~ m/^\s*$/) # if line is blank or white space discard
	{
	    $discard = 1;
	    last SWITCH;
	}

	if ($line =~ /^\s*\!/) # if line is a comment discard
	{
	    $discard = 1;
	    last SWITCH;
	}

    }
    return $discard;
}

sub _get_param_names {
  my $self = shift;
  my $line = shift;
  my $curr_operation = shift;

  my $param_str = "";
  if ($line =~ /\(([^)]+)\)(?:\s|\!)?\s*\!?/) {
    $param_str = $1;
    $curr_operation->{'all_params'} = 1;
  } elsif ($line =~ /\(([^!]+)\&/) {
    $param_str = $1;
  } elsif ($line =~ /\s*([^)]+)\)/) {
    $curr_operation->{'all_params'} = 1;
    $param_str = $1;
  } else {
    $line =~ s/&//;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    $param_str = $line;
  }
  my @param_names = split(/\s*,\s*/, $param_str);

  print "\t\tParams: ";	
  my $pc = 0;
  foreach my $name (@param_names) {
    print "$name ";

    push(@{$curr_operation->{'Param'}},
	 {
	  Name => $name,
	 });
    $pc++;
  }
  print "\n";

}

####-----

sub _is_package
  {
    my $self    = shift;
    my $package = shift;
    my $Diagram = $self->{Diagram};

    unless(ref $$package)
       {
	 my $filename = shift;
	 # create new class with name
	 $$package = Autodia::Diagram::Class->new($filename);
	 # add class to diagram
	 $Diagram->add_class($$package);
       }

    return;
  }



####-----

1;

