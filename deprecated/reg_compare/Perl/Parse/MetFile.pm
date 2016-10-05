####################################################################
#
#    This file was generated using Parse::Yapp version 1.05.
#
#        Don't edit this file, use source file instead.
#
#             ANY CHANGE MADE HERE WILL BE LOST !
#
####################################################################
package Parse::MetFile;
use vars qw ( @ISA );
use strict;

@ISA= qw ( Parse::Yapp::Driver );
#Included Parse/Yapp/Driver.pm file----------------------------------------
{
#
# Module Parse::Yapp::Driver
#
# This module is part of the Parse::Yapp package available on your
# nearest CPAN
#
# Any use of this module in a standalone parser make the included
# text under the same copyright as the Parse::Yapp module itself.
#
# This notice should remain unchanged.
#
# (c) Copyright 1998-2001 Francois Desarmenien, all rights reserved.
# (see the pod text in Parse::Yapp module for use and distribution rights)
#

package Parse::Yapp::Driver;

require 5.004;

use strict;

use vars qw ( $VERSION $COMPATIBLE $FILENAME );

$VERSION = '1.05';
$COMPATIBLE = '0.07';
$FILENAME=__FILE__;

use Carp;

#Known parameters, all starting with YY (leading YY will be discarded)
my(%params)=(YYLEX => 'CODE', 'YYERROR' => 'CODE', YYVERSION => '',
			 YYRULES => 'ARRAY', YYSTATES => 'ARRAY', YYDEBUG => '');
#Mandatory parameters
my(@params)=('LEX','RULES','STATES');

sub new {
    my($class)=shift;
	my($errst,$nberr,$token,$value,$check,$dotpos);
    my($self)={ ERROR => \&_Error,
				ERRST => \$errst,
                NBERR => \$nberr,
				TOKEN => \$token,
				VALUE => \$value,
				DOTPOS => \$dotpos,
				STACK => [],
				DEBUG => 0,
				CHECK => \$check };

	_CheckParams( [], \%params, \@_, $self );

		exists($$self{VERSION})
	and	$$self{VERSION} < $COMPATIBLE
	and	croak "Yapp driver version $VERSION ".
			  "incompatible with version $$self{VERSION}:\n".
			  "Please recompile parser module.";

        ref($class)
    and $class=ref($class);

    bless($self,$class);
}

sub YYParse {
    my($self)=shift;
    my($retval);

	_CheckParams( \@params, \%params, \@_, $self );

	if($$self{DEBUG}) {
		_DBLoad();
		$retval = eval '$self->_DBParse()';#Do not create stab entry on compile
        $@ and die $@;
	}
	else {
		$retval = $self->_Parse();
	}
    $retval
}

sub YYData {
	my($self)=shift;

		exists($$self{USER})
	or	$$self{USER}={};

	$$self{USER};
	
}

sub YYErrok {
	my($self)=shift;

	${$$self{ERRST}}=0;
    undef;
}

sub YYNberr {
	my($self)=shift;

	${$$self{NBERR}};
}

sub YYRecovering {
	my($self)=shift;

	${$$self{ERRST}} != 0;
}

sub YYAbort {
	my($self)=shift;

	${$$self{CHECK}}='ABORT';
    undef;
}

sub YYAccept {
	my($self)=shift;

	${$$self{CHECK}}='ACCEPT';
    undef;
}

sub YYError {
	my($self)=shift;

	${$$self{CHECK}}='ERROR';
    undef;
}

sub YYSemval {
	my($self)=shift;
	my($index)= $_[0] - ${$$self{DOTPOS}} - 1;

		$index < 0
	and	-$index <= @{$$self{STACK}}
	and	return $$self{STACK}[$index][1];

	undef;	#Invalid index
}

sub YYCurtok {
	my($self)=shift;

        @_
    and ${$$self{TOKEN}}=$_[0];
    ${$$self{TOKEN}};
}

sub YYCurval {
	my($self)=shift;

        @_
    and ${$$self{VALUE}}=$_[0];
    ${$$self{VALUE}};
}

sub YYExpect {
    my($self)=shift;

    keys %{$self->{STATES}[$self->{STACK}[-1][0]]{ACTIONS}}
}

sub YYLexer {
    my($self)=shift;

	$$self{LEX};
}


#################
# Private stuff #
#################


sub _CheckParams {
	my($mandatory,$checklist,$inarray,$outhash)=@_;
	my($prm,$value);
	my($prmlst)={};

	while(($prm,$value)=splice(@$inarray,0,2)) {
        $prm=uc($prm);
			exists($$checklist{$prm})
		or	croak("Unknow parameter '$prm'");
			ref($value) eq $$checklist{$prm}
		or	croak("Invalid value for parameter '$prm'");
        $prm=unpack('@2A*',$prm);
		$$outhash{$prm}=$value;
	}
	for (@$mandatory) {
			exists($$outhash{$_})
		or	croak("Missing mandatory parameter '".lc($_)."'");
	}
}

sub _Error {
	print "Parse error.\n";
}

sub _DBLoad {
	{
		no strict 'refs';

			exists(${__PACKAGE__.'::'}{_DBParse})#Already loaded ?
		and	return;
	}
	my($fname)=__FILE__;
	my(@drv);
	open(DRV,"<$fname") or die "Report this as a BUG: Cannot open $fname";
	while(<DRV>) {
                	/^\s*sub\s+_Parse\s*{\s*$/ .. /^\s*}\s*#\s*_Parse\s*$/
        	and     do {
                	s/^#DBG>//;
                	push(@drv,$_);
        	}
	}
	close(DRV);

	$drv[0]=~s/_P/_DBP/;
	eval join('',@drv);
}

#Note that for loading debugging version of the driver,
#this file will be parsed from 'sub _Parse' up to '}#_Parse' inclusive.
#So, DO NOT remove comment at end of sub !!!
sub _Parse {
    my($self)=shift;

	my($rules,$states,$lex,$error)
     = @$self{ 'RULES', 'STATES', 'LEX', 'ERROR' };
	my($errstatus,$nberror,$token,$value,$stack,$check,$dotpos)
     = @$self{ 'ERRST', 'NBERR', 'TOKEN', 'VALUE', 'STACK', 'CHECK', 'DOTPOS' };

#DBG>	my($debug)=$$self{DEBUG};
#DBG>	my($dbgerror)=0;

#DBG>	my($ShowCurToken) = sub {
#DBG>		my($tok)='>';
#DBG>		for (split('',$$token)) {
#DBG>			$tok.=		(ord($_) < 32 or ord($_) > 126)
#DBG>					?	sprintf('<%02X>',ord($_))
#DBG>					:	$_;
#DBG>		}
#DBG>		$tok.='<';
#DBG>	};

	$$errstatus=0;
	$$nberror=0;
	($$token,$$value)=(undef,undef);
	@$stack=( [ 0, undef ] );
	$$check='';

    while(1) {
        my($actions,$act,$stateno);

        $stateno=$$stack[-1][0];
        $actions=$$states[$stateno];

#DBG>	print STDERR ('-' x 40),"\n";
#DBG>		$debug & 0x2
#DBG>	and	print STDERR "In state $stateno:\n";
#DBG>		$debug & 0x08
#DBG>	and	print STDERR "Stack:[".
#DBG>					 join(',',map { $$_[0] } @$stack).
#DBG>					 "]\n";


        if  (exists($$actions{ACTIONS})) {

				defined($$token)
            or	do {
				($$token,$$value)=&$lex($self);
#DBG>				$debug & 0x01
#DBG>			and	print STDERR "Need token. Got ".&$ShowCurToken."\n";
			};

            $act=   exists($$actions{ACTIONS}{$$token})
                    ?   $$actions{ACTIONS}{$$token}
                    :   exists($$actions{DEFAULT})
                        ?   $$actions{DEFAULT}
                        :   undef;
        }
        else {
            $act=$$actions{DEFAULT};
#DBG>			$debug & 0x01
#DBG>		and	print STDERR "Don't need token.\n";
        }

            defined($act)
        and do {

                $act > 0
            and do {        #shift

#DBG>				$debug & 0x04
#DBG>			and	print STDERR "Shift and go to state $act.\n";

					$$errstatus
				and	do {
					--$$errstatus;

#DBG>					$debug & 0x10
#DBG>				and	$dbgerror
#DBG>				and	$$errstatus == 0
#DBG>				and	do {
#DBG>					print STDERR "**End of Error recovery.\n";
#DBG>					$dbgerror=0;
#DBG>				};
				};


                push(@$stack,[ $act, $$value ]);

					$$token ne ''	#Don't eat the eof
				and	$$token=$$value=undef;
                next;
            };

            #reduce
            my($lhs,$len,$code,@sempar,$semval);
            ($lhs,$len,$code)=@{$$rules[-$act]};

#DBG>			$debug & 0x04
#DBG>		and	$act
#DBG>		and	print STDERR "Reduce using rule ".-$act." ($lhs,$len): ";

                $act
            or  $self->YYAccept();

            $$dotpos=$len;

                unpack('A1',$lhs) eq '@'    #In line rule
            and do {
                    $lhs =~ /^\@[0-9]+\-([0-9]+)$/
                or  die "In line rule name '$lhs' ill formed: ".
                        "report it as a BUG.\n";
                $$dotpos = $1;
            };

            @sempar =       $$dotpos
                        ?   map { $$_[1] } @$stack[ -$$dotpos .. -1 ]
                        :   ();

            $semval = $code ? &$code( $self, @sempar )
                            : @sempar ? $sempar[0] : undef;

            splice(@$stack,-$len,$len);

                $$check eq 'ACCEPT'
            and do {

#DBG>			$debug & 0x04
#DBG>		and	print STDERR "Accept.\n";

				return($semval);
			};

                $$check eq 'ABORT'
            and	do {

#DBG>			$debug & 0x04
#DBG>		and	print STDERR "Abort.\n";

				return(undef);

			};

#DBG>			$debug & 0x04
#DBG>		and	print STDERR "Back to state $$stack[-1][0], then ";

                $$check eq 'ERROR'
            or  do {
#DBG>				$debug & 0x04
#DBG>			and	print STDERR 
#DBG>				    "go to state $$states[$$stack[-1][0]]{GOTOS}{$lhs}.\n";

#DBG>				$debug & 0x10
#DBG>			and	$dbgerror
#DBG>			and	$$errstatus == 0
#DBG>			and	do {
#DBG>				print STDERR "**End of Error recovery.\n";
#DBG>				$dbgerror=0;
#DBG>			};

			    push(@$stack,
                     [ $$states[$$stack[-1][0]]{GOTOS}{$lhs}, $semval ]);
                $$check='';
                next;
            };

#DBG>			$debug & 0x04
#DBG>		and	print STDERR "Forced Error recovery.\n";

            $$check='';

        };

        #Error
            $$errstatus
        or   do {

            $$errstatus = 1;
            &$error($self);
                $$errstatus # if 0, then YYErrok has been called
            or  next;       # so continue parsing

#DBG>			$debug & 0x10
#DBG>		and	do {
#DBG>			print STDERR "**Entering Error recovery.\n";
#DBG>			++$dbgerror;
#DBG>		};

            ++$$nberror;

        };

			$$errstatus == 3	#The next token is not valid: discard it
		and	do {
				$$token eq ''	# End of input: no hope
			and	do {
#DBG>				$debug & 0x10
#DBG>			and	print STDERR "**At eof: aborting.\n";
				return(undef);
			};

#DBG>			$debug & 0x10
#DBG>		and	print STDERR "**Dicard invalid token ".&$ShowCurToken.".\n";

			$$token=$$value=undef;
		};

        $$errstatus=3;

		while(	  @$stack
			  and (		not exists($$states[$$stack[-1][0]]{ACTIONS})
			        or  not exists($$states[$$stack[-1][0]]{ACTIONS}{error})
					or	$$states[$$stack[-1][0]]{ACTIONS}{error} <= 0)) {

#DBG>			$debug & 0x10
#DBG>		and	print STDERR "**Pop state $$stack[-1][0].\n";

			pop(@$stack);
		}

			@$stack
		or	do {

#DBG>			$debug & 0x10
#DBG>		and	print STDERR "**No state left on stack: aborting.\n";

			return(undef);
		};

		#shift the error token

#DBG>			$debug & 0x10
#DBG>		and	print STDERR "**Shift \$error token and go to state ".
#DBG>						 $$states[$$stack[-1][0]]{ACTIONS}{error}.
#DBG>						 ".\n";

		push(@$stack, [ $$states[$$stack[-1][0]]{ACTIONS}{error}, undef ]);

    }

    #never reached
	croak("Error in driver logic. Please, report it as a BUG");

}#_Parse
#DO NOT remove comment

1;

}
#End of include--------------------------------------------------




sub new {
        my($class)=shift;
        ref($class)
    and $class=ref($class);

    my($self)=$class->SUPER::new( yyversion => '1.05',
                                  yystates =>
[
	{#State 0
		ACTIONS => {
			'GROUP' => 1
		},
		GOTOS => {
			'group' => 2,
			'top_level' => 3,
			'group_start' => 4
		}
	},
	{#State 1
		ACTIONS => {
			"=" => 5
		}
	},
	{#State 2
		ACTIONS => {
			'END' => 6
		}
	},
	{#State 3
		ACTIONS => {
			'' => 7
		}
	},
	{#State 4
		ACTIONS => {
			'GROUP' => 1,
			'OBJECT' => 12
		},
		DEFAULT => -8,
		GOTOS => {
			'group_contents' => 8,
			'group' => 10,
			'object' => 9,
			'object_start' => 11,
			'group_start' => 4
		}
	},
	{#State 5
		ACTIONS => {
			'IDENT' => 13
		}
	},
	{#State 6
		DEFAULT => -1
	},
	{#State 7
		DEFAULT => 0
	},
	{#State 8
		ACTIONS => {
			'END_GROUP' => 15
		},
		GOTOS => {
			'group_end' => 14
		}
	},
	{#State 9
		ACTIONS => {
			'GROUP' => 1,
			'OBJECT' => 12
		},
		DEFAULT => -8,
		GOTOS => {
			'group_contents' => 16,
			'group' => 10,
			'object' => 9,
			'object_start' => 11,
			'group_start' => 4
		}
	},
	{#State 10
		ACTIONS => {
			'GROUP' => 1,
			'OBJECT' => 12
		},
		DEFAULT => -8,
		GOTOS => {
			'group_contents' => 17,
			'group' => 10,
			'object' => 9,
			'object_start' => 11,
			'group_start' => 4
		}
	},
	{#State 11
		ACTIONS => {
			'NUM_VAL' => 20,
			'CLASS' => 19
		},
		GOTOS => {
			'num_val' => 18,
			'object_contents' => 21,
			'class' => 22
		}
	},
	{#State 12
		ACTIONS => {
			"=" => 23
		}
	},
	{#State 13
		ACTIONS => {
			'GROUPTYPE' => 24
		},
		DEFAULT => -5,
		GOTOS => {
			'group_type' => 25
		}
	},
	{#State 14
		DEFAULT => -2
	},
	{#State 15
		ACTIONS => {
			"=" => 26
		}
	},
	{#State 16
		DEFAULT => -7
	},
	{#State 17
		DEFAULT => -6
	},
	{#State 18
		ACTIONS => {
			'VALUE' => 27,
			'CLASS' => 19
		},
		GOTOS => {
			'value' => 28,
			'class' => 29
		}
	},
	{#State 19
		ACTIONS => {
			"=" => 30
		}
	},
	{#State 20
		ACTIONS => {
			"=" => 31
		}
	},
	{#State 21
		ACTIONS => {
			'END_OBJECT' => 33
		},
		GOTOS => {
			'object_end' => 32
		}
	},
	{#State 22
		ACTIONS => {
			'GROUP' => 1,
			'NUM_VAL' => 20,
			'OBJECT' => 12
		},
		DEFAULT => -8,
		GOTOS => {
			'num_val' => 35,
			'group_contents' => 34,
			'group' => 10,
			'object' => 9,
			'object_start' => 11,
			'group_start' => 4
		}
	},
	{#State 23
		ACTIONS => {
			'IDENT' => 36
		}
	},
	{#State 24
		ACTIONS => {
			"=" => 37
		}
	},
	{#State 25
		ACTIONS => {
			'CLASS' => 19
		},
		DEFAULT => -23,
		GOTOS => {
			'class_opt' => 38,
			'class' => 39
		}
	},
	{#State 26
		ACTIONS => {
			'IDENT' => 40
		}
	},
	{#State 27
		ACTIONS => {
			"=" => 41
		}
	},
	{#State 28
		ACTIONS => {
			'CLASS' => 19
		},
		DEFAULT => -12,
		GOTOS => {
			'class' => 42
		}
	},
	{#State 29
		ACTIONS => {
			'VALUE' => 27
		},
		GOTOS => {
			'value' => 43
		}
	},
	{#State 30
		ACTIONS => {
			'STRING' => 44
		}
	},
	{#State 31
		ACTIONS => {
			'NUMBER' => 45
		}
	},
	{#State 32
		DEFAULT => -10
	},
	{#State 33
		ACTIONS => {
			"=" => 46
		}
	},
	{#State 34
		DEFAULT => -16
	},
	{#State 35
		ACTIONS => {
			'VALUE' => 27
		},
		GOTOS => {
			'value' => 47
		}
	},
	{#State 36
		DEFAULT => -11
	},
	{#State 37
		ACTIONS => {
			'IDENT' => 48
		}
	},
	{#State 38
		DEFAULT => -3
	},
	{#State 39
		DEFAULT => -22
	},
	{#State 40
		DEFAULT => -9
	},
	{#State 41
		ACTIONS => {
			'PAREN_STRING' => 50,
			'NUMBER' => 51,
			'STRING' => 49
		}
	},
	{#State 42
		DEFAULT => -13
	},
	{#State 43
		DEFAULT => -14
	},
	{#State 44
		DEFAULT => -24
	},
	{#State 45
		DEFAULT => -18
	},
	{#State 46
		ACTIONS => {
			'IDENT' => 52
		}
	},
	{#State 47
		DEFAULT => -15
	},
	{#State 48
		DEFAULT => -4
	},
	{#State 49
		DEFAULT => -19
	},
	{#State 50
		DEFAULT => -20
	},
	{#State 51
		DEFAULT => -21
	},
	{#State 52
		DEFAULT => -17
	}
],
                                  yyrules  =>
[
	[#Rule 0
		 '$start', 2, undef
	],
	[#Rule 1
		 'top_level', 2, undef
	],
	[#Rule 2
		 'group', 3,
sub
#line 25 "MetFile.yp"
{ { "$_[1]" => $_[2] } }
	],
	[#Rule 3
		 'group_start', 5,
sub
#line 29 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 4
		 'group_type', 3, undef
	],
	[#Rule 5
		 'group_type', 0, undef
	],
	[#Rule 6
		 'group_contents', 2,
sub
#line 37 "MetFile.yp"
{ my %a = (%{$_[1]}, %{$_[2]}); \%a }
	],
	[#Rule 7
		 'group_contents', 2,
sub
#line 39 "MetFile.yp"
{ my %a = (%{$_[1]}, %{$_[2]}); \%a }
	],
	[#Rule 8
		 'group_contents', 0,
sub
#line 41 "MetFile.yp"
{ {} }
	],
	[#Rule 9
		 'group_end', 3, undef
	],
	[#Rule 10
		 'object', 3,
sub
#line 48 "MetFile.yp"
{ { "$_[1]" => $_[2] } }
	],
	[#Rule 11
		 'object_start', 3,
sub
#line 52 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 12
		 'object_contents', 2,
sub
#line 56 "MetFile.yp"
{ $_[2] }
	],
	[#Rule 13
		 'object_contents', 3,
sub
#line 58 "MetFile.yp"
{ $_[2] }
	],
	[#Rule 14
		 'object_contents', 3,
sub
#line 60 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 15
		 'object_contents', 3,
sub
#line 62 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 16
		 'object_contents', 2,
sub
#line 64 "MetFile.yp"
{ $_[2] }
	],
	[#Rule 17
		 'object_end', 3, undef
	],
	[#Rule 18
		 'num_val', 3, undef
	],
	[#Rule 19
		 'value', 3,
sub
#line 74 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 20
		 'value', 3,
sub
#line 76 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 21
		 'value', 3,
sub
#line 78 "MetFile.yp"
{ $_[3] }
	],
	[#Rule 22
		 'class_opt', 1, undef
	],
	[#Rule 23
		 'class_opt', 0, undef
	],
	[#Rule 24
		 'class', 3, undef
	]
],
                                  @_);
    bless($self,$class);
}

#line 88 "MetFile.yp"


sub _Error {
    my($parser) = shift;

    $parser->{HAS_ERRORS} = 1;

    if(exists $parser->YYData->{ERRMSG}) {
	$parser->{ERR_PRINT} = $parser->YYData->{ERRMSG};
        delete $parser->YYData->{ERRMSG};
    } else {
	$parser->{ERR_PRINT} = "Syntax error at token: " . $parser->YYCurtok . " with value: " . $parser->YYCurval . "\n";

	my @expect = $parser->YYExpect;
	$parser->{ERR_PRINT} .= "Expected: @expect\n";
	
	my ($parse_location) = ($parser->YYData->{INPUT} =~ /((.|\n){1,80})/);
	$parser->{ERR_PRINT} .= "Parser halted at text:\n$parse_location\n";
    }

    #print $parser->{ERR_PRINT};
}

sub _Lexer {
    my($parser)=shift;

    $parser->YYData->{INPUT}
    or  $parser->YYData->{INPUT} = <STDIN>
    or  return('',undef);

    #print "------>" . $parser->YYData->{INPUT} . "<------\n"; # DEBUGGING

    # Get rid of spaces at the beginning of the current parsing location
    $parser->YYData->{INPUT} =~ s/^\n\s*//;
    $parser->YYData->{INPUT} =~ s/^[ \t]+//;

    for ($parser->YYData->{INPUT}) {

      s/^GROUPTYPE//
	and return('GROUPTYPE', 'GROUPTYPE');

      s/^GROUP//
	and return('GROUP', 'GROUP');

      s/^OBJECT//
	and return('OBJECT', 'OBJECT');

      s/^NUM_VAL//
	and return('NUM_VAL', 'NUM_VAL');

      s/^VALUE//
	and return('VALUE', 'VALUE');

      s/^CLASS//
	and return('CLASS', 'CLASS');

      s/^END_GROUP//
	and return('END_GROUP', 'END_GROUP');

      s/^END_OBJECT//
	and return('END_OBJECT', 'END_OBJECT');

      s/^END//
	and return('END', 'END');

      s/^\=//
	and return('=', '=');

      s/^([-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?)//
	and return('NUMBER',$1);

      s/^\(\"(((\\\")|[^\"])*)\"\)//
	and return('PAREN_STRING',$1);

      s/^\"(((\\\")|[^\"])*)\"//
	and return('STRING',$1);

      s/^([A-Za-z][A-Za-z0-9_]*)//
	and return('IDENT',$1);

    }
}

sub has_error {
  my $self = shift;
  return exists($self->{HAS_ERRORS}) && $self->{HAS_ERRORS};
}

sub get_error_msg {
  my $self = shift;
  return $self->{ERR_PRINT};
}

sub parse {
    my($self, $object_text) = @_;

    $self->YYData->{INPUT} = $object_text;
    my $result = $self->YYParse( yylex => \&_Lexer, yyerror => \&_Error );
}

# Debugging method used only for testing purposes

sub _test {
  use Data::Dumper;
  use strict;

  die "No file" unless -e $ARGV[0];
  open(FILE, $ARGV[0]);
  my @contents = <FILE>;
  close(FILE);

  my $j_contents = join("", @contents);

  my $met_parse = Parse::MetFile->new();
  my $result = $met_parse->parse($j_contents);

  print Dumper($result);
}

#_test();

1;
