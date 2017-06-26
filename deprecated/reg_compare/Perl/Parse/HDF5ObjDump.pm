####################################################################
#
#    This file was generated using Parse::Yapp version 1.05.
#
#        Don't edit this file, use source file instead.
#
#             ANY CHANGE MADE HERE WILL BE LOST !
#
####################################################################
package Parse::HDF5ObjDump;
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
			"HDF5" => 1
		},
		GOTOS => {
			'file' => 2
		}
	},
	{#State 1
		ACTIONS => {
			'STRING' => 4
		},
		GOTOS => {
			'file_name' => 3
		}
	},
	{#State 2
		ACTIONS => {
			'' => 5
		}
	},
	{#State 3
		ACTIONS => {
			"{" => 6
		}
	},
	{#State 4
		DEFAULT => -2
	},
	{#State 5
		DEFAULT => 0
	},
	{#State 6
		ACTIONS => {
			"DATASET" => 7
		},
		GOTOS => {
			'dataset' => 8
		}
	},
	{#State 7
		ACTIONS => {
			'STRING' => 10
		},
		GOTOS => {
			'dataset_name' => 9
		}
	},
	{#State 8
		ACTIONS => {
			"}" => 11
		}
	},
	{#State 9
		ACTIONS => {
			"{" => 12
		}
	},
	{#State 10
		DEFAULT => -4
	},
	{#State 11
		DEFAULT => -1
	},
	{#State 12
		ACTIONS => {
			"DATATYPE" => 15
		},
		GOTOS => {
			'dataset_info' => 13,
			'dataset_type' => 14
		}
	},
	{#State 13
		ACTIONS => {
			"}" => 16
		}
	},
	{#State 14
		ACTIONS => {
			"DATASPACE" => 18
		},
		GOTOS => {
			'dataset_space' => 17
		}
	},
	{#State 15
		ACTIONS => {
			"H5T_STD_U32LE" => 20,
			"H5T_NATIVE_ULONG" => 19,
			"H5T_STD_I16LE" => 21,
			"H5T_STD_I8LE" => 24,
			"H5T_IEEE_F32LE" => 25,
			"H5T_STD_U16LE" => 26,
			"H5T_STD_I32BE" => 27,
			"H5T_NATIVE_INT" => 28,
			"H5T_NATIVE_ULLONG" => 29,
			"H5T_NATIVE_UINT" => 30,
			"H5T_STD_U64BE" => 31,
			"H5T_STD_U32BE" => 32,
			"H5T_STD_U8LE" => 33,
			"H5T_STD_I64BE" => 34,
			"H5T_STD_U8BE" => 35,
			"H5T_IEEE_F32BE" => 36,
			"H5T_NATIVE_LLONG" => 37,
			"H5T_IEEE_F64BE" => 38,
			"H5T_STD_I8BE" => 39,
			"H5T_COMPOUND" => 40,
			"H5T_STD_I64LE" => 41,
			"H5T_NATIVE_UCHAR" => 42,
			"H5T_STD_U16BE" => 43,
			"H5T_NATIVE_CHAR" => 44,
			"H5T_NATIVE_USHORT" => 45,
			"H5T_NATIVE_DOUBLE" => 47,
			"H5T_NATIVE_SHORT" => 48,
			"H5T_STD_I32LE" => 50,
			"H5T_IEEE_F64LE" => 51,
			"H5T_NATIVE_LONG" => 53,
			"H5T_STRING" => 54,
			"H5T_NATIVE_LDOUBLE" => 55,
			"H5T_STD_U64LE" => 57,
			"H5T_STD_I16BE" => 56,
			"H5T_NATIVE_FLOAT" => 58
		},
		GOTOS => {
			'float' => 46,
			'integer' => 23,
			'atomic_type' => 52,
			'string' => 22,
			'datatype' => 59,
			'compound_type' => 49
		}
	},
	{#State 16
		DEFAULT => -3
	},
	{#State 17
		ACTIONS => {
			"DATA" => 61
		},
		GOTOS => {
			'data' => 60
		}
	},
	{#State 18
		ACTIONS => {
			"SIMPLE" => 64,
			"SCALAR" => 62
		},
		GOTOS => {
			'scalar_space' => 65,
			'simple_space' => 66,
			'dataspace' => 63
		}
	},
	{#State 19
		DEFAULT => -35
	},
	{#State 20
		DEFAULT => -25
	},
	{#State 21
		DEFAULT => -15
	},
	{#State 22
		DEFAULT => -11
	},
	{#State 23
		DEFAULT => -9
	},
	{#State 24
		DEFAULT => -13
	},
	{#State 25
		DEFAULT => -39
	},
	{#State 26
		DEFAULT => -23
	},
	{#State 27
		DEFAULT => -16
	},
	{#State 28
		DEFAULT => -32
	},
	{#State 29
		DEFAULT => -37
	},
	{#State 30
		DEFAULT => -33
	},
	{#State 31
		DEFAULT => -26
	},
	{#State 32
		DEFAULT => -24
	},
	{#State 33
		DEFAULT => -21
	},
	{#State 34
		DEFAULT => -18
	},
	{#State 35
		DEFAULT => -20
	},
	{#State 36
		DEFAULT => -38
	},
	{#State 37
		DEFAULT => -36
	},
	{#State 38
		DEFAULT => -40
	},
	{#State 39
		DEFAULT => -12
	},
	{#State 40
		ACTIONS => {
			"{" => 67
		}
	},
	{#State 41
		DEFAULT => -19
	},
	{#State 42
		DEFAULT => -29
	},
	{#State 43
		DEFAULT => -22
	},
	{#State 44
		DEFAULT => -28
	},
	{#State 45
		DEFAULT => -31
	},
	{#State 46
		DEFAULT => -10
	},
	{#State 47
		DEFAULT => -43
	},
	{#State 48
		DEFAULT => -30
	},
	{#State 49
		DEFAULT => -8
	},
	{#State 50
		DEFAULT => -17
	},
	{#State 51
		DEFAULT => -41
	},
	{#State 52
		DEFAULT => -7
	},
	{#State 53
		DEFAULT => -34
	},
	{#State 54
		ACTIONS => {
			"{" => 68
		}
	},
	{#State 55
		DEFAULT => -44
	},
	{#State 56
		DEFAULT => -14
	},
	{#State 57
		DEFAULT => -27
	},
	{#State 58
		DEFAULT => -42
	},
	{#State 59
		DEFAULT => -6
	},
	{#State 60
		DEFAULT => -5
	},
	{#State 61
		ACTIONS => {
			"{" => 69
		}
	},
	{#State 62
		DEFAULT => -60
	},
	{#State 63
		DEFAULT => -57
	},
	{#State 64
		ACTIONS => {
			"{" => 70
		}
	},
	{#State 65
		DEFAULT => -58
	},
	{#State 66
		DEFAULT => -59
	},
	{#State 67
		ACTIONS => {
			"H5T_STD_U32LE" => 20,
			"H5T_NATIVE_ULONG" => 19,
			"H5T_STD_I16LE" => 21,
			"H5T_STD_I8LE" => 24,
			"H5T_IEEE_F32LE" => 25,
			"H5T_STD_U16LE" => 26,
			"H5T_STD_I32BE" => 27,
			"H5T_NATIVE_INT" => 28,
			"H5T_NATIVE_ULLONG" => 29,
			"H5T_NATIVE_UINT" => 30,
			"H5T_STD_U64BE" => 31,
			"H5T_STD_U32BE" => 32,
			"H5T_STD_U8LE" => 33,
			"H5T_STD_U8BE" => 35,
			"H5T_STD_I64BE" => 34,
			"H5T_IEEE_F32BE" => 36,
			"H5T_NATIVE_LLONG" => 37,
			"H5T_IEEE_F64BE" => 38,
			"H5T_STD_I8BE" => 39,
			"H5T_COMPOUND" => 40,
			"H5T_STD_I64LE" => 41,
			"H5T_NATIVE_UCHAR" => 42,
			"H5T_NATIVE_CHAR" => 44,
			"H5T_STD_U16BE" => 43,
			"H5T_NATIVE_USHORT" => 45,
			"H5T_NATIVE_DOUBLE" => 47,
			"H5T_NATIVE_SHORT" => 48,
			"H5T_STD_I32LE" => 50,
			"H5T_IEEE_F64LE" => 51,
			"H5T_NATIVE_LONG" => 53,
			"H5T_STRING" => 54,
			"H5T_NATIVE_LDOUBLE" => 55,
			"H5T_STD_I16BE" => 56,
			"H5T_STD_U64LE" => 57,
			"H5T_NATIVE_FLOAT" => 58
		},
		DEFAULT => -55,
		GOTOS => {
			'float' => 46,
			'integer' => 23,
			'atomic_type' => 52,
			'string' => 22,
			'member_type_def' => 71,
			'datatype' => 72,
			'compound_type' => 49
		}
	},
	{#State 68
		ACTIONS => {
			"STRSIZE" => 73
		}
	},
	{#State 69
		ACTIONS => {
			"(" => 78,
			'string_data' => 79,
			"{" => 84,
			'REAL' => 83,
			'STRING' => 81
		},
		GOTOS => {
			'complex_space_data' => 85,
			'compound_element' => 82,
			'any_element' => 74,
			'simple_space_data' => 80,
			'any_data_seq' => 77,
			'atomic_element' => 76,
			'complex_element' => 75
		}
	},
	{#State 70
		ACTIONS => {
			"(" => 86
		},
		GOTOS => {
			'dims' => 87
		}
	},
	{#State 71
		ACTIONS => {
			"}" => 88
		}
	},
	{#State 72
		ACTIONS => {
			'STRING' => 90
		},
		GOTOS => {
			'field_name' => 89
		}
	},
	{#State 73
		ACTIONS => {
			'REAL' => 92
		},
		GOTOS => {
			'strsize' => 91
		}
	},
	{#State 74
		ACTIONS => {
			"," => 93
		},
		DEFAULT => -74
	},
	{#State 75
		ACTIONS => {
			"," => 94
		},
		DEFAULT => -71
	},
	{#State 76
		DEFAULT => -76
	},
	{#State 77
		DEFAULT => -70
	},
	{#State 78
		ACTIONS => {
			'REAL' => 83,
			'STRING' => 81
		},
		GOTOS => {
			'atomic_element' => 95
		}
	},
	{#State 79
		ACTIONS => {
			"}" => 96
		}
	},
	{#State 80
		ACTIONS => {
			"}" => 97
		}
	},
	{#State 81
		DEFAULT => -79
	},
	{#State 82
		DEFAULT => -77
	},
	{#State 83
		DEFAULT => -78
	},
	{#State 84
		ACTIONS => {
			"{" => 84,
			'REAL' => 83,
			'STRING' => 81
		},
		GOTOS => {
			'compound_element' => 82,
			'any_element' => 74,
			'any_data_seq' => 98,
			'atomic_element' => 76
		}
	},
	{#State 85
		ACTIONS => {
			"}" => 99
		}
	},
	{#State 86
		ACTIONS => {
			"H5S_UNLIMITED" => 103,
			'REAL' => 102
		},
		GOTOS => {
			'dims_values' => 101,
			'dimension' => 100
		}
	},
	{#State 87
		ACTIONS => {
			"/" => 104
		}
	},
	{#State 88
		DEFAULT => -53
	},
	{#State 89
		ACTIONS => {
			";" => 105
		}
	},
	{#State 90
		DEFAULT => -56
	},
	{#State 91
		ACTIONS => {
			";" => 106
		}
	},
	{#State 92
		DEFAULT => -46
	},
	{#State 93
		ACTIONS => {
			"{" => 84,
			'REAL' => 83,
			'STRING' => 81
		},
		GOTOS => {
			'compound_element' => 82,
			'any_element' => 74,
			'any_data_seq' => 107,
			'atomic_element' => 76
		}
	},
	{#State 94
		ACTIONS => {
			"(" => 78
		},
		GOTOS => {
			'complex_space_data' => 108,
			'complex_element' => 75
		}
	},
	{#State 95
		ACTIONS => {
			")" => 109
		}
	},
	{#State 96
		DEFAULT => -67
	},
	{#State 97
		DEFAULT => -68
	},
	{#State 98
		ACTIONS => {
			"}" => 110
		}
	},
	{#State 99
		DEFAULT => -69
	},
	{#State 100
		ACTIONS => {
			"," => 111
		},
		DEFAULT => -63
	},
	{#State 101
		ACTIONS => {
			")" => 112
		}
	},
	{#State 102
		DEFAULT => -65
	},
	{#State 103
		DEFAULT => -66
	},
	{#State 104
		ACTIONS => {
			"(" => 86
		},
		GOTOS => {
			'dims' => 113
		}
	},
	{#State 105
		ACTIONS => {
			"H5T_STD_U32LE" => 20,
			"H5T_NATIVE_ULONG" => 19,
			"H5T_STD_I16LE" => 21,
			"H5T_STD_I8LE" => 24,
			"H5T_IEEE_F32LE" => 25,
			"H5T_STD_U16LE" => 26,
			"H5T_STD_I32BE" => 27,
			"H5T_NATIVE_INT" => 28,
			"H5T_NATIVE_ULLONG" => 29,
			"H5T_NATIVE_UINT" => 30,
			"H5T_STD_U64BE" => 31,
			"H5T_STD_U32BE" => 32,
			"H5T_STD_U8LE" => 33,
			"H5T_STD_U8BE" => 35,
			"H5T_STD_I64BE" => 34,
			"H5T_IEEE_F32BE" => 36,
			"H5T_NATIVE_LLONG" => 37,
			"H5T_IEEE_F64BE" => 38,
			"H5T_STD_I8BE" => 39,
			"H5T_COMPOUND" => 40,
			"H5T_STD_I64LE" => 41,
			"H5T_NATIVE_UCHAR" => 42,
			"H5T_NATIVE_CHAR" => 44,
			"H5T_STD_U16BE" => 43,
			"H5T_NATIVE_USHORT" => 45,
			"H5T_NATIVE_DOUBLE" => 47,
			"H5T_NATIVE_SHORT" => 48,
			"H5T_STD_I32LE" => 50,
			"H5T_IEEE_F64LE" => 51,
			"H5T_NATIVE_LONG" => 53,
			"H5T_STRING" => 54,
			"H5T_NATIVE_LDOUBLE" => 55,
			"H5T_STD_I16BE" => 56,
			"H5T_STD_U64LE" => 57,
			"H5T_NATIVE_FLOAT" => 58
		},
		DEFAULT => -55,
		GOTOS => {
			'float' => 46,
			'integer' => 23,
			'atomic_type' => 52,
			'string' => 22,
			'member_type_def' => 114,
			'datatype' => 72,
			'compound_type' => 49
		}
	},
	{#State 106
		ACTIONS => {
			"STRPAD" => 115
		}
	},
	{#State 107
		DEFAULT => -75
	},
	{#State 108
		DEFAULT => -72
	},
	{#State 109
		ACTIONS => {
			":" => 116
		}
	},
	{#State 110
		DEFAULT => -80
	},
	{#State 111
		ACTIONS => {
			"H5S_UNLIMITED" => 103,
			'REAL' => 102
		},
		GOTOS => {
			'dims_values' => 117,
			'dimension' => 100
		}
	},
	{#State 112
		DEFAULT => -62
	},
	{#State 113
		ACTIONS => {
			"}" => 118
		}
	},
	{#State 114
		DEFAULT => -54
	},
	{#State 115
		ACTIONS => {
			"H5T_STR_NULLTERM" => 122,
			"H5T_STR_NULLPAD" => 120,
			"H5T_STR_SPACEPAD" => 121
		},
		GOTOS => {
			'strpad' => 119
		}
	},
	{#State 116
		ACTIONS => {
			"{" => 84,
			'REAL' => 83,
			'STRING' => 81
		},
		GOTOS => {
			'compound_element' => 82,
			'any_element' => 123,
			'atomic_element' => 76
		}
	},
	{#State 117
		DEFAULT => -64
	},
	{#State 118
		DEFAULT => -61
	},
	{#State 119
		ACTIONS => {
			";" => 124
		}
	},
	{#State 120
		DEFAULT => -48
	},
	{#State 121
		DEFAULT => -49
	},
	{#State 122
		DEFAULT => -47
	},
	{#State 123
		DEFAULT => -73
	},
	{#State 124
		ACTIONS => {
			"CSET" => 125
		}
	},
	{#State 125
		ACTIONS => {
			"H5T_CSET_ASCII" => 127
		},
		GOTOS => {
			'cset' => 126
		}
	},
	{#State 126
		ACTIONS => {
			";" => 128
		}
	},
	{#State 127
		DEFAULT => -50
	},
	{#State 128
		ACTIONS => {
			"CTYPE" => 129
		}
	},
	{#State 129
		ACTIONS => {
			"H5T_C_S1" => 131,
			"H5T_FORTRAN_S1" => 132
		},
		GOTOS => {
			'ctype' => 130
		}
	},
	{#State 130
		ACTIONS => {
			";" => 133
		}
	},
	{#State 131
		DEFAULT => -51
	},
	{#State 132
		DEFAULT => -52
	},
	{#State 133
		ACTIONS => {
			"}" => 134
		}
	},
	{#State 134
		DEFAULT => -45
	}
],
                                  yyrules  =>
[
	[#Rule 0
		 '$start', 2, undef
	],
	[#Rule 1
		 'file', 5,
sub
#line 25 "HDF5ObjDump.yp"
{ $_[4] }
	],
	[#Rule 2
		 'file_name', 1, undef
	],
	[#Rule 3
		 'dataset', 5,
sub
#line 32 "HDF5ObjDump.yp"
{ { "$_[2]" => $_[4] } }
	],
	[#Rule 4
		 'dataset_name', 1, undef
	],
	[#Rule 5
		 'dataset_info', 3,
sub
#line 39 "HDF5ObjDump.yp"
{ { "TYPE"  => $_[1],
		          "SPACE" => $_[2],
			  "DATA"  => $_[3] } }
	],
	[#Rule 6
		 'dataset_type', 2,
sub
#line 47 "HDF5ObjDump.yp"
{ $_[2] }
	],
	[#Rule 7
		 'datatype', 1, undef
	],
	[#Rule 8
		 'datatype', 1, undef
	],
	[#Rule 9
		 'atomic_type', 1, undef
	],
	[#Rule 10
		 'atomic_type', 1, undef
	],
	[#Rule 11
		 'atomic_type', 1, undef
	],
	[#Rule 12
		 'integer', 1, undef
	],
	[#Rule 13
		 'integer', 1, undef
	],
	[#Rule 14
		 'integer', 1, undef
	],
	[#Rule 15
		 'integer', 1, undef
	],
	[#Rule 16
		 'integer', 1, undef
	],
	[#Rule 17
		 'integer', 1, undef
	],
	[#Rule 18
		 'integer', 1, undef
	],
	[#Rule 19
		 'integer', 1, undef
	],
	[#Rule 20
		 'integer', 1, undef
	],
	[#Rule 21
		 'integer', 1, undef
	],
	[#Rule 22
		 'integer', 1, undef
	],
	[#Rule 23
		 'integer', 1, undef
	],
	[#Rule 24
		 'integer', 1, undef
	],
	[#Rule 25
		 'integer', 1, undef
	],
	[#Rule 26
		 'integer', 1, undef
	],
	[#Rule 27
		 'integer', 1, undef
	],
	[#Rule 28
		 'integer', 1, undef
	],
	[#Rule 29
		 'integer', 1, undef
	],
	[#Rule 30
		 'integer', 1, undef
	],
	[#Rule 31
		 'integer', 1, undef
	],
	[#Rule 32
		 'integer', 1, undef
	],
	[#Rule 33
		 'integer', 1, undef
	],
	[#Rule 34
		 'integer', 1, undef
	],
	[#Rule 35
		 'integer', 1, undef
	],
	[#Rule 36
		 'integer', 1, undef
	],
	[#Rule 37
		 'integer', 1, undef
	],
	[#Rule 38
		 'float', 1, undef
	],
	[#Rule 39
		 'float', 1, undef
	],
	[#Rule 40
		 'float', 1, undef
	],
	[#Rule 41
		 'float', 1, undef
	],
	[#Rule 42
		 'float', 1, undef
	],
	[#Rule 43
		 'float', 1, undef
	],
	[#Rule 44
		 'float', 1, undef
	],
	[#Rule 45
		 'string', 15,
sub
#line 84 "HDF5ObjDump.yp"
{ { "STRING" => { 
		                        "STRSIZE" => $_[4],
					"STRPAD"  => $_[7],
					"CSET"    => $_[10],
					"CTYPE"   => $_[13] 
			              }
                        } 
                      }
	],
	[#Rule 46
		 'strsize', 1, undef
	],
	[#Rule 47
		 'strpad', 1, undef
	],
	[#Rule 48
		 'strpad', 1, undef
	],
	[#Rule 49
		 'strpad', 1, undef
	],
	[#Rule 50
		 'cset', 1, undef
	],
	[#Rule 51
		 'ctype', 1, undef
	],
	[#Rule 52
		 'ctype', 1, undef
	],
	[#Rule 53
		 'compound_type', 4,
sub
#line 112 "HDF5ObjDump.yp"
{ { "$_[1]" => $_[3] } }
	],
	[#Rule 54
		 'member_type_def', 4,
sub
#line 116 "HDF5ObjDump.yp"
{ $_[4]->{$_[2]} = $_[1] ; $_[4] }
	],
	[#Rule 55
		 'member_type_def', 0,
sub
#line 118 "HDF5ObjDump.yp"
{ {} }
	],
	[#Rule 56
		 'field_name', 1, undef
	],
	[#Rule 57
		 'dataset_space', 2,
sub
#line 128 "HDF5ObjDump.yp"
{ $_[2] }
	],
	[#Rule 58
		 'dataspace', 1, undef
	],
	[#Rule 59
		 'dataspace', 1, undef
	],
	[#Rule 60
		 'scalar_space', 1, undef
	],
	[#Rule 61
		 'simple_space', 6,
sub
#line 139 "HDF5ObjDump.yp"
{ { "$_[1]" => { CURRENT_DIMS => $_[3],
				       MAX_DIMS     => $_[5] 
                                     } 
                        }
                      }
	],
	[#Rule 62
		 'dims', 3,
sub
#line 147 "HDF5ObjDump.yp"
{ $_[2] }
	],
	[#Rule 63
		 'dims_values', 1,
sub
#line 151 "HDF5ObjDump.yp"
{ [$_[1]] }
	],
	[#Rule 64
		 'dims_values', 3,
sub
#line 153 "HDF5ObjDump.yp"
{ push(@{$_[3]}, $_[1]); $_[3] }
	],
	[#Rule 65
		 'dimension', 1, undef
	],
	[#Rule 66
		 'dimension', 1, undef
	],
	[#Rule 67
		 'data', 4,
sub
#line 163 "HDF5ObjDump.yp"
{ $_[3] }
	],
	[#Rule 68
		 'data', 4,
sub
#line 165 "HDF5ObjDump.yp"
{ $_[3] }
	],
	[#Rule 69
		 'data', 4,
sub
#line 167 "HDF5ObjDump.yp"
{ $_[3] }
	],
	[#Rule 70
		 'simple_space_data', 1, undef
	],
	[#Rule 71
		 'complex_space_data', 1,
sub
#line 174 "HDF5ObjDump.yp"
{ [$_[1]] }
	],
	[#Rule 72
		 'complex_space_data', 3,
sub
#line 176 "HDF5ObjDump.yp"
{ push(@{$_[3]}, $_[1]); $_[3] }
	],
	[#Rule 73
		 'complex_element', 5,
sub
#line 180 "HDF5ObjDump.yp"
{ { $_[2] => $_[5] } }
	],
	[#Rule 74
		 'any_data_seq', 1,
sub
#line 184 "HDF5ObjDump.yp"
{ [$_[1]] }
	],
	[#Rule 75
		 'any_data_seq', 3,
sub
#line 186 "HDF5ObjDump.yp"
{ push(@{$_[3]}, $_[1]); $_[3] }
	],
	[#Rule 76
		 'any_element', 1, undef
	],
	[#Rule 77
		 'any_element', 1, undef
	],
	[#Rule 78
		 'atomic_element', 1, undef
	],
	[#Rule 79
		 'atomic_element', 1, undef
	],
	[#Rule 80
		 'compound_element', 3,
sub
#line 198 "HDF5ObjDump.yp"
{ $_[2] }
	]
],
                                  @_);
    bless($self,$class);
}

#line 202 "HDF5ObjDump.yp"


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

    ##print $parser->{ERR_PRINT}; # DEBUGGING
}

sub _Lexer {
    my($parser)=shift;

    $parser->YYData->{INPUT}
    or  $parser->YYData->{INPUT} = <STDIN>
    or  return('',undef);

    #print "------>" . $parser->YYData->{INPUT} . "<------\n"; # DEBUGGING

    # Get rid of spaces at the beginning of the current parsing location
    $parser->YYData->{INPUT} =~ s/^\n//;
    $parser->YYData->{INPUT} =~ s/^[ \t]+//;

    for ($parser->YYData->{INPUT}) {

	s/^{//
	    and return('{', '{');
	s/^}//
	    and return('}', '}');
	s/^;//
	    and return(';', ';');

        s/^([-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?)//
	    and return('REAL',$1);
	s/^\"(((\\\")|[^\"])*)\"//
	    and return('STRING',$1);

        s/^([A-Za-z][A-Za-z0-9_]*)//
	    and return($1,$1);

        s/^(.)//s
	    and return($1,$1);
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

  die "No file" unless -e $ARGV[0];
  open(FILE, $ARGV[0]);
  my @contents = <FILE>;
  close(FILE);

  my $j_contents = join("", @contents);

  my $hdf5_parse = Parse::HDF5ObjDump->new();
  my $result = $hdf5_parse->parse($j_contents);

  print Dumper($result);
  my @a = keys(%{$result});
}

#_test();

1;
