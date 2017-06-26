from nose.tools import *
from full_physics.docopt_simple_m import *

def test_docopt_simple():
    '''Test of docopt_simple'''
    usage= '''
Usage: my_test [options] <in1> <in2> <out>
       my_test -h | --help
       my_test --version

This is a simple test.

Options:
  -h --help
     Print out help
  --my-int=n
     A integer argument [default: 1]
  --my_underscore
     A argument with a underscore
  --my-string=s
     A string argument [default: ]
  --my-float=f
     A float argument
  -v --version
     Print program version
'''
    d = docopt_simple(usage, argv=["a1", "a2", "o3"])
    assert "in1" in d
    assert "foo" not in d
    assert d.in1 == "a1"
    assert d.in2 == "a2"
    assert d.out == "o3"
    assert d.my_underscore == False
    d2 = docopt_simple(usage, argv=["--my_underscore", "a1", "a2", "o3"])
    assert d2.my_underscore == True
    assert d.my_string == ""
    d2 = docopt_simple(usage, argv=["--my-string=foo", "a1", "a2", "o3"])
    assert d2.my_string == "foo"
    assert d.my_int == 1
    d2 = docopt_simple(usage, argv=["--my-int=+10", "a1", "a2", "o3"])
    assert d2.my_int == 10
    d2 = docopt_simple(usage, argv=["--my-int=-21", "a1", "a2", "o3"])
    assert d2.my_int == -21
    d2 = docopt_simple(usage, argv=["--my-float=1.2", "a1", "a2", "o3"])
    assert d2.my_float == 1.2
    d2 = docopt_simple(usage, argv=["--my-float=+1.2", "a1", "a2", "o3"])
    assert d2.my_float == 1.2
    d2 = docopt_simple(usage, argv=["--my-float=-1.2", "a1", "a2", "o3"])
    assert d2.my_float == -1.2
    d2 = docopt_simple(usage, argv=["--my-float=-1.", "a1", "a2", "o3"])
    assert d2.my_float == -1.
    d2 = docopt_simple(usage, argv=["--my-float=.3", "a1", "a2", "o3"])
    assert d2.my_float == .3
    d2 = docopt_simple(usage, argv=["--my-float=1e3", "a1", "a2", "o3"])
    assert d2.my_float == 1e3

    


