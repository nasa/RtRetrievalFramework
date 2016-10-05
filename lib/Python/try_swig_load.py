# Depending on how we are built, we might or might not have SWIG available.
# The python modules need to work either way, so we have this simple test
# and load if found.
try:
    from full_physics_swig import *
    have_full_physics_swig = True
except ImportError:
    have_full_physics_swig = False
