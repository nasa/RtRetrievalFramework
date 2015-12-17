.. _example-forward-model:

Calculating a Forward Model (A bit dated)
=========================================

**This example is a bit dated, but we will leave  this in place. You should
look at the example** :ref:`example-l2-run`

Here is an example that uses a basic forward model calculation.

First we import the modules we will be using.

.. ipython::

   In [1]: from full_physics import *

   In [3]: import os

   In [4]: import imp

Then, we set up to use one of the TCCON sounding we use in our unit testing.
You can use other data if you like, this just happens to be convenient data
that we have around.

Note that we communicate the files and sounding id through environment 
variables. The standard configuration file for production does this because
this is a nice interface when writing scripts. Since we are using this
standard config file rather than writing out own, you need to use the
same interface.

An alternative is to write your own configuration file that has these
values hard wired. Use whichever way is the most convenient.

.. ipython::

   In [6]: tccon_path = os.environ["L2_TCCON_SMALL_SET_PATH"]

   In [7]: os.environ["spectrum_file"] = tccon_path + "/acos_L1bB2900_tccon_5_good_qual.h5"

   In [8]: os.environ["ecmwf_file"] = tccon_path + "/acos_EcmB2900_tccon_5_good_qual.h5"

   In [9]: os.environ["cloud_file"] = tccon_path + "/acos_CldB2900_tccon_5_good_qual.h5"

   In [10]: os.environ["sounding_id"] = "20090827005603"

Now we load in the configuration file. This creates all the objects.

.. ipython::

   In [17]: config = imp.load_source('config', os.environ["L2_INPUT_PATH"] + "/gosat/config/config.py")

Run the forward model for the first band:
   
.. ipython::

   In [24]: spec_0 = config.lua_config.forward_model.radiance(0)

You can get the wavenumbers like so:

.. ipython::

   In [25]: spec_0.wavenumber

The high resolution radiance:

.. ipython::

   In [26]: spec_0.value

The high resolution jacobian:

.. ipython::

    In [27]: spec_0.spectral_range.data_ad.jacobian
