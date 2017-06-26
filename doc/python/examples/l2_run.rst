.. _example-l2-run:

Working with Level 2 Runs
============================

A common use case is wanting to work with an existing Level 2 run, to either
get more diagnostic information out than we put in the HDF file, or to attempt
to determine a problem with a run.

While you can work directly with all the C++ classes wrapped in python, 
there is a convenience class to make this a bit easier. Some of the 
interesting classes are a bit nested. This makes sense in the C++ code, but
it can be inconvenient in python.

Setting up for a run, we need to specify the Lua configuration file, 
sounding id, ECMWF file, spectrum file, and if present cloud file.

We have an example here where we set up a forward model run to work
as close to the OCO-2 simulator as possible, and then compare the L2
results with the OCO-2 simulator. This duplicates the test done in
test/oco2_simulator/compare_simulator.py.

.. ipython::

   In [1]: from full_physics import *

   In [1]: import os

   In [1]: import numpy as np
   
   In [1]: from matplotlib.pylab import *

   In [3]: bdir = os.environ["L2_TEST_PATH"] + "/oco2_simulator/"
  
   In [7]: sid_name = "2010090912004075"

   In [7]: ecmwf_name = "../input/oco2_sim_met.h5"

   In [7]: l1b_name = "../input/oco2_sim_l1b.h5"

   In [7]: scene_file = "../input/oco2_sim_scene.h5"

   In [7]: r = L2Run(bdir + "config/config_orbit_sim_match.lua", sid_name,
      ...: ecmwf_name, l1b_name, scene_file = scene_file)

Note that the file names l1b_name etc. are relative to the location of
the Lua config file. Often this can be convenient, since usually the input
is right next to wherever the config file is. But you can also just give
the full path if you find that confusing.

Once we have the L2Run object, we can do simple things like look
at the full configuration of the forward model:

.. ipython::

   In [1]: print r.forward_model

At this point, the state is the initial state set by all the apriori values:

.. ipython::

   In [1]: print r.state_vector

You can look at the current XCO2 value based on this initial state (not the
final value you would have after running a retrieval):

.. ipython::

   In [1]: print r.xco2

We can calculate the high resolution spectrum, and compare to the results
from the OCO-2 simulator. This example does just one band, but you can as 
easily do all 3:

.. ipython::
   
   In [14]: sdir = bdir + "simulator_result/"

   In [14]: hres_ref = np.genfromtxt(sdir + "hires_10774_1.txt", skip_header = 1)

   In [14]: hres_sd = SpectralDomain(hres_ref[:,0], Unit("nanometer"))

   In [14]: refl = hres_ref[:,2]

Here we read in the full text file output, skipping the first line as 
a header. The first column ("0" in zero-based notation) has units added
to create a SpectralDomain. The SpectralDomain is somewhat poor name
for the x-axis of the spectrum - we couldn't find a better name for
"might be wavelength, or might be wavenumber". 

We then calculate the reflectance from the radiative transfer. A couple
of things to note:

* The reflectance function takes a final argument "skip jacobian". We
  set this True here since we don't want the Jacobian. Leave this the
  default value of False if you do want the Jacobian
* The returns a Spectrum object. This has both the wavenumber/wavelength
  with units (a SpectralDomain), and the value with units and possibly
  a Jacobian (a SpectralRange).  This is a nested object in C++, but you
  can just grab the "value", "wavenumber" or "wavelength" as a array by
  using those python helper function.

.. ipython::
   
   In [14]: rcalc = r.radiative_transfer.reflectance(hres_sd, 0, True).value

   In [15]: title("Difference High res spectra band 1");
   
   In [16]: xlabel("Wavelength");

   In [17]: ylabel("Difference %");

   @savefig l2_run_hres_diff.png width=4in
   In [18]: plot(hres_sd.wavelength(), (refl - rcalc) / max(refl) * 100.0);








