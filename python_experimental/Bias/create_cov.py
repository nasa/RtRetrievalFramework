from full_physics import *
import shelve
import glob
import re
import os

# Short script to create a covariance matrix to use in the initial step.
# This is currently hard coded for a specific set of data, since this is
# just a prototype
basedir = "/groups/algorithm/l2_fp/oco2_simulator_test/automated/tags/B3.5.00_aerosol_testing_0/Orbit022/"

for f in glob.glob(basedir + "output/l2_*.h5"):
    m = re.search("l2_(.*).h5", os.path.basename(f))
    sid = m.group(1)
    l2run = L2Run.create_apriori_state(basedir + "oco_Orbit022.config", sid, 0, "config_quicklook_aerosol_test_0.lua")
    # Gives us x_a and the apriori covariance matrix
    x_a = np.copy(l2run.state_vector.state)
    cov_a = np.copy(l2run.state_vector.state_covariance)
    cov_a_original = np.copy(l2run.state_vector.state_covariance)
    cov_p = l2run.output_file.read_double_3d("/RetrievalResults/aposteriori_covariance_matrix")[0,:,:]
    # For now, only work with land data. Determine this by looking at length
    # of x_a
    if(len(x_a) != 46):
        continue
    print sid
    # CO2 VMR ok
    # H2O scaling uncertainty change to 0.1
    cov_a[20, 20] = 0.1 ** 2
    # Surface pressure ok
    # Temperature uncertainty 2
    cov_a[22, 22] = 2.0 ** 2
    # Aerosol ok
    # Albedo is 2 * aposteriori. For simplicity, we only use the diagonal 
    # elements
    for i in range(35,41):
        cov_a[i, i] = cov_p[i, i] * 2

    # Hold instrument dispersion and fluorescence constant by giving really
    # small covariance
    cov_a[41,41] = 1e-8**2
    cov_a[42,42] = 1e-8**2
    cov_a[43,43] = 1e-8**2
    cov_a[44,44] = 1e-8**2
    cov_a[45,45] = 1e-8**2

    os.mkdir(sid)
    d = shelve.open(sid + "/cov_initial.shlv")
    d["x_a"] = x_a
    d["cov_a"] = cov_a
    d["cov_a_original"] = cov_a_original
    d.close()
