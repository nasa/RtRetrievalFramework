from full_physics import *
import shelve

version = "September 13, 2013"
usage = '''Usage:
  compare_err_save.py [options] <sounding_id>
  compare_err_save.py -h | --help
  compare_err_save.py -v | --version

Script to do a retrieval from a draw, using a supplied x_a and
covariance matrix. We only produce output if the retrieval is 
successful

Options:
  -h --help         
     Print this message

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
sounding_id = str(args.sounding_id)
basedir = "/groups/algorithm/l2_fp/oco2_simulator_test/automated/tags/B3.5.00_aerosol_testing_0/Orbit022/"
print "Doing sounding", sounding_id
l2run = L2Run.create_from_existing_run(basedir + "oco_Orbit022.config", 
                                       sounding_id, 0,
                                       "config_quicklook_aerosol_test_0.lua")
l2run.forward_model.level_1b = Level1bCache(l2run.forward_model.level_1b)

numdraw = 700
xco2 = []
xco2_uncer = []
p = []
p_uncer = []
for dindex in range(numdraw):
    fname1 = sounding_id + "/retrieval_second_%d.shlv" % dindex
    fname2 = sounding_id + "/draw_%d.shlv" % dindex
    if(os.path.isfile(fname1)):
        din1 = shelve.open(fname1, "r")
        din2 = shelve.open(fname2, "r")
        # Restore state (including solver state) to final solution
        rad = din2["radiance_simulated"]
        xco2_true = din2["xco2_true"]
        p_true = din2["x_true"][21]
        for i in range(3):
            l2run.forward_model.level_1b.set_radiance(i, rad[i])
        l2run.solver.state = din1["solver_state"]
        l2run.forward_model.state_vector.update_state(l2run.solver.x_solution,
                                                      l2run.solver.aposteriori_covariance)
        xco2_uncer.append(l2run.error_analysis.xco2_uncertainty * 1e6)
        xco2.append(l2run.xco2 - xco2_true)
        p.append(l2run.atmosphere.pressure.surface_pressure.value.value - p_true)
        p_uncer.append(np.sqrt(l2run.solver.aposteriori_covariance[21,21]))
        print "Done with %d" % dindex
dout = shelve.open(sounding_id + "/stat.shlv")
dout["xco2"] = xco2
dout["xco2_uncer"] = xco2_uncer
dout["p"] = p
dout["p_uncer"] = p_uncer
