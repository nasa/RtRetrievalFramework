from full_physics import *
import shelve

version = "August 29, 2013"
usage = '''Usage:
  do_retrieval.py [options] <sounding_id> <covariance> <draw> <output>
  do_retrieval.py -h | --help
  do_retrieval.py -v | --version

Script to do a retrieval from a draw, using a supplied x_a and
covariance matrix. We only produce output if the retrieval is 
successful

Options:
  -h --help         
     Print this message

  --xtrue-initial-guess
     If selected, used x_true as the initial guess in the retrieval

  --ptrue-initial-guess
     If selected, use the true surface pressure in the initial guess

  --atrue-initial-guess
     If selected, use the true aerosol in the initial guess

  --noel-initial-guess
     Use an initial guess idea that Noel had (called REG, didn't actually
     work so well).

  --noel-initial-guess2
     Use an initial guess idea that Noel had called RID.

  --lam=f
     Lambda value to use for noel's initial guess 2 (called lam rather than
     lambda because that is a reserved python keyword). [default: 1.0]

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
din = shelve.open(args.covariance, "r")
x_a = din["x_a"]
if("cov_a_original" in din):
    cov_a = din["cov_a_original"]
else:
    cov_a = din["cov_a"]
din = shelve.open(args.draw, "r")
x_true = din["x_true"]
xco2_true = din["xco2_true"]
rad = din["radiance_simulated"]

if(args.xtrue_initial_guess):
    x_initial = np.copy(x_true)
else:
    x_initial = np.copy(x_a)
if(args.ptrue_initial_guess):
    x_initial[21] = x_true[21]
if(args.atrue_initial_guess):
    x_initial[23:35] = x_true[23:35]

basedir = "/groups/algorithm/l2_fp/oco2_simulator_test/automated/tags/B3.5.00_aerosol_testing_0/Orbit022/"

l2run = L2Run.create_from_existing_run(basedir + "oco_Orbit022.config", 
                                       str(args.sounding_id), 0,
                                       "config_quicklook_aerosol_test_0.lua")

# Update what the solver is looking at to match the simulated L1b 
# radiance data
l2run.forward_model.level_1b = Level1bCache(l2run.forward_model.level_1b)
slc = []
sindex = 0
for i in range(3):
    plist = l2run.forward_model.spectral_grid.pixel_list(i)
    l2run.forward_model.level_1b.set_radiance(i, rad[i], plist)
    slc.append(SliceType(sindex, sindex + rad[i].data.shape[0]))
    sindex = sindex + rad[i].data.shape[0]

# Use Noel's idea for initial guess if requested.
if(args.noel_initial_guess):
    residual, se, k = l2run.cost_function.cost_function(x_a)
    z = (residual + np.dot(k,x_a)) / np.sqrt(se)
    d = np.empty(k.shape)
    for i in range(d.shape[1]):
        d[:,i] = (1/np.sqrt(se)) * k[:,i] * np.sqrt(cov_a[i,i])
    beta = np.linalg.lstsq(d, z)[0]
    x_reg = np.empty(beta.shape)
    for i in range(x_reg.shape[0]):
        x_reg[i] = x_a[i] + beta[i] * np.sqrt(cov_a[i,i])
    x_initial = x_reg

# Use Noel's idea for initial guess if requested.
if(args.noel_initial_guess2):
    residual, se, k = l2run.cost_function.cost_function(x_a)
    z = (residual + np.dot(k,x_a)) / np.sqrt(se)
    zcor = np.copy(z)
    for s in slc:
        zcor[s] -= np.mean(z[s])
    k = np.matrix(k)
    d = np.matrix(np.empty(k.shape))
    for i in range(d.shape[1]):
        d[:,i] = (1/np.sqrt(se)) * k[:,i] * np.sqrt(cov_a[i,i])
    lam = args.lam
    beta = np.dot((d.T * d + lam * np.eye(d.shape[1])).I * d.T, zcor)
    x_rid = np.empty(beta.shape[1])
    for i in range(x_rid.shape[0]):
        x_rid[i] = x_a[i] + beta[0,i] * np.sqrt(cov_a[i,i])
    x_initial = x_rid

# Retrieve
have_sol = l2run.solver.solve(x_initial, x_a, cov_a)
if(have_sol):
    l2run.forward_model.state_vector.update_state(l2run.solver.x_solution,
                                                  l2run.solver.aposteriori_covariance)
    res = shelve.open(args.output)
    res["x_initial"] = x_initial
    res["x_a"] = x_a
    res["cov_a"] = cov_a
    res["x_true"] = x_true
    res["xco2_true"] = xco2_true
    res["x_retrieved"] = np.copy(l2run.state_vector.state)
    res["xco2_uncer"] = l2run.error_analysis.xco2_uncertainty * 1e6
    res["xco2_retrieved"] = l2run.xco2
    res["p_uncer"] = np.sqrt(l2run.solver.aposteriori_covariance[21,21])
    res["p_retrieved"] = l2run.atmosphere.pressure.surface_pressure.value.value
    res["solver_state"] = l2run.solver.state
    res.close()
