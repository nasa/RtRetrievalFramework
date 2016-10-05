# This aggregates all the results together
from full_physics import *
import cPickle
import scipy.io

x_a = []
x_true = []
xco2_true = []
x_ret = []
xco2_ret = []
flist = glob.glob("job*.pkl")
flist.sort()
for fname in flist:
    with open(fname) as f:
        t = cPickle.load(f)
    xco2_a = t["XCO2_a"]
    x_a.append(t["X_a"])
    x_true.append(t["X_true"])
    xco2_true.append(t["XCO2_true"])
    x_ret.append(t["X_retrieved"])
    xco2_ret.append(t["XCO2_retrieved"])

with open("my_cov.pkl") as f:
    dist_cov = cPickle.load(f)

res = {}
res["Distribution_covariance"] =  dist_cov
res["X_a"] = np.vstack(x_a)
res["XCO2_a"] = xco2_a
res["X_true"] = np.vstack(x_true)
res["XCO2_true"] = np.vstack(xco2_true)
res["X_retrieved"] = np.vstack(x_ret)
res["XCO2_retrieved"] = np.vstack(xco2_ret)
scipy.io.savemat("l2_simulate.mat", res)

