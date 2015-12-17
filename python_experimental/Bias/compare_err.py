from full_physics import *
import glob
import os
import shelve
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt

pdf = PdfPages("err_plot.pdf")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
for sounding_id in glob.glob("2*"):
    din = shelve.open(sounding_id + "/stat.shlv")
    xco2 = din["xco2"]
    xco2_uncer = din["xco2_uncer"]
    p = din["p"]
    p_uncer = din["p_uncer"]

    plt.figure()
    plt.hist(xco2, 50)
    plt.title("XCO2 Retrieved - True Sounding %s" % sounding_id)
    plt.xlabel("PPM")
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(xco2), np.sqrt(np.cov(xco2)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()

    plt.figure()
    plt.hist(xco2_uncer, 50)
    plt.title("XCO2 Uncertainty Sounding %s" % sounding_id)
    plt.xlabel("PPM")
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(xco2_uncer), np.sqrt(np.cov(xco2_uncer)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()

    t = np.array(xco2) / np.array(xco2_uncer)
    plt.figure()
    plt.hist(t, 50)
    plt.title("(XCO2 Retrieved - True) / Uncertainty\nSounding %s" % sounding_id)
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(t), np.sqrt(np.cov(t)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()


    plt.figure()
    p = np.array(p) / 100
    plt.hist(p, 50)
    plt.title("Surface Pressure Retrieved - True Sounding %s" % sounding_id)
    plt.xlabel("hPa")
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(p), np.sqrt(np.cov(p)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()

    plt.figure()
    p_uncer = np.array(p_uncer) / 100
    plt.hist(p_uncer, 50)
    plt.title("Surface Pressure Uncertainty Sounding %s" % sounding_id)
    plt.xlabel("hPa")
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(p_uncer), np.sqrt(np.cov(p_uncer)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()

    t = np.array(p) / np.array(p_uncer)
    plt.figure()
    plt.hist(t, 50)
    plt.title("(Surface Pressure Retrieved - True) / Uncertainty\nSounding %s" % sounding_id)
    txtstr = '$\mu=%f$\n$\sigma=%f$' % (np.mean(t), np.sqrt(np.cov(t)))
    plt.gca().text(0.05, 0.95, txtstr, transform=plt.gca().transAxes, verticalalignment='top',bbox=props)
    pdf.savefig()

pdf.close()



