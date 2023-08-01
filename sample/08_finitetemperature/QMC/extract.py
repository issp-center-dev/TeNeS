import numpy as np
import pyalps
import pyalps.plot as plot

#         Observable name : output_prefix
names = {
    "Energy Density": "ene",
    "Specific Heat": "spec",
}

xnames = ["T"]
foreachs = [["L", "Gamma"]]
fe_types = [[np.int64, np.float64]]


def extract(data, xname, names, foreach, fe_types):
    if np.isscalar(foreach):
        foreach = [foreach]
    if np.isscalar(fe_types):
        fe_types = [fe_types]
    for name in names:
        for obs in pyalps.collectXY(data, xname, name, foreach=foreach):
            vals = [typ(obs.props[sym]) for sym, typ in zip(foreach, fe_types)]
            filename = names[name]
            for sym, val in zip(foreach, vals):
                filename += "-{}{}".format(sym, val)
            filename += ".dat"
            with open(filename, "w") as f:
                f.write(plot.convertToText([obs]).replace(" +/- ", " "))


result_files = pyalps.getResultFiles(prefix="params")

data = pyalps.loadMeasurements(result_files, names.keys())

for xname, fe, fet in zip(xnames, foreachs, fe_types):
    extract(data, xname, names, fe, fet)
