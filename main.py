import math
import numpy as np

MASSES = {
    "H":  1.007825032,
    "O": 15.994914619,
    }
AVOGADRO = 6.02214076e+23
c        = 2.99792458e+08
BOHR     = 5.29177210903e-11

def load_fort15(filename):
    """load a fort.15 file and return it as a matrix"""
    ret = []
    with open(filename) as f:
        f.readline()
        for line in f:
            ret.extend([float(x) for x in line.split()])
    r = int(math.sqrt(len(ret)))
    return np.reshape(ret, (r, r))

def mass_weight(fcs, ws):
    """Mass-weight the force constant matrix using F^M = m^-1/2 F
    m^-1/2"""
    d = np.diag(1/np.sqrt(ws))
    return d @ fcs @ d

crawford = False
if crawford:
    atoms = ["O", "H", "H"]
    fc2 = load_fort15("crawdad.15")
else:
    atoms = ["H", "O", "H"]
    fc2 = load_fort15("h2o.spectro.15")

weights = [3*[MASSES[x]] for x in atoms]
weights = sum(weights, [])
fc2 = mass_weight(fc2, weights)
vals, _ = np.linalg.eigh(fc2)
vals = [0 if abs(x) < 1e-10 else x for x in vals]

# fcs in Eh/ao²
if crawford:
    conv = 2625.498413 * 1000 / (AVOGADRO * BOHR * BOHR * 1.66054e-27)
    for v in vals:
        print("%.4f" % (math.sqrt(v*conv)/(100*c*2*math.pi)))
# fcs in attojoules
else:
    # 5.034e+3 from
    # http://laser.chem.olemiss.edu/~nhammer/constants.html
    # 1.02 is a magic number to match spectro output
    a = np.sqrt(vals)[-3:]
    h = 6.62607015e-34
    print(a*1e-21/(h*c))
    print(a[2]/a[0])
    # vals = np.sqrt(vals) * 5.034e+3 * 1.02115
    # print(vals)

