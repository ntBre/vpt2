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
    """load a fort.15 file and return it as a list"""
    ret = []
    with open(filename) as f:
        f.readline()
        for line in f:
            ret.extend([float(x) for x in line.split()])
    return ret

def mass_weight(fcs, ws):
    ret = np.zeros_like(fcs)
    for i, r in enumerate(fcs):
        for j, c in enumerate(r):
            ret[i, j] = fcs[i, j] / math.sqrt(ws[i]*ws[j])
    return ret

atoms = ["O", "H", "H"]
weights = [3*[MASSES[x]] for x in atoms]
# fold
weights = sum(weights, [])
fc2 = load_fort15("crawdad.15")
a = np.reshape(f, (r, r))
fc2 = mass_weight(a, weights)
vals, _ = np.linalg.eigh(fc2)
vals = [0 if abs(x) < 1e-10 else x for x in vals]
conv = 2625.498413 * 1000 / (AVOGADRO * BOHR * BOHR * 1.66054e-27)
for v in vals:
    print(math.sqrt(v*conv)/(100*c*2*math.pi))
# r = int(math.sqrt(len(f)))
# eye = np.diag(np.sqrt(1.0/masses))

# v = eye @ a @ eye
# vals, _ = np.linalg.eigh(v)
# # probably I want to take the eigenvectors corresponding to the
# # non-zero eigenvalues instead of the values themselves
# print(vals[-3:])
