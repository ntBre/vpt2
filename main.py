import numpy as np


def load_fort15():
    """load a fort.15 file and return it as a list"""
    ret = []
    with open("fort.15") as f:
        f.readline()
        for line in f:
            ret.extend([float(x) for x in line.split()])
    return ret


H = 1.007825032
O = 15.994914619

masses = np.array([H, H, H, O, O, O, H, H, H])

f = load_fort15()
a = np.reshape(f, (9, 9))
eye = np.diag(np.sqrt(1.0/masses))

v = eye @ a @ eye
vals, _ = np.linalg.eigh(v)
# probably I want to take the eigenvectors corresponding to the
# non-zero eigenvalues instead of the values themselves
print(vals[-3:])
