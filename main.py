import numpy as np


def load_fort15():
    """load a fort.15 file and return it as a list"""
    ret = []
    with open("fort.15") as f:
        f.readline()
        for line in f:
            ret.extend([float(x) for x in line.split()])
    return ret


f = load_fort15()
a = np.reshape(f, (9, 9))
vals, _ = np.linalg.eigh(a)
print(vals[-3:])
