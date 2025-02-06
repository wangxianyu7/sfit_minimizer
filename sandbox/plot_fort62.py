import numpy as np
import matplotlib.pyplot as plt

sfit_62 = np.genfromtxt(
    'fort.62', dtype=None,
    names=['nob', 'k', 'time', 'dft0', 'dfu0', 'dftE', 'dfrho'])

for i, key in enumerate(['dft0', 'dfu0', 'dftE', 'dfrho']):
    plt.figure()
    plt.title('{0} {1}'.format(i, key))
    plt.scatter(sfit_62['time'], sfit_62[key])

plt.show()