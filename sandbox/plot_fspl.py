import numpy as np
import matplotlib.pyplot as plt
import MulensModel as mm

data = mm.MulensData(
    file_name='/Users/jyee/PycharmProjects/sfit_minimizer/data/MMTest/FSPL_Obs_1_I.pho',
    phot_fmt='mag')

t_0 = 2451697.19995
t_E = 25.
rho = 0.0065
t_star = t_E * rho
n = 100
print(t_0-n*t_star, t_0+n*t_star)

data.plot()
plt.axvline(t_0-n*t_star)
plt.axvline(t_0+n*t_star)
plt.show()