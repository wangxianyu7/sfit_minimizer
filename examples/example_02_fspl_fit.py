import sfit_minimizer

import MulensModel as mm
import os.path
import matplotlib.pyplot as plt
from matplotlib import gridspec


mm.utils.MAG_ZEROPOINT = 18.

data_path = os.path.join(sfit_minimizer.DATA_PATH, 'MMTest')
datafiles = ['FSPL_Obs_1_I.pho', 'FSPL_Obs_2_V.pho']
datasets = []
for i, filename in enumerate(datafiles):
    if i == 0:
        band = 'I'
    else:
        band = 'V'

    data = mm.MulensData(
        file_name=os.path.join(data_path, filename), phot_fmt='mag',
        bandpass=band)

    datasets.append(data)

gammas = {'I': 0.44,
          'V': 0.72,
          'H': 0.26}
model = mm.Model({'t_0': 2451697.2, 'u_0': 0.00600, 't_E': 25.00000, 'rho': 0.0065})
n_t_star = 10.
t_star = model.parameters.rho * model.parameters.t_E
model.set_magnification_methods([
    model.parameters.t_0 - n_t_star * t_star,
    'finite_source_LD_Yoo04',
    model.parameters.t_0 + n_t_star * t_star])
for band, value in gammas.items():
    model.set_limb_coeff_gamma(band, value)

event = mm.Event(datasets=datasets, model=model)
print(event)
print(model)
parameters_to_fit = ['t_0', 'u_0', 't_E', 'rho']
initial_guess = [2451697.19995,  0.00600, 25.00000,  0.00650, 1.30000,  0.00000,
                 1.00000,  0.00000]

my_func = sfit_minimizer.mm_funcs.PSPLFunction(event, parameters_to_fit)

result = sfit_minimizer.minimize(
    my_func, x0=initial_guess, tol=1e-4,
    options={'step': 'adaptive'}, verbose=True)

print('Full Results:')
print(result)

values = result.x
sigmas = result.sigmas
print('results: ')
print(values)
print('+/-')
print(sigmas)

my_func.update_all(values)
print('chi2: ', my_func.chi2)

# Plot results
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
plt.figure()
ax11 = plt.subplot(gs[0])
plt.title('Data and Fitted Model (Default)')
my_func.event.plot_model(
    t_range=[(model.parameters.t_0 - 3*t_star), (model.parameters.t_0 + 3*t_star)],
    subtract_2450000=True, bandpass='I', label='I', color='red')
my_func.event.plot_model(
    t_range=[(model.parameters.t_0 - 3*t_star), (model.parameters.t_0 + 3*t_star)],
    subtract_2450000=True, bandpass='V', label='V', color='blue')
my_func.event.plot_data(subtract_2450000=True)
plt.ylim(13.5,11.5)
# Plot the residuals
plt.subplot(gs[1], sharex=ax11)
my_func.event.plot_residuals(subtract_2450000=True)
plt.xlim((model.parameters.t_0 - 2.*t_star - 2450000.),
         (model.parameters.t_0 + 2.*t_star - 2450000.))


plt.show()
