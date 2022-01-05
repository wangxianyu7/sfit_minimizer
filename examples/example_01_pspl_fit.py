import sfit_minimizer
import matplotlib.pyplot as plt
import MulensModel as mm
from matplotlib import gridspec

datasets = [None]  # some data TBD
model = mm.Model({})  # some model TBD
event = mm.Event(datasets=datasets, model=model)

parameters_to_fit = ['t_0', 'u_0', 't_E']
initial_guess = []  # Wrong initial condition
my_func = sfit_minimizer.mm_funcs.PSPLFunction(event, parameters_to_fit)

result = sfit_minimizer.minimize(
    my_func, x0=initial_guess, tol=1e-3, 
    options={'step': 'adaptive'})

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
my_func.event.plot_model(subtract_2450000=True)
my_func.event.plot_data(subtract_2450000=True)
plt.title('Data and Fitted Model (Default)')
# Plot the residuals
plt.subplot(gs[1], sharex=ax11)
my_func.event.plot_residuals(subtract_2450000=True)

plt.show()
