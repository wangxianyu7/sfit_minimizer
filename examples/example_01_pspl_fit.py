import sfit_minimizer
import matplotlib.pyplot as plt
import numpy as np
import MulensModel as mm
from matplotlib import gridspec

class PSPLFunction(sfit_minimizer.SFitFunction):

    def __init__(self, event, parameters_to_fit):
        self.event = event
        self.parameters_to_fit = parameters_to_fit
        flattened_data = self.flatten_data()
        print(np.array(flattened_data).shape)
        sfit_minimizer.SFitFunction.__init__(self, data=flattened_data)

    def flatten_data(self):
        """ Concatenate good points for all datasets into a single array with 
        columns: Date, flux, err.
        """
        flattened_data = []
        for dataset in self.event.datasets():
            data = [dataset.time[dataset.good], dataset.flux[dataset.good], 
                    dataset.err_flux[dataset.good]]
            flattened_data.append(data)

        print(np.array(flattened_data).shape)

        return np.flatten(flattened_data, axis=0)

    def update_all(self, theta0):
        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, val, theta0[key])

        sfit_minimizer.SFitFunction.update_all(self, theta0)

    def calc_res(self):
        """Calculate expected values of the residuals"""
        # Need to add a thing to update the model parameters using theta.
        
        res = []
        for fit in self.event.fits():
            res_dataset = fit.get_residuals(phot_fmt='flux', bad=False)
            res.append(res_dataset[0])

        self.res = np.flatten(res, axis=0)
        print(self.res.shape)

    def calc_df(self):
        """
        Calculate the derivatives of the fitting function and store as
        self.df.

        """
        #df = self.event.df_per_point()
        df = []
        for fit in self.event.fits():
            df_dataset = fit.get_df_mulens_per_point(bad=False)
            # Need to add flux factors
            df.append(df_dataset)  # Probably need to replace with a numpy concatenation/stack function.

        # Need to also add partials wrt flux
        self.df = np.flatten(np.array(df), axis=0)
        print(self.df.shape)


datasets = [None]  # some data TBD
model = mm.Model({})  # some model TBD
event = mm.Event(datasets=datasets, model=model)

parameters_to_fit = ['t_0', 'u_0', 't_E']
initial_guess = [] # Wrong initial condition
my_func = PSPLFunction(event, parameters_to_fit)

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
