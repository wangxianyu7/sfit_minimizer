import sfit_minimizer
import matplotlib.pyplot as plt


class PSPLFunction(sfit_minimizer.SFitFunction):

    def __init__(self, event, parameters_to_fit):
        self.event = event
        self.parameters_to_fit = parameters_to_fit
        flattened_data = self.flatten_data()
        print(np.array(flattened_data).shape)
        sfit_minimizer.SFitFunction.__init__(data=flattened_data, theta=theta)

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
            setattr(event.model.parameters, val, theta0[key])

        sfit_minimizer.SFitFunction.update_all(theta0)

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
        """Calculate the derivatives of the fitting function and store as 
        self.df."""
        df = self.event.df_per_point()
        for fit in self.event.fits():
            df_dataset = fit.get_df_per_point(bad=False)
            df.append(df_dataset)

        self.df = np.flatten(np.array(df), axis=0)
        print(self.df.shape)


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

my_func.theta = values
my_func.update_all()
print('chi2: ', my_func.chi2)

plt.figure()
plt.title('Values')
plt.errorbar(
    my_func.data[:, 0], my_func.data[:, 1], 
    yerr=my_func.data[:, 2], fmt='o')
x = np.arange(0, 100)
plt.plot(x, theta_new[0] + theta_new[1] * x, color='red', zorder=5)

plt.figure()
plt.title('Residuals')
plt.errorbar(
    my_func.data[:, 0], my_func.res, 
    yerr=my_func.data[:, 2], fmt='o')
plt.axhline(0, color='red', zorder=5)

plt.show()
