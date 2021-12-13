import sfit_minimize
import matplotlib.pyplot as plt


class PSPLFunction(sfit_minimize.SFitFunction):

    def __init__(self, event):
        self.event = event
        flattened_data = self.flatten_data()
        print(np.array(flattened_data).shape)
        sfit_minimize.SFitFunction.__init__(data=flattened_data, theta=theta)

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


initial_guess = [] # Wrong initial condition
my_func = PSPLFunction(event)

result = sfit_minimize.minimize(
    my_func, x0=initial_guess, args=(parameters_to_fit), tol=1e-3, 
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
