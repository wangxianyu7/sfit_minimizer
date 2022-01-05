import numpy as np
import sfit_minimizer


class PSPLFunction(sfit_minimizer.SFitFunction):
    """
    A class for fitting a point-source point-lens microlensing light curve to observed data using
    :py:func:`sfit_minimizer.sfit_minimize.minimize()`. Simultaneously fits microlensing parameters and source and
    blend fluxes for each dataset.

    Arguments:
        event: *MulensModel.Event()* object
            event contains datasets and an initial model.

        parameters_to_fit: *list* of *str*
            list of the named model parameters to be fit. (Not including the fluxes.)

    """
    # Need to think about what to do to fix fluxes:
    # - how to pass that information to this class?
    # - how to implement fixing the fluxes.

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

    def update_all(self, theta0=None):
        for (key, val) in enumerate(self.parameters_to_fit):
            print(key)
            setattr(self.event.model.parameters, val, theta0[key])

        for i in range(len(self.event.datasets)):
            print(len(self.parameters_to_fit) + 2 * i, len(self.parameters_to_fit) + 2 * i + 1)
            self.event.fit[i].fix_source_flux = theta0[len(self.parameters_to_fit) + 2 * i]
            self.event.fit[i].fix_blend_flux = theta0[len(self.parameters_to_fit) + 2 * i + 1]

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
        dfunc = None
        for i, fit in enumerate(self.event.fits()):
            dfunc_dataset = fit.source_flux * fit.get_df_mulens_per_point(bad=False)
            # M x N
            dfunc_df_source = fit.get_data_magnification(bad=False)  # 1 x N
            dfunc_df_blend = np.zeros((1, np.sum(fit.dataset.good)))  # 1 x N
            print(dfunc_dataset.shape, dfunc_df_source.shape, dfunc_df_blend.shape)

            # resulting shape should be M+2 x N
            dfunc_dataset = np.hstack((dfunc_dataset, dfunc_df_source, dfunc_df_blend))
            print(dfunc_dataset.shape)
            # I think this is wrong, because I need two separate elements of dfunc for each dataset.

            if dfunc is None:
                dfunc = dfunc_dataset
            else:
                dfunc = np.vstack((dfunc, dfunc_dataset))

        self.dfunc = dfunc
        print(self.dfunc.shape)