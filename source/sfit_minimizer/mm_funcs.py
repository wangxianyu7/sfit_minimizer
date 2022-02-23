import numpy as np
import sfit_minimizer


class PSPLFunction(sfit_minimizer.SFitFunction):
    """
    A class for fitting a point-source point-lens microlensing light curve to
    observed data using :py:func:`sfit_minimizer.sfit_minimize.minimize()`.
    Simultaneously fits microlensing parameters and source and blend fluxes for
    each dataset.

    Arguments:
        event: *MulensModel.Event()* object
            event contains datasets and an initial model.

        parameters_to_fit: *list* of *str*
            list of the named model parameters to be fit. (Not including the
            fluxes.)

    """
    # Need to think about what to do to fix fluxes:
    # - how to pass that information to this class?
    # - how to implement fixing the fluxes.

    def __init__(self, event, parameters_to_fit, estimate_fluxes=False):
        self.event = event
        fix_source_flux = {}
        fix_blend_flux = {}
        if estimate_fluxes:
            self.event.fit_fluxes()
            for dataset in self.event.datasets:
                (source_flux, blend_flux) = self.event.get_flux_for_dataset(
                    dataset)
                fix_source_flux[dataset] = source_flux
                fix_blend_flux[dataset] = blend_flux
        else:
            for dataset in self.event.datasets:
                fix_source_flux[dataset] = 1.0
                fix_blend_flux[dataset] = 0.0

        self.event.fix_source_flux = fix_source_flux
        self.event.fix_blend_flux = fix_blend_flux

        self.parameters_to_fit = parameters_to_fit
        flattened_data = self.flatten_data()
        sfit_minimizer.SFitFunction.__init__(self, data=flattened_data)

    def flatten_data(self):
        """ Concatenate good points for all datasets into a single array with
        columns: Date, flux, err.
        """
        for i, dataset in enumerate(self.event.datasets):
            data = [dataset.time[dataset.good], dataset.flux[dataset.good],
                    dataset.err_flux[dataset.good]]
            if i == 0:
                flattened_data = np.array(data)
            else:
                flattened_data = np.hstack((flattened_data, data))

        return flattened_data.transpose()

    def update_all(self, theta0=None, verbose=False):
        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, val, theta0[key])

        for i, dataset in enumerate(self.event.datasets):
            self.event.fix_source_flux[dataset] = theta0[
                len(self.parameters_to_fit) + 2 * i]
            self.event.fix_blend_flux[dataset] = theta0[
                len(self.parameters_to_fit) + 2 * i + 1]

        self.event.fit_fluxes(bad=False)
        if verbose:
            print('new value:', theta0)
            print('fluxes:', self.event.fluxes)

        sfit_minimizer.SFitFunction.update_all(self, theta0, verbose=verbose)

    def calc_residuals(self):
        """Calculate expected values of the residuals"""
        for i, fit in enumerate(self.event.fits):
            res_dataset = fit.get_residuals(phot_fmt='flux', bad=False)
            if i == 0:
                res = np.array(res_dataset[0][fit.dataset.good])
            else:
                res = np.hstack(
                    (res, res_dataset[0][fit.dataset.good]))

        self.residuals = res

    def calc_df(self):
        """
        Calculate the derivatives of the fitting function and store as
        self.df.

        """
        dfunc = None
        for i, fit in enumerate(self.event.fits):
            dA_dparm = fit.get_d_A_d_params_for_point_lens_model(
                self.parameters_to_fit)
            for j, key in enumerate(self.parameters_to_fit):
                x = fit.source_flux * dA_dparm[key][fit.dataset.good]
                if j == 0:
                    dfunc_dataset = x
                else:
                    dfunc_dataset = np.vstack((dfunc_dataset, x))

            # M x N
            dfunc_df_source = np.array(
                [fit.get_data_magnification(bad=False)[fit.dataset.good]])  # 1 x N
            dfunc_df_blend = np.ones((1, np.sum(fit.dataset.good)))  # 1 x N

            # resulting shape should be M+2 x N
            dfunc_dataset = np.vstack(
                (dfunc_dataset, dfunc_df_source, dfunc_df_blend))

            if dfunc is None:
                dfunc = dfunc_dataset
            else:
                x = np.vstack(
                    (dfunc_dataset[:len(self.parameters_to_fit), :],
                        np.zeros( (i * 2, np.sum(fit.dataset.good)) )) )
                y = np.hstack(
                    (np.zeros( (2, dfunc.shape[1]) ), dfunc_dataset[-2:, :],) )
                dfunc = np.hstack((dfunc, x))
                dfunc = np.vstack( (dfunc, y) )

        self.df = dfunc
