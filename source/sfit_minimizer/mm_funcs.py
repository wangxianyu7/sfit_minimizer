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

    Note: if you want to fix the source or blend flux for a particular dataset,
    use the *fix_source_flux* or *fix_blend_flux* keywords in *event* as usual.

    """

    def __init__(self, event, parameters_to_fit, estimate_fluxes=False):
        self.event = event
        self.parameters_to_fit = parameters_to_fit
        self._set_flux_indices()
        self._initialize_fluxes(estimate_fluxes)
        flattened_data = self.flatten_data()
        sfit_minimizer.SFitFunction.__init__(self, data=flattened_data)

    def _set_flux_indices(self):
        # must be carried out before initialize_fluxes()
        self.fs_indices = []
        self.fb_indices = []
        for dataset in self.event.datasets:
            if len(self.fs_indices) == 0:
                fs_index = len(self.parameters_to_fit)
            else:
                if self.fb_indices[-1] is not None:
                    n = 2
                else:
                    n = 1

                ind_num = -1
                for j in self.fs_indices:
                    if j is not None:
                        if j > ind_num:
                            ind_num = j

                if ind_num < 0:
                    if self.fb_indices[-1] is not None:
                        fs_index = self.fb_indices[-1] + 1
                    else:
                        fs_index = len(self.parameters_to_fit)

                else:
                    fs_index = ind_num + n

            if dataset in self.event.fix_source_flux.keys():
                if self.event.fix_source_flux[dataset] is not False:
                    fs_index = None

            self.fs_indices.append(fs_index)
            if len(self.fb_indices) == 0:
                fb_index = len(self.parameters_to_fit)
                if self.fs_indices[-1] is not None:
                    fb_index += 1

            else:
                if self.fs_indices[-1] is not None:
                    n = 2
                else:
                    n = 1

                ind_num = -1
                for j in self.fb_indices:
                    if j is not None:
                        if j > ind_num:
                            ind_num = j

                if ind_num < 0:
                    if self.fs_indices[-1] is not None:
                        fb_index = self.fs_indices[-1] + 1
                    else:
                        fb_index = len(self.parameters_to_fit)

                else:
                    fb_index = ind_num + n

            if dataset in self.event.fix_blend_flux.keys():
                if self.event.fix_blend_flux[dataset] is not False:
                    fb_index = None

            self.fb_indices.append(fb_index)

    def _initialize_fluxes(self, estimate_fluxes):
        # Set initial source fluxes
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
                if dataset in self.event.fix_source_flux.keys():
                    fix_source_flux[dataset] = self.event.fix_source_flux[
                        dataset]
                else:
                    fix_source_flux[dataset] = 1.0

                if dataset in self.event.fix_blend_flux.keys():
                    fix_blend_flux[dataset] = self.event.fix_blend_flux[
                        dataset]
                else:
                    fix_blend_flux[dataset] = 0.0

        self.event.fix_source_flux = fix_source_flux
        self.event.fix_blend_flux = fix_blend_flux

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
        if theta0 is None:
            raise ValueError('theta0 must be passed to update_all()')

        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, val, theta0[key])

        for i, dataset in enumerate(self.event.datasets):
            if self.fs_indices[i] is not None:
                self.event.fix_source_flux[dataset] = theta0[self.fs_indices[i]]

            if self.fb_indices[i] is not None:
                self.event.fix_blend_flux[dataset] = theta0[self.fb_indices[i]]

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
        n_fluxes = 0
        for i, fit in enumerate(self.event.fits):
            n_fluxes_dataset = 0
            dA_dparm = fit.get_d_A_d_params_for_point_lens_model(
                self.parameters_to_fit)
            for j, key in enumerate(self.parameters_to_fit):
                x = fit.source_flux * dA_dparm[key][fit.dataset.good]
                if j == 0:
                    dfunc_dataset = x
                else:
                    dfunc_dataset = np.vstack((dfunc_dataset, x))

            # 1 x N
            if self.fs_indices[i] is not None:
                dfunc_df_source = np.array(
                    [fit.get_data_magnification(bad=False)[fit.dataset.good]])
                dfunc_dataset = np.vstack(
                    (dfunc_dataset, dfunc_df_source))
                n_fluxes_dataset += 1

            if self.fb_indices[i] is not None:
                dfunc_df_blend = np.ones((1, np.sum(fit.dataset.good)))
                dfunc_dataset = np.vstack(
                    (dfunc_dataset, dfunc_df_blend))
                n_fluxes_dataset += 1

            # resulting shape should be M+2 x N
            if dfunc is None:
                dfunc = dfunc_dataset
            else:
                # add columns of zeros for flux parameters for other datasets.
                x = np.vstack(
                    (dfunc_dataset[:len(self.parameters_to_fit), :],
                        np.zeros( (n_fluxes, np.sum(fit.dataset.good)) )) )
                # pad flux columns with zeros for datapoints from other datasets
                y = np.hstack(
                    (np.zeros( (n_fluxes_dataset, dfunc.shape[1]) ),
                     dfunc_dataset[-n_fluxes_dataset:, :],) )
                dfunc = np.hstack((dfunc, x))  # goes below existing
                dfunc = np.vstack( (dfunc, y) )  # goes next to existing

            n_fluxes += n_fluxes_dataset

        self.df = dfunc
