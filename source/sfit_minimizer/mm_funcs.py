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
    # Need to think about what to do to fix fluxes:
    # - how to pass that information to this class?
    # - how to implement fixing the fluxes.

    def __init__(self, event, parameters_to_fit, estimate_fluxes=False):

        def initialize_fluxes(self):
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

        self.event = event
        self.parameters_to_fit = parameters_to_fit
        self._set_flux_indices()
        initialize_fluxes(self)

        # Store information about fixed fluxes
        # self.fix_source_flux = {}
        # self.fix_blend_flux = {}
        # for dataset in self.event.datasets:
        #     if dataset in self.event.fix_source_flux.keys():
        #         self.fix_source_flux[dataset] = True
        #     else:
        #         self.fix_source_flux[dataset] = False
        #
        #     if dataset in self.event.fix_blend_flux.keys():
        #         self.fix_blend_flux[dataset] = True
        #     else:
        #         self.fix_blend_flux[dataset] = False

        # Figure out flux indices



        # print('fix_flux parameters')
        # print('event', self.event.fix_source_flux, self.event.fix_blend_flux)
        # event.fit_fluxes()
        # print('fluxes', event.fluxes)
        # print('check fits')
        # for fit in event.fits:
        #     print(
        #         '(True, True)', fit.fix_source_flux is not False,
        #         fit.fix_blend_flux is not False)
        #print('func', self.fix_source_flux, self.fix_blend_flux)

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
                    fs_index = len(self.parameters_to_fit)
                else:
                    fs_index = ind_num + n

            if dataset in self.event.fix_source_flux.keys():
                if self.event.fix_source_flux[dataset] is not False:
                    fs_index = None

            self.fs_indices.append(fs_index)
            # print('fs_indices')
            # print(self.event.fix_source_flux)
            # print(self.fs_indices)

            if len(self.fb_indices) == 0:
                fb_index = len(self.parameters_to_fit)
                if self.fs_indices[-1] is not None:
                    fb_index += 1

            else:
                if self.fs_indices[-1] is not None:
                    n = 2
                else:
                    n = 1

                # #condition = (np.array(self.fb_indices) is None)
                # #print(condition)
                # #ind_num = np.where(condition)
                # print('ind_num')
                # print(self.fb_indices)
                # print(ind_num)
                # print(len(ind_num[0]))
                # if len(ind_num[0]) == 0:
                #     fb_index = len(self.parameters_to_fit)
                #     if self.fs_indices[-1] is not None:
                #         fb_index += 1
                # else:
                #     fb_index = np.nanmax(self.fb_indices) + n

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
            # print(self.event.fix_blend_flux)
            # print(self.fb_indices)
            # print('fb_indices')

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
        # This leads to all kinds of problems if theta0 is not defined.

        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, val, theta0[key])

        for i, dataset in enumerate(self.event.datasets):
            #if not self.fix_source_flux[dataset]:
            #    self.event.fix_source_flux[dataset] = theta0[
            #        len(self.parameters_to_fit) + 2 * i]
            if self.fs_indices[i] is not None:
                self.event.fix_source_flux[dataset] = theta0[self.fs_indices[i]]

            #if not self.fix_blend_flux[dataset]:
            #    self.event.fix_blend_flux[dataset] = theta0[
            #        len(self.parameters_to_fit) + 2 * i + 1]

            if self.fb_indices[i] is not None:
                self.event.fix_blend_flux[dataset] = theta0[self.fb_indices[i]]

        # print('check fits 2')
        # for fit in self.event.fits:
        #     print(
        #         '(True, True)', fit.fix_source_flux is not False,
        #         fit.fix_blend_flux is not False)

        self.event.fit_fluxes(bad=False)
        #self.event.get_chi2()
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
                # pad flux columns with zeros for data points from other datasets
                y = np.hstack(
                    (np.zeros( (n_fluxes_dataset, dfunc.shape[1]) ),
                     dfunc_dataset[-n_fluxes_dataset:, :],) )
                dfunc = np.hstack((dfunc, x)) # goes below existing
                dfunc = np.vstack( (dfunc, y) ) # goes next to existing

            n_fluxes += n_fluxes_dataset

        print('df checks')
        print('theta0', len(self.theta))
        print('dfunc', dfunc.shape)
        self.df = dfunc

    # def calc_bmat(self):
    #     sfit_minimizer.SFitFunction.calc_bmat(self)
    #     for i, dataset in enumerate(self.event.datasets):
    #         if self.fix_source_flux[dataset] is not False:
    #             self.bmat[len(self.parameters_to_fit) + 2 * i,
    #                       len(self.parameters_to_fit) + 2 * i] = 10e15
    #
    #         if self.fix_blend_flux[dataset] is not False:
    #             self.bmat[len(self.parameters_to_fit) + 2 * i + 1,
    #                       len(self.parameters_to_fit) + 2 * i + 1] = 10e15
