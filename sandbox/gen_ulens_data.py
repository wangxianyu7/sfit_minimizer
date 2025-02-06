"""
Generate some datasets from a microlensing model.
"""
import MulensModel as mm
import os.path
import numpy as np
import matplotlib.pyplot as plt
import sfit_minimizer


def generate_fake_data(
        model=None, source_flux=None, blend_flux=None, band=None,
        t_range=None,
        times=None, n_pts=100, randomize=False,
        fractional_uncertainty=0.01, output_file=None, plot=False):
    """
    Generate a fake dataset from a microlensing model.

    Keyword Arguments:
        model: *MulensModel.Model*, required
            The microlensing model from which to generate the data.

        source_flux: *float*, required
            flux of the source.

        blend_flux: *float*, required
            flux of the blend.

        band: *str*, optional
            bandpass of the dataset, e.g. "I".

        t_range: 2-element *list* or *tuple*, optional
            time range over which to generate the fake data.

        times: *np.array*
            list of times for generating the data.

        n_pts: *int*, optional
            Number of data points to generate.

        randomize: *bool*, optional
            Whether or not to add random noise to the data. Default is False; i.e., perfect data.

        fractional_uncertainty: *float*, optional
            fractional measurement uncertainty in the fluxes. Used to calculate error bars.

        output_file: *str*, optional
            File to output the fake dataset to. If None (default), data are not saved.

        plot: *bool*, optional
            If True, make a plot of the generated data. (Default is False).

    :return: *None*
    """

    # Setup
    _check_inputs(model=model, source_flux=source_flux, blend_flux=blend_flux)
    if t_range is None:
        t_range = model.parameters.t_0 + 2.5 * model.parameters.t_E * np.array(
            [-1, 1])

    if times is None:
        times = np.linspace(t_range[0], t_range[1], n_pts)

    if band == 'I':
        gamma = 0.44
    elif band == 'V':
        gamma = 0.72
    else:
        gamma = None

    fluxes = source_flux * model.get_magnification(times, gamma=gamma) + blend_flux
    errors = fractional_uncertainty * fluxes

    if randomize:
        new_fluxes = fluxes + np.random.normal(0., errors, len(fluxes))
        errors *= new_fluxes / fluxes
        fluxes = new_fluxes

    (mags, mag_errs) = mm.Utils.get_mag_and_err_from_flux(
        fluxes, errors, zeropoint=18.)

    if output_file is not None:
        data = np.vstack((times, mags, mag_errs)).transpose()
        print(data.shape)
        header = '{0}\n'.format(model)
        header += 'source flux: {0}'.format(source_flux)
        header += 'blend flux: {0}'.format(blend_flux)
        if randomize:
            header += 'Data have been randomized.'
        else:
            header += 'Data are perfect.'

        np.savetxt(output_file, data, header=header)

    if plot:
        mod_source_flux = mm.Utils.get_flux_from_mag(
            mm.Utils.get_mag_from_flux(source_flux, zeropoint=18.))
        mod_blend_flux = mm.Utils.get_flux_from_mag(
            mm.Utils.get_mag_from_flux(blend_flux, zeropoint=18.))

        plt.figure()
        model.plot_lc(
            t_range=t_range, source_flux=mod_source_flux,
            blend_flux=mod_blend_flux)
        plt.errorbar(times, mags, yerr=mag_errs, fmt='o')
        plt.gca().minorticks_on()
        plt.show()


def _check_inputs(model=None, source_flux=None, blend_flux=None):
    if model is None:
        raise ValueError('model must be defined.')

    if source_flux is None:
        raise ValueError('source_flux must be defined.')

    if blend_flux is None:
        raise ValueError('blend_flux must be defined.')

    if not (isinstance(model, (mm.Model))):
        raise TypeError('model must be an MulensModel.Model object')

    if not (isinstance(source_flux, (float))):
        raise TypeError('source_flux must be a float.')

    if not (isinstance(blend_flux, (float))):
        raise TypeError('blend_flux must be a float.')


if __name__ == '__main__':
    parameter_sets = [
        {'fluxes': [mm.Utils.get_flux_from_mag(19.9155+4), 0.1], 'n_pts': 3000, 'n_tE': 10.,
         'fractional_uncertainty': 0.01, 'band': 'I'},
        {'fluxes': [2.35, 0.67], 'n_pts': 300, 'n_tE': 0.1,
         'fractional_uncertainty': 0.05, 'band': 'V'}]
    #parameter_sets = [
    #    {'fluxes': [1.2, 0.1], 'n_pts': 250, 'n_tE': 5., 'fractional_uncertainty': 0.01}]
    #model = mm.Model(
    #    {'t_0': 2451697., 'u_0': 0.005, 't_E': 25.2, 'rho': 0.0067},
    #    coords='18:00:00 -30:00:00')
    model = mm.Model(
        {'t_0': 2458310.777, 'u_0': 0.006689, 't_E': 15.886, 'rho': 0.0001},
        coords='17:59:10.26 -27:50:06.3') # OB181185
    model.set_magnification_methods([
        model.parameters.t_0 - 1.5, 'finite_source_LD_Yoo04', model.parameters.t_0 + 1.5])
    #dir = os.path.join(sfit_minimizer.DATA_PATH, 'MMTest')
    dir = './'
    for i, parameters in enumerate(parameter_sets):
        output_file = os.path.join(dir, 'phot_1185_{0}.dat'.format(i + 1))
        #output_file = None
        t_range = model.parameters.t_0 + parameters['n_tE'] * model.parameters.t_E * np.array([-1, 1])
        generate_fake_data(
            model=model, source_flux=parameters['fluxes'][0],
            blend_flux=parameters['fluxes'][1], band=parameters['band'],
            t_range=t_range, n_pts=parameters['n_pts'],
            fractional_uncertainty=parameters['fractional_uncertainty'],
            randomize=False, output_file=output_file, plot=True)
