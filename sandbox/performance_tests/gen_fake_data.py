"""
Generate a fake dataset tha scipy.minimize fails to fit.
"""
import MulensModel as mm
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

kmtdata = np.genfromtxt(
    'KMTC14_I.pysis.txt', dtype=None, names=['times'], usecols=[0])
times = kmtdata['times']

ulens_keys = ['t_0', 'u_0', 't_E', 'rho', 'pi_E_N', 'pi_E_E']
def make_event(params):
    ulens_params = {}
    for key in ulens_keys:
        if key in params.keys():
            ulens_params[key] = params[key]

    model = mm.Model(ulens_params)

    magnification = model.get_magnification(times)

    fluxes = params['f_source'] * magnification + params['f_blend']
    mags = mm.Utils.get_mag_from_flux(fluxes)
    fractional_uncertainty = 0.01 * (mags - 15.)
    fractional_uncertainty[fractional_uncertainty < 0.001] = 0.001
    errors = fractional_uncertainty * fluxes
    new_fluxes = fluxes + np.random.normal(0., errors, len(fluxes))
    errors *= new_fluxes / fluxes
    fluxes = new_fluxes
    (mags, mag_errs) = mm.Utils.get_mag_and_err_from_flux(
        fluxes, errors, zeropoint=18.)
    np.savetxt('fake_data.dat', np.array([times, mags, mag_errs]).transpose())
    fake_data = mm.MulensData(
        [times, mags, mag_errs], phot_fmt='mag')

    event = mm.Event(datasets=fake_data, model=model)
    #event.plot()
    #fake_data.plot()
    plt.show()

    return event


def get_initial_guess(params):
    trial = []
    for key in ulens_keys:
        if key in params.keys():
            if key == 't_0':
                trial.append(params[key] + (np.random.normal(0, 3.)))
            else:
                trial.append(params[key] * (1 + np.random.normal(0, 0.5)))

    return trial

def get_parameters_to_fit(params):
    ulens_params = []
    for key in ulens_keys:
        if key in params.keys():
            ulens_params.append(key)

    return ulens_params


def evaluate_results(results):
    print(results)

def chi2_func(theta, parameters_to_fit, event):
    """
    Calculate chi2 for given values of parameters

    Keywords :
        theta: *np.ndarray*
            Vector of parameter values, e.g.,
            `np.array([5380., 0.5, 20.])`.

        parameters_to_fit: *list* of *str*
            List of names of parameters corresponding to theta, e.g.,
            `['t_0', 'u_0', 't_E']`.

        event: *MulensModel.Event*
            Event which has datasets for which chi2 will be calculated.

    Returns :
        chi2: *float*
            Chi2 value for given model parameters.
    """
    # First we have to change the values of parameters in
    # event.model.parameters to values given by theta.
    for (parameter, value) in zip(parameters_to_fit, theta):
        setattr(event.model.parameters, parameter, value)

    # After that, calculating chi2 is trivial:
    return event.get_chi2()

params = {'t_0': 8303., 'u_0': 0.6, 't_E': 10.,
          'f_source': 1.0, 'f_blend': 2.0}

print('input', params)
event = make_event(params)
initial_guess = get_initial_guess(params)
parameters_to_fit = get_parameters_to_fit(params)
print('initial guess', initial_guess)
results = minimize(chi2_func, x0=initial_guess, args=(parameters_to_fit, event),
    method='Nelder-Mead')
print('final')
print(event.model)
print(event.get_chi2())
evaluate_results(results)
event.plot()
plt.show()