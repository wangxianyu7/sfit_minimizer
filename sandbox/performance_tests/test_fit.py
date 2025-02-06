import MulensModel as mm
import scipy.optimize as op
import os.path

def chi2_fun(theta, event, parameters_to_fit):
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
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])

    # After that, calculating chi2 is trivial:
    return event.get_chi2()


def jacobian(theta, event, parameters_to_fit):
    """
    Calculate chi^2 gradient (also called Jacobian).

    Note: this implementation is robust but possibly inefficient. If
    chi2_fun() is ALWAYS called before jacobian with the same parameters,
    there is no need to set the parameters in event.model; also,
    event.calculate_chi2_gradient() can be used instead (which avoids fitting
    for the fluxes twice).
    """
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])

    return event.get_chi2_gradient(parameters_to_fit)

# Read in the data file
dir = 'KB171219'
field = 31
datasets = []
for site in ['C', 'A', 'S']:
    data = mm.MulensData(
        file_name=os.path.join(
            dir, 'KMT{0}{1:02}_I.pysis'.format(site, field)),
    usecols=[0,3,4])
    datasets.append(data)

# Initialize the fit
parameters_to_fit = ["t_0", "u_0", "t_E"]
# Approximate values of the parameters are needed:
t_0 = 7951.3
u_0 = 0.1
t_E = 18.
model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Link the data and the model
event = mm.Event(datasets=datasets, model=model)
print('Initial Trial\n{0}'.format(event.model.parameters))
print('Initial chi2: {0:12.2f}'.format(event.get_chi2()))

# Find the best-fit parameters
for method in ['Nelder-Mead', 'Newton-CG', 'BFGS']:
    initial_guess = [t_0, u_0, t_E]
    result = op.minimize(
        chi2_fun, x0=initial_guess, args=(event, parameters_to_fit),
        jac=jacobian,
        method=method)
    print(result.x)
    (fit_t_0, fit_u_0, fit_t_E) = result.x

    # Save the best-fit parameters
    chi2 = chi2_fun(result.x, event, parameters_to_fit)

    # Output the fit parameters
    msg = 'Best Fit: t_0 = {0:12.5f}, u_0 = {1:6.4f}, t_E = {2:8.3f}'
    print(msg.format(fit_t_0, fit_u_0, fit_t_E))
    print('Chi2 = {0:12.2f}'.format(chi2))
    print('scipy.optimize.minimize result:')
    print(result)
