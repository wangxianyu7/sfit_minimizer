import sfit_minimize
import LinearFunction
import scipy.optimize as op
import matplotlib.pyplot as plt

def chi2_fun(theta, my_func):
    """
    for given event set attributes from parameters_to_fit (list of
    str) to values from theta list
    """
    my_func.theta = theta
    my_func.reset_all()

    return event.get_chi2()

def partials_fun(theta, my_func):
    # Is this enough? How does minimize access the rest of the information 
    # from my_func? e.g., res
    my_func.theta = theta
    return my_func.df

data = np.loadtxt('../data/test_data_10000pts_Poisson.txt', skiprows=2)
initial_guess = [] # Wrong initial condition
my_func = LinearFunction(data=data, theta=theta)

result = op.minimize(
    chi2_fun, method=sfit_minimize.minimize, x0=initial_guess,
    args=(ev, parameters_to_fit),
    partials_fun=partials_fun, tol=1e-3, options={'step': 'adaptive'})

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
