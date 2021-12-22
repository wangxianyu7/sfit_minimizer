"""
Assuming a linear function of the form "f = a0 + a1 * x + a2 * x^2 ...", check
that the minimizer works as expected.
"""
import numpy as np
import os.path
import sfit_minimizer


class LinearFunction(sfit_minimizer.SFitFunction):

    def __init__(self, data=None, theta=None):
        sfit_minimizer.SFitFunction.__init__(self, data=data, theta=theta)

    def calc_model(self):
        """Calculate expected values of the model"""
        ymod = []
        for i in range(len(self.theta)):
            ymod.append(self.theta[i] * self.data[:, 0] ** i)

        ymod = np.array(ymod)
        self.ymod = np.sum(ymod, axis=0)

    def calc_df(self):
        """Calculate the derivatives of the fitting function and store as
        self.df."""
        df = []
        for i in range(len(self.theta)):
            df.append(self.data[:, 0] ** i)

        self.df = np.array(df)


def test_perfect_chi2():
    """chi2 should be 0"""
    data = np.loadtxt(os.path.join(sfit_minimizer.DATA_PATH, 'test_data_10perfect.txt'), skiprows=2)
    theta = [3, 2]
    my_func = LinearFunction(data=data, theta=theta)
    my_func.calc_chi2()

    tol = 1e-6
    assert np.sum(my_func.res) < tol
    assert my_func.chi2 < tol


# Test the SFitFunction calculations produce the expected results
def test_fit():
    """
    Fit a line to perfect data: y = m * x + b
    :return:
    """
    # Fake data
    data = np.loadtxt(os.path.join(sfit_minimizer.DATA_PATH, 'test_data_10000pts_Poisson.txt'), skiprows=2)

    # Wrong initial condition
    theta = [4, 2.1]

    # Compare results
    my_func = LinearFunction(data=data, theta=theta)
    my_func.update_all()
    print('step', my_func.step)
    print('df.shape', my_func.df.shape)

    # Check ymod, res
    y = data[:, 0] * theta[1] + theta[0]
    res = y - data[:, 1]
    np.testing.assert_almost_equal(my_func.ymod, y)
    np.testing.assert_almost_equal(my_func.res, res)

    # Check b, c, d matrices
    b = np.zeros((2, 2))
    b[0, 0] = np.sum(1. / data[:, 2] ** 2)
    b[0, 1] = np.sum(data[:, 0] / data[:, 2] ** 2)
    b[1, 0] = b[0, 1]
    b[1, 1] = np.sum(data[:, 0] ** 2 / data[:, 2] ** 2)
    c = np.linalg.inv(b)
    d = np.zeros(2)
    d[0] = np.sum(-res / data[:, 2] ** 2)
    d[1] = np.sum(-res * data[:, 0] / data[:, 2] ** 2)

    np.testing.assert_almost_equal(my_func.bmat, b)
    np.testing.assert_almost_equal(my_func.cmat, c)
    np.testing.assert_almost_equal(my_func.dvec, d)

    # Check true parameters are recovered.
    theta_new = theta + my_func.step
    sigmas = my_func.get_sigmas()
    try:
        assert np.abs(theta_new[0] - 3.) < 3. * sigmas[0]
        assert np.abs(theta_new[1] - 2.) < 3. * sigmas[1]
    except AssertionError:
        return my_func, theta_new


# Test that the minimize function produces the expected results
def test_minimize():
    data = np.loadtxt(os.path.join(sfit_minimizer.DATA_PATH, 'test_data_10000pts_Poisson.txt'), skiprows=2)
    initial_guess = [4, 2.1]  # Wrong initial condition
    my_func = LinearFunction(data=data)

    result = sfit_minimizer.minimize(
        my_func, x0=initial_guess, tol=1e-7,
        options={'step': 'adaptive'}, verbose=False)

    my_func.theta = result.x
    np.testing.assert_almost_equal(my_func.get_chi2(), 10865.52999, decimal=3)
