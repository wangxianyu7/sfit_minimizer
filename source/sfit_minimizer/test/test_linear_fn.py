"""
Assuming a linear function of the form "f = a0 + a1 * x + a2 * x^2 ...", check
that the minimizer works as expected.
"""
import numpy as np


class LinearFunction():
    """
    Calculate the properties of a linear function with respect to the data using
    the sfit minimizer algorithm --> will calculate the next step, which can be
    used in a minimizer.

    Arguments:
        data: array-like with shape (N, 3)
            data[:, 0] = x values
            data[:, 1] = y values
            data[:, 2] = y errors

        theta: tuple, list, np.array
            vector with initial values of parameters to be fitted.
            F = sum_i( theta[i] * x**i )

    """

    def __init__(self, data, theta):
        self.data = data
        self.theta = theta
        self._reset_all()

    def _reset_all(self):
        self._ymod = None
        self._res = None
        self._chi2 = None
        self._df = None
        self._dchi2 = None
        self._chi2_gradient = None
        self._dvec = None
        self._bmat = None
        self._cmat = None
        self._sigmas = None
        self._step = None

    def update_all(self, theta0):
        """
        For a new set of model parameters theta0, recalculate all of the
        data properties with respect to the new model.

        :param theta0:
        :return:
        """
        self.theta = theta0
        self._reset_all()
        self.calc_model()
        self.calc_res()
        self.calc_chi2()
        self.calc_df()
        self.calc_dchi2()
        self.calc_dvec()
        self.calc_bmat()
        self.calc_cmat()
        self.calc_step()

    def calc_model(self):
        pass

    def calc_res(self):
        pass

    def calc_chi2(self):
        pass

    def calc_df(self):
        pass

    def calc_dchi2(self):
        pass

    def calc_dvec(self):
        pass

    def calc_bmat(self):
        pass

    def calc_cmat(self):
        pass

    def calc_step(self):
        pass

    def calc_sigmas(self):
        pass

    @property
    def ymod(self):
        """ model values for data[:, 0]. """
        return self._ymod

    @property
    def res(self):
        """ residuals of the model from the data data = y - ymod"""
        return self._res

    @property
    def chi2(self):
        """ chi2 of the model with respect to the data"""
        return self._chi2

    @property
    def df(self):
        """numerical partial derivatives of the *fitted function F* with respect
        to the parameters, calculated at each data point.
        shape = (len(theta), len(data))"""
        return self._df

    @property
    def dchi2(self):
        """numerical partial derivatives of the *chi2* with respect
        to the parameters, calculated at each data point.
        shape = (len(theta), len(data))"""
        return self._dchi2

    @property
    def chi2_gradient(self):
        """value of the chi2 gradient. shape = ( len(theta) ) """
        return self._chi2_gradient

    @property
    def dvec(self):
        """
        d vector from sfit, used to calculate the step size.
        shape = ( len(theta) )
        """
        return self._dvec

    @property
    def bmat(self):
        """
        b matrix from sfit. invert to find c matrix used for stepping and
        finding sigmas.
        shape = ( len(theta), len(theta) )
        """
        return self._bmat

    @property
    def cmat(self):
        """c matrix from sfit. shape = ( len(theta), len(theta) )"""
        return self._cmat

    @property
    def sigmas(self):
        """Uncertainty in each parameter theta"""
        if self._sigmas is None:
            self.calc_sigmas()

        return self._sigmas

    @property
    def step(self):
        """Step for each parameter theta"""
        return self._step

    # Chi2 of the line
    def get_chi2(self, theta):
        ymod = []
        for i in range(len(theta)):
            ymod.append(theta[i] * self.data[:, 0] ** i)

        ymod = np.array(ymod)
        print(ymod.shape)
        ymod = np.sum(ymod, axis=1)
        print(ymod.shape)

        self.res = ymod - self.data[:, 1]
        self.chi2 = np.sum(self.res ** 2 / self.data[:, 2] ** 2)

        return self.chi2

    # Derivative of a line
    def get_partials(self, theta):
        df = []
        for i in range(len(theta)):
            df.append(self.data[:, 0] ** i)

        self.df = np.array(df)
        print(self.df.shape)

        return np.sum(self.df, axis=1)

    # Derivative of the chi2 of a line
    def get_chi2_gradient(self, theta):
        chi2_gradient = []
        for i in range(len(theta)):
            chi2_gradient.append(
                -2 * self.res * self.df[i] / self.data[:, 2] ** 2)

        self.chi2_gradient = np.array(chi2_gradient)
        print(self.chi2_gradient.shape)

        return np.sum(self.chi2_gradient, axis=1)

    def get_dvec(self):
        self.dvec = np.sum(self.chi2_gradient, axis=0) / 2.
        print(self.dvec.shape)

        return self.dvec

    def get_bmat(self):
        self.bmat = np.zeros((len(self.theta), len(self.theta)))


def test_1():
    """
    Fit a line to perfect data: y = m * x + b
    :return:
    """


    # Fake data

    # Wrong initial condition

    # Compare results
    # Check b, c, d matrices
    # Check true parameters are recovered.