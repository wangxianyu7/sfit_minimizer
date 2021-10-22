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

    def update_all(self, theta0=None):
        """
        For a new set of model parameters theta0, recalculate all of the
        data properties with respect to the new model.

        :param theta0:
        :return:
        """
        if theta0 is not None:
            self._reset_all()
            self.theta = theta0

        self.calc_model()
        self.calc_res()
        self.calc_chi2()
        self.calc_df()
        self.calc_dchi2()
        self.calc_dvec()
        self.calc_bmat()
        self.calc_cmat()
        self.calc_step()

    # Calculation functions
    def calc_model(self):
        """Calculate expected values of the model"""
        ymod = []
        for i in range(len(self.theta)):
            ymod.append(self.theta[i] * self.data[:, 0] ** i)

        ymod = np.array(ymod)
        #print('calc_model(): shape of ymod (expect len(theta) x N)', ymod.shape)
        self._ymod = np.sum(ymod, axis=0)
        #print('calc_model(): shape of ymod after sum (expect N)', self.ymod.shape)

    def calc_res(self):
        if self.ymod is None:
            self.calc_model()

        self._res = self.ymod - self.data[:, 1]

    def calc_chi2(self):
        if self.res is None:
            self.calc_res()

        self._chi2 = np.sum(self.res ** 2 / self.data[:, 2] ** 2)

    def calc_df(self):
        """Calculate the derivatives of the fitting function and store as self.df."""
        df = []
        for i in range(len(self.theta)):
            df.append(self.data[:, 0] ** i)

        self._df = np.array(df)
        #print('calc_df(): shape of df (expect len(theta) x N)', self.df.shape)

    def calc_dchi2(self):
        """Calculate the gradient of chi2 at each datapoint"""
        if self.df is None:
            self.calc_df()

        if self.res is None:
            self.calc_res()

        chi2_gradient = []
        for i in range(len(self.theta)):
            chi2_gradient.append(
                -2 * self.res * self.df[i] / self.data[:, 2]**2)

        self._dchi2 = np.array(chi2_gradient)
        #print(
        #    'calc_dchi2(): shape of dchi2 (expect len(theta) x N)',
        #    self.dchi2.shape)

    def calc_chi2_gradient(self):
        """Calculate the gradient of chi2 summed over ALL datapoints"""
        if self.dchi2 is None:
            self.calc_dchi2()

        self._chi2_gradient = np.sum(self.dchi2, axis=1)
        #print(
        #    'calc_chi2_gradient(): shape of chi2_gradient, expect len(theta)',
        #    self.chi2_gradient.shape)

    def calc_dvec(self):
        """
        Calculate the d vector and store it as self.dvec.
        d_i = d_i + (partial chi2/partial a_i)/2
        """
        if self.dchi2 is None:
            self.calc_dchi2()

        self._dvec = np.sum(self.dchi2, axis=1) / 2.
        #print('calc_dvec(): shape of dvec, expect len(theta)', self.dvec.shape)

    def calc_bmat(self):
        """
        Calculate the b matrix and store it as self.bmat.
        b_ij = b_ij + (partial F/partial a_i)(partial F/partial a_j)/error^2
        """
        if self.df is None:
            self.calc_df()

        self._bmat = np.zeros((len(self.theta), len(self.theta)))
        for i in range(len(self.theta)):
            for j in range(len(self.theta)):
                self._bmat[i, j] = np.sum(
                    self.df[i] * self.df[j] / self.data[:, 2]**2)

    def calc_cmat(self):
        """Invert the b matrix to find the c matrix."""
        if self.bmat is None:
            self.calc_bmat()

        self._cmat = np.invert(self.bmat)

    def calc_step(self):
        """
        Calculate the full proposed step size.
        Delta a_i = sum_j c_ij * d_j.
        """
        if self.cmat is None:
            self.calc_cmat()

        if self.dvec is None:
            self.calc_dvec()

        self._step = np.zeros(len(self.theta))
        for i in range(len(self.theta)):
            for j in range(len(self.theta)):
                self._step[i] += self.cmat[i, j] * self.dvec[j]

    def calc_sigmas(self):
        """
        Calculate the uncertainty in each parameter theta.
        sigma_a_i = sqrt(c_ii)
        """
        if self.cmat is None:
            self.calc_cmat()

        self._sigmas = np.zeros(len(self.theta))
        for i in range(len(self.theta)):
            self._sigmas[i] = np.sqrt(self.cmat[i, i])

    # Properties
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
        return self._sigmas

    @property
    def step(self):
        """Step for each parameter theta"""
        return self._step

    # Retrieval Functions
    def get_chi2(self, theta):
        """Calculate and return the chi2 of the data w.r.t. the model"""
        self.calc_chi2()

        return self.chi2


    def get_partials(self, theta):
        """Calculate and return the partial derivatives of the fitting function"""
        self.calc_df()

        return np.sum(self.df, axis=1)

    def get_chi2_gradient(self, theta):
        """Calculate and return the derivative of the chi2 of a line."""
        self.calc_chi2_gradient()

        return np.sum(self.chi2_gradient, axis=1)

    def get_dvec(self):
        """Calculate and return the d vector"""
        self.calc_dvec()

        return self.dvec

    def get_bmat(self):
        """Calculate and return the b matrix"""
        self.calc_bmat()

        return self.bmat

    def get_cmat(self):
        """Calculate and return the c matrix"""
        self.calc_cmat()

        return self.cmat

    def get_step(self):
        """Calculate and return the full step for theta_new = theta_old + step."""
        self.calc_step()

        return self.step

    def get_sigmas(self):
        """Calculate and return the uncertainty in each parameter theta."""
        self.calc_sigmas()

        return self.sigmas


def test_perfect_chi2():
    """chi2 should be 0"""
    data = np.loadtxt('test_data_10perfect.txt', skiprows=2)
    theta = [3, 2]
    my_func = LinearFunction(data=data, theta=theta)
    my_func.update_all()

    tol = 1e-6
    assert np.sum(my_func.res) < tol
    assert my_func.chi2 < tol

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

def test_2():
    """
    Fit a line to imperfect data: y = m * x + b
    :return:
    """
    # Fake data

    # Wrong initial condition

    # Compare results
    # Check b, c, d matrices
    # Check true parameters are recovered.

if __name__ == '__main__':
    test_perfect_chi2()