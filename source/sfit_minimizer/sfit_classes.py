import numpy as np


class SFitResults(object):
    """
    Results of the minimization algorithm. Modeled after scipy.optimize.OptimizeResult.

    Arguments:
        func: :py:class:`~SFitFunction`

    Keywords:
        success: *bool*
            see :py:attr:`~success`

        msg: *str*
            see :py:attr:`~msg`

    """

    def __init__(self, func, success=None, msg=None):
        self._x = func.theta
        self._sigmas = func.get_sigmas()
        self._success = success
        self._msg = msg

    def __repr__(self):
        output_str = 'x = {0}\n'.format(self.x)
        output_str += 'sigmas = {0}\n'.format(self.sigmas)
        output_str += 'success = {0}\n'.format(self.success)
        output_str += 'msg = {0}'.format(self.msg)
        return output_str

    @property
    def x(self):
        """ Best-fit (last calculated) value of the parameters. """
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @property
    def sigmas(self):
        """ Uncertainties for :py:attr:`~x`"""
        return self._sigmas

    @sigmas.setter
    def sigmas(self, value):
        self._sigmas = value

    @property
    def success(self):
        """ True = algorithm completed successfully. False = algorithm ran into a limit or failure mode."""
        return self._success

    @success.setter
    def success(self, value):
        self._success = value

    @property
    def msg(self):
        """ If :py:attr:`~success` == False, some information about the failure mode."""
        return self._msg

    @msg.setter
    def msg(self, value):
        self._msg = value


class SFitFunction(object):
    """
    Main class for sfit minimization routine. Establishes all the necessary functions for executing the A. Gould
    algorithm. At a minimum, the user should specify:

      EITHER :py:func:`~calc_model` OR :py:func:`~calc_residuals`

      AND :py:func:`~calc_df`.

    Arguments:
        data: *np.array* with shape (N, 3)
            an array of the data to be fitted with columns (x, y, yerr).

        theta: *list* of length M
            trial values for the parameters of the function to be fit.

    """

    def __init__(self, data=None, theta=None):
        if isinstance(data, np.ndarray):
            if data.shape[1] == 3:
                self.data = data
            else:
                raise ValueError('data must have shape (N, 3). Shape: {0}'.format(data.shape))

        else:
            raise TypeError(
                'data must be an np.array object with shape (N, 3). Type: {0}'.format(type(data))
            )

        print('data checker implemented')

        self.theta = theta
        self._reset_all()

    def _reset_all(self):
        self._ymod = None
        self._residuals = None
        self._chi2 = None
        self._df = None
        self._dchi2 = None
        self._dvec = None
        self._bmat = None
        self._cmat = None
        self._sigmas = None
        self._step = None

    def update_all(self, theta0=None):
        """
        For a new set of model parameters theta0, recalculate all of the
        data properties with respect to the new model.
        """
        if theta0 is not None:
            self._reset_all()
            self.theta = theta0

        self.calc_model()
        self.calc_residuals()
        self.calc_chi2()
        self.calc_df()
        self.calc_dchi2()
        self.calc_dvec()
        self.calc_bmat()
        self.calc_cmat()
        self.calc_step()

    # Calculation functions
    def calc_model(self):
        """
        FUNCTION SPECIFIC: Calculate expected values of the model. May
        be explicitly defined for a specific class that inherits
        SFitFunction.

        sets self.ymod =
            *np.array* of shape (N), where each element k is the value of the model evaluated at data[k, 0].

        """
        raise NotImplementedError('User must define calc_model() for a child class of SFitFunction.')

    def calc_residuals(self):
        """
        Difference between the data and the model. Either this or
        calc_model() should be explicitly defined for a specific
        class that inherits SFitFunction.

        sets self.residuals =
            *np.array* of shape (N), where each element k = model[k] - data[k, 1].

        """
        if self.ymod is None:
            self.calc_model()

        self.residuals = self.ymod - self.data[:, 1]

    def calc_chi2(self):
        """
        Calculate the chi2 of the data relative to the model.

        sets self.chi2 = Sum_k (residuals[k]^2) / data[k, 2]^2
        """
        if self.residuals is None:
            self.calc_residuals()

        self._chi2 = np.sum(self.residuals ** 2 / self.data[:, 2] ** 2)

    def calc_df(self):
        """
        FUNCTION SPECIFIC: Calculate the derivatives of the fitting
        function relative to each parameter theta and store as
        self.df.  Should be explicitly defined for a specific class
        that inherits SFitFunction.

        sets self.df =
            *np.array* of shape (M, N) where df[i, k] returns the derivative of the fitting function with respect to
            parameter i evaluated at point k.

        """
        raise NotImplementedError('User must define calc_model() for a child class of SFitFunction.')

    def calc_dchi2(self):
        """Calculate the gradient of chi2 at each datapoint"""

        if self.df is None:
            self.calc_df()

        if self.residuals is None:
            self.calc_residuals()

        chi2_gradient = []
        for i in range(len(self.theta)):
            chi2_gradient.append(
                -2 * self.residuals * self.df[i] / self.data[:, 2]**2)

        self._dchi2 = np.array(chi2_gradient)

    def calc_dvec(self):
        """
        Calculate the d vector and store it as self.dvec.
        d_i = d_i + (partial chi2/partial a_i)/2
        """
        if self.dchi2 is None:
            self.calc_dchi2()

        self._dvec = np.sum(self.dchi2, axis=1) / 2.

    def calc_bmat(self):
        """
        Calculate the b matrix and store it as self.bmat.
        b_ij = b_ij + (partial F/partial a_i)(partial F/partial a_j)/error^2
        """
        if self.df is None:
            self.calc_df()

        self._bmat = np.inner(self.df, self.df / self.data[:, 2]**2)
        # self._bmat = np.zeros((len(self.theta), len(self.theta)))
        # for i in range(len(self.theta)):
        #     for j in range(len(self.theta)):
        #         self._bmat[i, j] = np.sum(
        #             self.df[i] * self.df[j] / self.data[:, 2]**2)

    def calc_cmat(self):
        """Invert the b matrix to find the c matrix."""
        if self.bmat is None:
            self.calc_bmat()

        self._cmat = np.linalg.inv(self.bmat)

    def calc_step(self):
        """
        Calculate the full proposed step size.
        Delta a_i = sum_j c_ij * d_j.
        """
        if self.cmat is None:
            self.calc_cmat()

        if self.dvec is None:
            self.calc_dvec()

        self._step = np.sum(self.cmat * self.dvec, axis=1)
        # self._step = np.zeros(len(self.theta))
        # for i in range(len(self.theta)):
        #     for j in range(len(self.theta)):
        #         self._step[i] += self.cmat[i, j] * self.dvec[j]

    def calc_sigmas(self):
        """
        Calculate the uncertainty in each parameter theta.
        sigma_a_i = sqrt(c_ii)
        """
        if self.cmat is None:
            self.calc_cmat()

        # self._sigmas = np.zeros(len(self.theta))
        # for i in range(len(self.theta)):
        #     self._sigmas[i] = np.sqrt(self.cmat[i, i])

        self._sigmas = np.sqrt(np.diagonal(self.cmat))

    # Properties
    @property
    def theta(self):
        """
        *list* or *np.array*

        Parameters of the function being fit to the data.
        """
        return self._theta

    @theta.setter
    def theta(self, value):
        self._reset_all()
        if value is None:
            self._theta = value
        else:
            if isinstance(value, (list, np.ndarray)):
                self._theta = value
            else:
                raise TypeError(
                    'theta must be either list or np.array. Type: {0}'.format(type(value))
                )

    @property
    def ymod(self):
        """
        *np.array* of shape (N)

        model values for data[:, 0].
        """
        return self._ymod

    @ymod.setter
    def ymod(self, value):
        if isinstance(value, np.ndarray):
            if len(value) == len(self.data):
                self._ymod = value
            else:
                raise ValueError(
                    'ymod should have the same length as data. data.shape: {0}, ymod.shape: {1}'.format(
                        self.data.shape, value.shape))

        else:
            raise TypeError('ymod must be an np.array object. Type: {0}'.format(type(value)))

    @property
    def residuals(self):
        """
        *np.array* of shape (N)

        residuals of the model from the data: residuals = y - ymod
        """
        return self._residuals

    @residuals.setter
    def residuals(self, value):
        if isinstance(value, np.ndarray):
            if len(value) == len(self.data):
                self._residuals = value
            else:
                raise ValueError(
                    'residuals should have the same length as data. data.shape: {0}, residuals.shape: {1}'.format(
                        self.data.shape, value.shape))
        else:
            raise TypeError('residuals must be an np.array object. Type: {0}'.format(type(value)))

    @property
    def chi2(self):
        """ chi2 of the model with respect to the data"""
        return self._chi2

    @property
    def df(self):
        """
        *np.array* of shape (M, N)

        Partial derivatives of the *fitted function F* with respect
        to the parameters, calculated at each data point.
        shape = (len(theta), len(data))"""
        if self._df is None:
            self.calc_df()

        return self._df

    @df.setter
    def df(self, value):
        if isinstance(value, np.ndarray):
            if value.shape == (len(self.theta), len(self.data)):
                self._df = value
            else:
                raise ValueError(
                    'df should have the shape (M, N) where M is the number of parameters and N is the number of data' +
                    ' points. len(theta): {2}, data.shape: {0}, df.shape: {1}'.format(
                        self.data.shape, value.shape, len(self.theta)))
        else:
            raise TypeError('df must be an np.array object. Type: {0}'.format(type(value)))

    @property
    def dchi2(self):
        """Partial derivatives of the *chi2* with respect
        to the parameters, calculated at each data point.
        shape = (len(theta), len(data))"""
        return self._dchi2

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
    def get_chi2(self):
        """Calculate and return the chi2 of the data w.r.t. the model"""
        self.calc_chi2()

        return self.chi2

    def get_partials(self):
        """Calculate and return the partial derivatives of the fitting function"""
        self.calc_df()

        return np.sum(self.df, axis=1)

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
