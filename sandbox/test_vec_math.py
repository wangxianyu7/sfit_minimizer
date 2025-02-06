import numpy as np
import timeit


class TestVec(object):
    def __init__(self):
        self.n = 1000
        self.npar = 5
        self.iterations = 1000
        self.data_k = np.random.uniform(size=(self.npar, self.n))
        self.err_k = np.random.uniform(size=(self.n))

    def loop_version(self):
        raise NotImplementedError()

    def numpy_version(self):
        raise NotImplementedError()

    def test(self):
        mat_loop = self.loop_version()
        mat_numpy = self.numpy_version()

        print('{0}:'.format(type(self)))
        print('compare shape: ')
        print(mat_loop.shape, mat_numpy.shape)
        np.testing.assert_almost_equal(mat_loop, mat_numpy)
        print('loop: ')
        print(timeit.timeit(self.loop_version, number=self.iterations))
        print('numpy: ')
        print(timeit.timeit(self.numpy_version, number=self.iterations))


class BMatTest(TestVec):

    def __init__(self):
        self.n = 1000
        self.npar = 5
        self.iterations = 1000
        self.data_k = np.random.uniform(size=(self.npar, self.n))
        self.err_k = np.random.uniform(size=(self.n))

    def loop_version(self):
        bmat = np.zeros(shape=(self.npar, self.npar))
        for i in range(self.npar):
            for j in range(self.npar):
                bmat[i, j] = np.sum(self.data_k[i] * self.data_k[j] / self.err_k ** 2)

        return bmat

    def numpy_version(self):
        bmat = np.inner(self.data_k, self.data_k / self.err_k**2)

        return bmat

class StepTest(TestVec):

    def __init__(self):
        TestVec.__init__(self)
        self.theta = range(self.npar)
        self.cmat = np.random.uniform(size=(self.npar, self.npar))
        self.dvec = np.random.uniform(size=self.npar)

    def loop_version(self):
        self.step = np.zeros(len(self.theta))
        for i in range(len(self.theta)):
            for j in range(len(self.theta)):
                self.step[i] += self.cmat[i, j] * self.dvec[j]

        return self.step

    def numpy_version(self):
        self.step = np.sum(self.cmat * self.dvec, axis=1)

        return self.step

if __name__ == '__main__':
    bmat = StepTest()
    bmat.test()
