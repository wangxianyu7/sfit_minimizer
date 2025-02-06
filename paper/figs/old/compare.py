import matplotlib.pyplot as plt
import numpy as np
from plotting import default_settings


class FitResults():

    def __init__(self, filename):
        self.sfit = []
        self.sfit_sigmas = []
        self.neldar_mead = []
        self.bfgs = []
        self.newton_cg = []
        results = open(filename, 'r')
        for line in results.readlines():
            if line[0] == 'K':
                name = line
            else:
                elements = line.split()
                values = [float(value) for value in elements[1:]]
                new_entry = np.hstack((name, values))
                if elements[0].strip() == 'sigmas':
                    self.sfit_sigmas.append(new_entry)
                elif elements[0] == 'SFit':
                    self.sfit.append(new_entry)
                elif elements[0] == 'Neldar-Mead':
                    self.neldar_mead.append(new_entry)
                elif elements[0] == 'Newton-CG':
                    self.newton_cg.append(new_entry)
                elif elements[0] == 'BFGS':
                    self.bfgs.append(new_entry)

        results.close()

        self.sfit = np.array(self.sfit)
        self.sfit_sigmas = np.array(self.sfit_sigmas)
        self.bfgs = np.array(self.bfgs)
        self.neldar_mead = np.array(self.neldar_mead)
        self.newton_cg = np.array(self.newton_cg)
        #print(sfit.shape, sfit_sigmas.shape)

        self.all_fits = np.dstack((
            self.sfit, self.neldar_mead, self.newton_cg, self.bfgs))
        self.fit_type = ['SFit', 'Neldar-Mead', 'Newton-CG', 'BFGS']
        #print(self.all_fits.shape)
        print('Number of events:', self.all_fits.shape[0])
        print('Number of fits:', self.all_fits.shape[2])
        self._best_index = None
        self._dchi2 = None
        self._reported_success = None
        self._true_success = None


    def count_epic_failures(self, verbose=False):
        count = 0
        for event in self.all_fits:
            #print(event[1, :])
            #print(event.shape)
            success = np.sum(event[1, :].astype(float))
            if success == 0:
                count += 1
                if verbose:
                    print('{0} all fits failed'.format(event[0, 0]))

        print('total complete failures', count)

    def count_false_postives(self, verbose=False):
        # Defined as fit was successful but not best chi2.
        false_positives = np.array([0, 0, 0, 0])
        for i, event in enumerate(self.all_fits):
            dchi2 = self.dchi2[i, :]
            bad = (dchi2 > 0.001) & (event[1, :].astype(float) > 0)
            false_positives += bad.astype(int)

            if verbose:
                print('best', event[:, self.best_index[i]])
                print('success', event[1, :])
                print('dchi2', dchi2)
                print('bad', bad)
                print('false_positives', false_positives)


        print('False positives:')
        for i in range(len(self.fit_type)):
            print(self.fit_type[i], false_positives[i])

    def get_nev_statistics(self):
        print('N func evals:')
        print('{0:12} {1:6} {2:6} {3:6} {4:6} {5:6}'.format(
            'type', 'notnan', 'mean', 'median', 'std', 'max'))
        for i in range(len(self.fit_type)):
            nfev = self.all_fits[:, 2, i].astype(float)
            print('{0:12} {1:6} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f}'.format(
                self.fit_type[i], np.sum(np.isfinite(nfev)),
                np.nanmean(nfev), np.nanmedian(nfev), np.nanstd(nfev),
                np.nanmax(nfev)))

        print('N jac evals:')
        print('{0:12} {1:6} {2:6} {3:6} {4:6} {5:6}'.format(
            'type', 'notnan', 'mean', 'median', 'std', 'max'))
        for i in range(len(self.fit_type)):
            njev = self.all_fits[:, 7, i].astype(float)
            print('{0:12} {1:6} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f}'.format(
                self.fit_type[i], np.sum(np.isfinite(njev)),
                np.nanmean(njev), np.nanmedian(njev), np.nanstd(njev),
                np.nanmax(njev)))

    def calc_dchi2s(self):
        best_chi2 = []
        for i, index in enumerate(self.best_index):
            best_chi2.append(self.all_fits[i, 3, index].astype(float))

        best_chi2 = np.array(best_chi2)
        chi2s = self.all_fits[:, 3, :].astype(float)
        self._dchi2 = np.subtract(chi2s, best_chi2[:, np.newaxis])

    def plot_dchi2(self):
        for max in [10, 100, 1000]:
            plt.figure()
            for i in range(len(self.fit_type)):
                plt.subplot(2, 2, i+1)
                plt.title(self.fit_type[i])
                plt.hist(self.dchi2[:, i], bins=10, range=[0, max])
                plt.yscale('log')
                plt.minorticks_on()

        #plt.tight_layout()

    def tabulate_dchi2(self):
        def print_fits(fits):
            outstr = ''
            for i, fit in enumerate(self.fit_type):
                outstr += '{0:15}'.format(fit)
                for item in fits[i, :]:
                    outstr += ' {0:4}'.format(item.astype(int))

                outstr += '\n'

            print(outstr)

        print('N >Dchi2:')

        tol = [0., 0.1, 1., 10., 100.]
        count_all = np.zeros((len(self.fit_type), len(tol)))
        count_success = np.zeros((len(self.fit_type), len(tol)))
        count_fail = np.zeros((len(self.fit_type), len(tol)))
        for i, event in enumerate(self.all_fits):
            dchi2 = self.dchi2[i, :]
            #print(event.transpose())
            #print(dchi2)
            for j, tolerance in enumerate(tol):
                count_all[:, j] += (dchi2 >= tolerance).astype(int)
                count_success[:, j] += ((dchi2 >= tolerance) &
                                        (event[1, :].astype(float) > 0.5)).astype(int)
                count_fail[:, j] += ((dchi2 >= tolerance) &
                                       (event[1, :].astype(float) < 0.5)).astype(int)

            #print_fits(count_all)


        print('All', tol)
        print_fits(count_all)
        print('Successful Fits', tol)
        print_fits(count_success)
        print('Failed Fits', tol)
        print_fits(count_fail)

    def set_best_index(self, verbose=False):
        index = []
        for event in self.all_fits:
            index.append(np.argmin(event[3, :].astype(float)))
            if verbose:
                print('')
                print(event.transpose())
                print('best:')
                print(event[:, index[-1]])

        self._best_index = index

    def plot_hist(self, fit_type=None, parameter=None, exclude=5, frac=False, **hist_kwargs):
        if parameter == 't_0':
            p = 4
        elif parameter == 'u_0':
            p = 5
        elif parameter == 't_E':
            p = 6

        best_value = []
        for i, index in enumerate(self.best_index):
            best_value.append(self.all_fits[i, p, index].astype(float))

        best_value = np.array(best_value)

        if fit_type == 'SFit':
            data = self.sfit
            f = 0
        elif fit_type == 'BFGS':
            data = self.bfgs
            f = 3
        else:
            raise NotImplementedError(
                'your fit_type has not been implemented:', fit_type)

        true_positives = (self.reported_success & self.true_success[:, f])
        false_positives = (self.reported_success & (~self.true_success[:, f]))
        false_negatives = (self.true_success[:, f] & (~self.reported_success))
        true_negatives = (~self.true_success[:, f]) & (~self.reported_success)
        for name, index, edgecolor in zip(
                ['False Positives', 'False Negatives'],
                [false_positives, false_negatives],
                ['red', 'black']):
            #print(data[index, p])
            #print(best_value[index])
            if np.sum(index) > 0:
                delta = np.subtract(data[index, p].astype(float), best_value[index])
                print(parameter, name, np.min(delta), np.max(delta), np.percentile(delta, [exclude, 100-exclude]))
                #bin_max = np.max(np.abs(np.percentile(delta, [exclude, 100-exclude])))
                if frac:
                    plt.hist(
                        delta / best_value[index], weights=np.ones(delta.shape[0]) / np.sum(index), label=name, facecolor='none',
                        edgecolor=edgecolor,
                        **hist_kwargs)
                else:
                    plt.hist(
                        delta, weights=np.ones(delta.shape[0]) / np.sum(index), label=name, facecolor='none', edgecolor=edgecolor,
                        **hist_kwargs)#, cumulative=True, density=True)

    @property
    def best_index(self):
        if self._best_index is None:
            self.set_best_index()

        return self._best_index

    @property
    def dchi2(self):
        if self._dchi2 is None:
            self.calc_dchi2s()

        return self._dchi2

    @property
    def true_success(self):
        if self._true_success is None:
            self._true_success = self.dchi2 < 0.1

        return self._true_success

    @property
    def reported_success(self):
        if self._reported_success is None:
            self._reported_success = self.all_fits[:, 2, 0].astype(float) > 0.5

        return self._reported_success




if __name__ == '__main__':
    #count_epic_failures()
    results = FitResults('fits_results.txt')
    #results.count_epic_failures()
    #results.count_false_postives()
    print(results.sfit[0, :])
    index = results.sfit[:, 1].astype(float) > 0.5
    results.all_fits = results.all_fits[index, :, :]
    results.tabulate_dchi2()
    #results.get_nev_statistics()
    #results.plot_dchi2()
    #plt.show()