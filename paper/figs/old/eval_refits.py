import numpy as np


class ReFitResults():

    def __init__(self, filename):
        self.sfit = []
        self.other = []
        results = open(filename, 'r')
        for line in results.readlines():
            if line[0] == 'K':
                name = line
            else:
                elements = line.split()
                if elements[1] == 'Error:':
                    values = [np.inf for i in range(6)]
                    values = np.hstack((0, values))
                else:
                    values = [float(value) for value in elements[1:]]
                new_entry = np.hstack((name, elements[0], values))
                if elements[0] == 'SFit':
                    self.sfit.append(new_entry)
                else:
                    self.other.append(new_entry)

        results.close()

        self.all_fits = np.dstack((self.sfit, self.other))
        print('Number of events:', self.all_fits.shape[0])
        print('Number of fits:', self.all_fits.shape[2])

        self._dchi2 = None
        self._reported_success = None
        self._true_success = None

    def calc_dchi2s(self):
        chi2s = self.all_fits[:, 4, 0].astype(float)
        self._dchi2 = np.subtract(chi2s, self.all_fits[:, 4, 1].astype(float))

    def count_success(self):
        reported_success = self.reported_success
        true_success = self.true_success
        print('')
        print('True positives', np.sum(reported_success & true_success))
        print('False positives', np.sum(reported_success & (~true_success)))
        print('True negatives', np.sum(
            (~reported_success) & (~true_success)))
        print('False negatives', np.sum(
            (~reported_success) & true_success))

    def print_failures(self):
        print('')
        print('True Negatives:')
        for event in self.all_fits[
            (~self.reported_success) & (~self.true_success), 0, 0]:
            print(event.strip())

        print('')
        print('False Positives:')
        for event in self.all_fits[
            self.reported_success & (~self.true_success), 0, 0]:
            print(event.strip())

        print('')
        print('False Negatives:')
        for event in self.all_fits[
            (~self.reported_success) & (self.true_success), 0, 0]:
            print(event.strip())


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

    @property
    def dchi2(self):
        if self._dchi2 is None:
            self.calc_dchi2s()

        return self._dchi2


if __name__ == '__main__':
    results = ReFitResults('refit_results.txt')
    #for event in results.all_fits:
    #    print(event.transpose())

    results.print_failures()
    results.count_success()