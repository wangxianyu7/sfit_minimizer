import glob, os.path
import numpy as np
import scipy.optimize as op
import MulensModel as mm
import sfit_minimizer


data_dir = os.path.join('/data', 'mulensdata', 'jyee', 'KMTdata', '2018')

KMTdata = np.genfromtxt(
    '../KMT2018.dat', dtype=None, encoding='utf-8',
    names=['Name', 'StarID', 'EF', 'AF', 'RA', 'Dec', 't0', 'tE', 'u0'],
    usecols=range(9))

anomalous_high = np.genfromtxt(
    '../2018_high.list', dtype=None, encoding='utf-8',
    names=['Name'], usecols=[0])
anomalous_low = np.genfromtxt(
    '../2018_low.list', dtype=None, encoding='utf-8',
    names=['Name'], usecols=[0])

def check_anomalous(event):
    num = '{0}:'.format(event[4:8])

    if num in anomalous_high['Name']:
        return True
    elif num in anomalous_low['Name']:
        return True

    else:
        return False


def get_longname(name):
    longname = 'KMT-20{0}-BLG-{1}'.format(name[2:4], name[4:8])
    return longname


def set_parameters(theta, event, parameters_to_fit):
    for (i, key) in enumerate(parameters_to_fit):
        if ((key == 't_E') or (key == 'rho')) and (theta[i] < 0.):
            return np.inf
        else:
            setattr(event.model.parameters, key, theta[i])


def chi2_fun(theta, event, parameters_to_fit):
    """
    for given event set attributes from parameters_to_fit (list of
    str) to values from theta list
    """
    set_parameters(theta, event, parameters_to_fit)
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
    set_parameters(theta, event, parameters_to_fit)
    return event.get_chi2_gradient(parameters_to_fit)


class FitComparer(object):

    def __init__(self, event):
        self.name = event
        self.datasets = self.setup_datasets()
        self.event = None

        self.initial_params = None
        self.initial_guess = None

        self.parameters_to_fit = ['t_0', 'u_0', 't_E']
        self.tol = 1e-5

        self.sfit_results = None
        self.neldar_mead_results = None
        self.newton_cg_results = None
        self.bfgs_results = None

    def _get_bad(self, data):
        # Flag bad data
        # Sky
        mean_sky = np.mean(data['sky'])
        sky_limit = mean_sky + 1. * np.std(data['sky'])
        bad_sky =(data['sky'] > sky_limit)

        # Chi2
        #chi2_limit = 20.
        #bad_chi2 = np.where(data['chi2'] > chi2_limit)

        # Seeing
        bad_see = (data['see'] >
                   (np.mean(data['see']) + 3. * np.std(data['see']))
                   )

        # Concatenate
        bad_data = bad_sky | bad_see
        return bad_data

    def setup_datasets(self):
        files = glob.glob(os.path.join(data_dir, self.name, '*I.pysis'))
        #print(data_dir)
        #print(self.name)
        #print(files)
        datasets = []
        for file_ in files:
            #print(file_)
            data = np.genfromtxt(
                file_, dtype=None, encoding='utf-8',
                names=['HJD', 'Delta_flux', 'flux_err', 'mag', 'mag_err',
                       'see', 'sky', 'secz'])
            bad = self._get_bad(data)
            dataset = mm.MulensData(
                [data['HJD'], -data['Delta_flux'], data['flux_err']],
                phot_fmt='flux')
            dataset.bad = bad
            datasets.append(dataset)
            #print('dataset', dataset)

        #print('datasets', datasets)
        return datasets

    def do_fit(self, verbose=False):
        self.set_initial_params(verbose=verbose)
        self.do_sfit_fit(verbose=verbose)
        self.do_neldar_mead_fit(verbose=verbose)
        self.do_newton_cg_fit(verbose=verbose)
        self.do_bfgs_fit(verbose=verbose)

    def set_initial_params(self, verbose=False):
        index = KMTdata['Name'] == get_longname(self.name)
        kmt_fit = {'t_0': float(KMTdata['t0'][index][0]) ,
                   'u_0': float(KMTdata['u0'][index][0]) ,
                   't_E': float(KMTdata['tE'][index][0]) }
        if verbose:
            print('kmt_fit', kmt_fit)

        def get_params(u0):
            params = {}
            params['t_0'] = kmt_fit['t_0']
            params['u_0'] = u0
            params['t_E'] = kmt_fit['t_E'] * kmt_fit['u_0'] / u0
            return params

        test_u0s = [0.01, 0.3, 0.7, 1.0, 1.5]
        chi2s = []
        for u0 in test_u0s:
            params = get_params(u0)
            event = mm.Event(
                datasets=self.datasets, model=mm.Model(params))
            chi2s.append(event.get_chi2())
            if verbose:
                print('u0, chi2, params', u0, event.chi2, params)

        best = np.argmin(chi2s)
        self.initial_params = get_params(test_u0s[best])
        self.initial_guess = [self.initial_params[key] for key in
                              self.parameters_to_fit]

    def initialize_event(self, verbose=False):
        self.event = mm.Event(
            datasets=self.datasets, model=mm.Model(self.initial_params))
        if verbose:
            print(self.event)

    def do_sfit_fit(self, verbose=False):
        self.initialize_event(verbose=verbose)
        self.sfit_results = sfit_minimizer.fit_mulens_event(
            self.event, parameters_to_fit=self.parameters_to_fit,
        tol=self.tol, verbose=verbose)

    def do_neldar_mead_fit(self, verbose=False):
        self.initialize_event(verbose=verbose)
        self.neldar_mead_results = op.minimize(
            chi2_fun, x0=self.initial_guess, method='nelder-mead',
            args=(self.event, self.parameters_to_fit),
        tol=self.tol)

    def do_newton_cg_fit(self, verbose=False):
        self.initialize_event(verbose=verbose)
        self.newton_cg_results = op.minimize(
            chi2_fun, x0=self.initial_guess,
            args=(self.event, self.parameters_to_fit),
            method='Newton-CG', jac=jacobian,
        tol=self.tol)

    def do_bfgs_fit(self, verbose=False):
        self.initialize_event(verbose=verbose)
        self.bfgs_results = op.minimize(
            chi2_fun, x0=self.initial_guess,
            args=(self.event, self.parameters_to_fit),
            method='BFGS', jac=jacobian,
        tol=self.tol)

    def get_results_header(self):
        outstr = ''
        outstr += '{0:15}'.format('method')
        outstr += ' {0:>6}'.format('scss?')
        outstr += ' {0:>8}'.format('nfev')
        outstr += ' {0:>10}'.format('chi2')
        outstr += ' {0:>10} {1:>9} {2:>8}'.format('t_0', 'u_0', 't_E')
        outstr += ' {0:>8}'.format('njev')
        outstr += '\n'
        return outstr

    def get_nev(self, method, results):
        if method == 'SFit':
            njev = 'nan'
            if results.nit is None:
                nfev = 'nan'
            else:
                nfev = '{0}'.format(results.nit)
        else:
            nfev = '{0}'.format(results.nfev)
            # nfev = '{0}+{1}'.format(results.nfev, results.njev)
            if method == 'Neldar-Mead':
                njev = 'nan'
            else:
                njev = results.njev

        return nfev, njev

    def format_results(self, method, results):
        nfev, njev = self.get_nev(method, results)
        outstr = ''
        outstr += '{0:15}'.format(method)
        outstr += ' {0:6}'.format(results.success)
        outstr += ' {0:>8}'.format(nfev)
        outstr += ' {0:10.2f}'.format(results.fun)
        outstr += ' {0:10.5f} {1:9.5f} {2:8.3f}'.format(
            results.x[0], results.x[1], results.x[2])
        outstr += ' {0:>8}'.format(njev)
        outstr += '\n'

        return outstr

    def format_all_results(self):
        outstr = ''
        for method, results in zip(
                ['SFit', 'Neldar-Mead', 'Newton-CG', 'BFGS'],
                [self.sfit_results, self.neldar_mead_results,
                 self.newton_cg_results, self.bfgs_results]):
            #print(method, results)
            if results is not None:
                outstr += self.format_results(method, results)
                if method == 'SFit':
                    outstr += '{0:>42} {1:10.5f} {2:9.5f} {3:8.3f}\n'.format(
                        'sigmas', results.sigmas[0], results.sigmas[1],
                        results.sigmas[2])


        return outstr

    def print_results(self):
        print(self.get_results_header())
        print(self.format_all_results())


if __name__ == '__main__':
    #events = ['KB180021', 'KB180244', 'KB180774', 'KB181112', 'KB181247']
    n_start = 526
    n_stop = 2781
    events = ['KB18{0:04}'.format(i) for i in np.arange(n_start, n_stop)]
    outfile = open('fit_results_{0:04}_{1:04}.txt'.format(n_start, n_stop), 'w')
    for event in events:
        if check_anomalous(event):
            print('anomalous', event)
        else:
            index = KMTdata['Name'] == get_longname(event)
            if (KMTdata['EF'][index][0] == '1') and (KMTdata['t0'][index][0] != '-'):
                print(event)
                outfile.write('{0}\n'.format(event))
                comparison_fit = FitComparer(event)
                comparison_fit.do_fit(verbose=False)
                comparison_fit.print_results()
                outfile.write(comparison_fit.format_all_results())
            else:
                print('skipping', KMTdata[index])

    outfile.close()