import numpy.testing
from compare import FitResults
from fit_event import FitComparer

results = FitResults('fits_results.txt')
for i, event in enumerate(results.all_fits):
    if (event[1, 0].astype(float) < 0.5) and (event[1, 3].astype(float) < 0.5):
        print(event[0, 0].strip())
        best = event[:, results.best_index[i]]
        method = results.fit_type[results.best_index[i]]
        outstr = ''
        outstr += '{0:15}'.format(method)
        outstr += ' {0:>6}'.format(best[1])
        outstr += ' {0:>8}'.format(best[2])
        outstr += ' {0:>10}'.format(best[3])
        outstr += ' {0:>10} {1:>9} {2:>8}'.format(best[4], best[5], best[6])
        outstr += ' {0:>8}'.format(best[7])
        print(outstr)

        fitter = FitComparer(event=event[0, 0].strip())
        fitter.initial_params = {'t_0': event[4, 3].astype(float),
                                 'u_0': event[5, 3].astype(float),
                                 't_E': event[6, 3].astype(float)}
        try:
            fitter.do_sfit_fit()
            print(fitter.format_results('SFit', fitter.sfit_results).strip())
            #print(fitter.sfit_results.fun,
            #      event[3, results.best_index[i]].astype(float))
            # try:
            #     numpy.testing.assert_almost_equal(
            #         fitter.sfit_results.fun,
            #         event[3, results.best_index[i]].astype(float), decimal=2)
            #     print('Success')
            # except Exception:
            #     print('Failure')
        except Exception as e:
            print('{0:15} Error: {1}'.format('SFit', e))
            #print('Failure')
