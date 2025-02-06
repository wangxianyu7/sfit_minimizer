from compare import FitResults
import matplotlib.pyplot as plt

results = FitResults('fits_results.txt')

hist_kwargs = {
    't_0': {'bins': 21, 'range': [-0.005, 0.005]},
    'u_0': {'bins': 21, 'range': [-1, 1]},
    't_E': {'bins': 21, 'range': [-5, 5]}
}

parameters = ['t_0', 'u_0', 't_E']
for fit_type in ['SFit', 'BFGS']:
    print(fit_type)
    plt.figure(figsize=(12, 4))
    plt.suptitle(fit_type)
    for i, parameter in enumerate(parameters):
        plt.subplot(1, len(parameters), i+1)
        plt.minorticks_on()
        plt.yscale('log')
        if parameter == 't_0':
            plt.xlabel(r'$\Delta {0}$'.format(parameter))
            results.plot_hist(fit_type=fit_type, parameter=parameter, **hist_kwargs[parameter])
        else:
            plt.xlabel(r'$\Delta {0} / (\mathrm{{best}}\, {0})$'.format(parameter))
            results.plot_hist(fit_type=fit_type, parameter=parameter, frac=True, **hist_kwargs[parameter])

        if i == 1:
            plt.legend(bbox_to_anchor=[0.5, 1.05], ncol=3, loc='lower center', borderaxespad=0.)

    plt.tight_layout()
    plt.subplots_adjust(top=0.75)
    plt.savefig('hist_{0}.png'.format(fit_type), dpi=300)

#plt.savefig('test_hist.png', dpi=300)
plt.show()