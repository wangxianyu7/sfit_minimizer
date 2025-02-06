import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatterSciNotation

from plotting import default_settings


class FitResults():
    _dchi2s = [0.1, 1., 10., 100.]

    def __init__(self, filename='fits_results_formatted.txt'):
        self.fits = pd.read_csv(filename, delim_whitespace=True)
        self.fits['Dchi2'] = self.fits['Chi2'] - self.fits['BestChi2']
        self.n_events = len(np.unique(self.fits['Event']))
        self.fit_types = np.unique(self.fits['Method'])

    def classify(self, fit_type, dchi2, index=None):
        if index is None:
            fit_data = self.fits[self.fits['Method'] == fit_type]
        else:
            fits = self.fits.iloc[index]
            fit_data = fits[fits['Method'] == fit_type]

        reported_success = fit_data['Success'] == 1
        true_success = fit_data['Dchi2'] <= dchi2

        true_positives = np.sum(reported_success & true_success)
        false_positives = np.sum(reported_success & (~true_success))
        true_negatives = np.sum((~reported_success) & (~true_success))
        false_negatives = np.sum((~reported_success) & (true_success))

        return (true_positives, false_positives, true_negatives, false_negatives)

    def output_classifications(self, dchi2=0.1):
        print('Delta chi2 limit: ', dchi2)
        for fit in self.fit_types:
            (true_positives, false_positives, true_negatives, false_negatives) = self.classify(fit, dchi2)

            print(fit)
            print('True Positives:  {0:4} {1:3.0f}'.format(true_positives, 100 * true_positives / self.n_events))
            print('False Positives: {0:4} {1:3.0f}'.format(false_positives, 100 * false_positives / self.n_events))
            print('True Negatives:  {0:4} {1:3.0f}'.format(true_negatives, 100 * true_negatives / self.n_events))
            print('False Negatives: {0:4} {1:3.0f}'.format(false_negatives, 100 * false_negatives / self.n_events))
            print('')

    def make_table_1(self, dchi2s=[0.1, 1., 10., 100.], outfile='table_1.tex'):

        header = '\\tablehead{\n'
        header += '\\colhead{} & \\colhead{} & \\multicolumn{4}{c}{$\\Delta\\chi^2 <$}\\\\\n'
        header += '\\colhead{Algorithm} & \\colhead{Total} '
        header += ''.join([''.join(['& \\multicolumn{2}{|c}{', '{0}'.format(dchi2), '} ']) for dchi2 in dchi2s]) #& \\multicolumn{2}{|c}{1.0}& \\multicolumn{|2}{c}{10.0} & \\multicolumn{2}{|c}{100.0}\\\\\n'
        header += '\\\\'
        header += '\\colhead{} & \\colhead{} '
        header += ''.join(['& \\multicolumn{1}{|c}{N} & \\colhead{\%} ' for i in range(len(dchi2s))])
        #\\colhead{N} & \\colhead{\%} & \\colhead{N} & \\colhead{\%} & \\colhead{N} & \\colhead{\%} & \\colhead{N} & \\colhead{\%} \n'
        header += '\n}\n'
        header += '\\startdata'
        print(header)

        def make_section(index=None):
            if index is None:
                fit_types = self.fit_types
                fits = self.fits
            else:
                fit_types = np.setdiff1d(self.fit_types, ['SFit'])
                #print('fit types', fit_types)
                fits = self.fits.iloc[index]

            print('\\multicolumn{10}{l}{Algorithm Reported Success:}\\\\')
            for fit in fit_types:
                n_success = np.sum((fits['Method'] == fit) & (fits['Success'] == 1))
                fit_str = '\\{0:20} & {1:4} '.format(fit.lower().replace('-', ''), n_success)

                for dchi2 in dchi2s:
                    (true_positives, false_positives, true_negatives, false_negatives) = self.classify(fit, dchi2, index=index)
                    fit_str += '& {0:4} & {1:3.0f} '.format(true_positives, 100 * true_positives / n_success)

                fit_str += '\\\\'
                print(fit_str)

            print('\\hline')
            print('\\multicolumn{10}{l}{Algorithm Reported Failure:}\\\\')
            for fit in fit_types:
                n_failure = np.sum((fits['Method'] == fit) & (fits['Success'] == 0))
                fit_str = '\\{0:20} & {1:4} '.format(fit.lower().replace('-', ''), n_failure)

                for dchi2 in dchi2s:
                    (true_positives, false_positives, true_negatives, false_negatives) = self.classify(fit, dchi2, index=index)
                    fit_str += '& {0:4} & {1:3.0f} '.format(false_negatives, 100 * false_negatives / n_failure)

                fit_str += '\\\\'
                print(fit_str)

        section_1_header = '\\hline\\hline\n'
        section_1_header += ''.join(['\\multicolumn{10}{l}{All', ' {0:4} '.format(int(self.fits.shape[0] / len(self.fit_types))), 'Events:}\\\\\n'])
        section_1_header += '\\hline\\hline'
        print(section_1_header)
        make_section()

        sfit_success = (self.fits['Success'][self.fits['Method'] == 'SFit'] == 1)
        sfit_success_index = np.where(np.repeat(sfit_success, len(self.fit_types)))
        section_2_header = '\\hline\\hline\n'
        section_2_header += ''.join(
            ['\\multicolumn{10}{l}{', '{0}'.format(np.sum(sfit_success)),
             ' Events for which \\sfit\\, reported success:}\\\\\n'])
        section_2_header += '\\hline\\hline'
        print(section_2_header)
        make_section(sfit_success_index)

        sfit_failure = (self.fits['Success'][self.fits['Method'] == 'SFit'] == 0)
        sfit_failure_index = np.where(np.repeat(sfit_failure, len(self.fit_types)))
        section_3_header = '\\hline\\hline\n'
        section_3_header += ''.join(
            ['\\multicolumn{10}{l}{', '{0}'.format(np.sum(sfit_failure)),
             ' Events for which \\sfit\\, reported failure:}\\\\\n'])
        section_3_header += '\\hline\\hline'
        print(section_3_header)
        make_section(sfit_failure_index)

    def make_table_2(self):

        for fit in self.fit_types:
            fit_str = '\\{0:20} '.format(fit.lower().replace('-', ''))
            n_fev = self.fits['nfev'][self.fits['Method'] == fit]
            fit_str += '& {0:5.1f} '.format(np.nanmean(n_fev))
            fit_str += '& {0:3} '.format(int(np.nanmedian(n_fev)))
            fit_str += '& {0:5.1f} '.format(np.nanstd(n_fev))
            fit_str += '& {0:3} '.format(int(np.nanmax(n_fev)))
            fit_str += '\\\\'
            print(fit_str)


    def make_dchi2_plot(self):
        n_zoom = 2
        fits = ['BFGS', 'Neldar-Mead', 'Newton-CG']

        fig = plt.figure(figsize=(10, 6))
        #plt.suptitle(r'$\Delta\chi^2$')
        #fig, axs = plt.subplots
        gs = fig.add_gridspec(ncols=len(fits), nrows=n_zoom)
        axs = gs.subplots()
        print(axs.shape)
        for i, fit_type in enumerate(fits):
            for j in range(n_zoom):
                print(i, j, j * n_zoom + i+1)
                #plt.subplot(n_zoom, 3, j * len(fits) + i+1)
                #ax = fig.add_subplot(gs[n_zoom - j - 1, i])
                plt.sca(axs[n_zoom - j - 1, i])
                plt.gca().set_aspect('equal')
                plt.plot([0, 10**12], [0, 10**12], color='black', linestyle='--')
                plt.scatter(
                    self.fits['Dchi2'][self.fits['Method'] == 'SFit'],
                    self.fits['Dchi2'][self.fits['Method'] == fit_type],
                    marker='o', s=20, zorder=5,
                    c=-self.fits['Success'][self.fits['Method'] == fit_type], cmap='rainbow')

                #plt.xlabel('SFit')
                plt.ylabel(r'$\Delta\chi^2$ {0}'.format(fit_type))
                if j == 0:
                    plt.xlabel(r'$\Delta\chi^2$ SFit')
                if j == n_zoom - 1:
                    #plt.title(fit_type)
                    plt.xscale('log')
                    plt.xlim(10, 10.**5)
                    plt.gca().set_xticks(10.**np.arange(1, 6))
                    plt.gca().get_xaxis().set_major_formatter(
                        LogFormatterSciNotation(base=10, labelOnlyBase=True))
                    plt.yscale('log')
                    plt.ylim(10, 10.**5)
                    plt.gca().set_yticks(10.**np.arange(1, 6))
                    plt.gca().get_yaxis().set_major_formatter(
                        LogFormatterSciNotation(base=10, labelOnlyBase=True))
                else:
                    plt.xlim(-1, 10.)
                    plt.ylim(-1, 10.)

                plt.minorticks_on()

        plt.tight_layout()
        plt.savefig('Dchi2.png', dpi=300)
        plt.show()

    def make_dchi2_plots_separate(self):
        fits = ['BFGS', 'Neldar-Mead', 'Newton-CG']

        def plot_points(fit_type):
            plt.plot([0, 10 ** 12], [0, 10 ** 12], color='black', linestyle='--')
            plt.scatter(
                self.fits['Dchi2'][self.fits['Method'] == 'SFit'],
                self.fits['Dchi2'][self.fits['Method'] == fit_type],
                marker='o', s=20, zorder=5,
                c=-self.fits['Success'][self.fits['Method'] == fit_type], cmap='rainbow')

        for fit_type in fits:
            fig = plt.figure(figsize=(8, 8))
            fig.text(0.05, 0.5, r'$\Delta\chi^2$ {0}'.format(fit_type), ha='center', va='center', rotation=90)
            fig.text(0.5, 0.085, r'$\Delta\chi^2$ SFit', ha='center', va='center')

            gs = fig.add_gridspec(ncols=2, nrows=2)
            gs.update(left=0.15, right=0.9, top=0.9, bottom=0.15, wspace=0.075, hspace=0.075)
            axs = gs.subplots()
            for i in range(2):
                for j in range(2):
                    plt.sca(axs[i, j])
                    plot_points(fit_type)

                    if j == 0:
                        plt.xlim(-1, 10.)
                        plt.minorticks_on()
                    else:
                        plt.xscale('log')
                        plt.xlim(10, 10. ** 5)
                        plt.gca().set_xticks(10. ** np.arange(2, 6))
                        plt.gca().get_xaxis().set_major_formatter(
                            LogFormatterSciNotation(base=10, labelOnlyBase=True))
                        plt.gca().tick_params(right=True, labelright=True, left=True, labelleft=False)

                    if i == 1:
                        plt.ylim(-1, 10.)
                        plt.minorticks_on()
                    else:
                        plt.yscale('log')
                        plt.ylim(10, 10. ** 5)
                        plt.gca().set_yticks(10. ** np.arange(1, 6))
                        plt.gca().get_yaxis().set_major_formatter(
                            LogFormatterSciNotation(base=10, labelOnlyBase=True))
                        plt.gca().tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)

            plt.tight_layout()
            plt.savefig('Dchi2_{0}.png'.format(fit_type), dpi=300)

        plt.show()


if __name__ == '__main__':
    fits = FitResults()
    print('Total Events', fits.n_events)
    print(fits.fits.dtypes)
    #fits.make_table_1()
    #fits.make_table_2()
    fits.make_dchi2_plots_separate()
    #fits.output_classifications()
    #fits.output_classifications(dchi2=1.)